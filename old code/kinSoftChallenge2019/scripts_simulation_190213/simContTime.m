%% Version 190212


%% Tested With
%R2011b on Windows 10
%R2017b on Ubuntu 16.04.4 LTS


%% How It Works
%discretization:
%for each frame and for each intensity channel
%	get the ratio of how long each state is populated during this frame
%	draw a random intensity from a gaussian distribution for each state
%	sum intensities, weight according to ratios of each state


%% Fields in struct param
%param.dimLabels
%param.numDims
%param.rmLastFrame
%param.configFileName
%param.numStates
%param.transitionRateMatrix
%param.frameRate
%param.bleachRate
%param.minTraceLen
%param.lookupCollapse
%param.numDegStates
%param.mu
%param.covar
%param.perTraceVariability
%param.numTraces


function simContTime()

%% Preamble

rng('shuffle')  %seeds the random number generator based on the current time

param.dimLabels = {'Idd'; 'Ida'; 'Iaa'};	%actually not used anywhere
param.numDims = size(param.dimLabels, 1);
param.rmLastFrame = 1;
verbose = 0;



%% Import

[param, success] = importSettings(param);

tic

if (success ~= 1)
	disp('Import failed')
	return
end


%Cholesky factorization, formula from help on randn
%sqrt(covar) on main diagonal
R = zeros(param.numDims,param.numDims,param.numDegStates);
for i = 1:param.numDegStates
	R(:,:,i) = chol(param.covar(:,:,i));
end



%% Sanity Check import

%check minTraceLen vs bleachRate
if (param.minTraceLen*10 > 1/param.bleachRate)
	disp('The supplied bleach rate is large compared to the minimal trace length.');
end

%check if rows in trans rate matrix sum to zero
testSum = sum(param.transitionRateMatrix, 2);
testSum = round(1e10*testSum)/1e10;
F = find(testSum ~= 0);
if (size(F, 1) > 0)
	g = sprintf('%d ', F);
	fprintf('The supplied transition rate matrix does not sum to zero in row(s): %s\n', g);
	disp('Aborting.');
	return
end

%check that mu_DA + mu_DD = const
testSum = param.mu(1,1,:) + param.mu(1,2,:);
testSum = testSum(:);	%reshape into a single column vector
testSum = testSum - testSum(1);
testSum = sum(testSum, 1);
if (testSum > 0)
	disp('The supplied mean fluorescence intensities for Idd+Ida are not constant for all states');
	disp('Aborting.');
	return
end

%check that per trace variability parameter is within reasonable bounds
trVar = param.perTraceVariability;
if (trVar < 0 || trVar > 5)
	disp('The supplied per trace variablity parameter is not within reasonable bounds.');
	disp('Aborting.');
	return
end



%% Simulation of all traces

cellAllTraces(param.numTraces,1:2) = {[]};
cellDwellTimes(param.numStates,1) = {[]};

for i = 1:param.numTraces
	[currTrace, mStateTime] = getTrace(param, R, verbose);
	
	cellAllTraces(i,1:2) = {currTrace; mStateTime};
	
	for j = 1:(size(mStateTime,1)-1)		%omit last dwell, as it is inf
		currState = mStateTime(j,1);
		cellDwellTimes(currState) = {[cellDwellTimes{currState}; mStateTime(j,4)]};
	end
end



%% Export

folderNameBase = 'sim_';
formatOut = 'yymmdd_HHMMSS';
folderNameDate = datestr(now,formatOut);

folderName = [folderNameBase folderNameDate];
[status,msg] = mkdir(folderName);

if ~status
	disp(msg)
	return
end

saveCDF = pwd;
cd(folderName);

%loop over all traces and save them
fileNameTrace = 'trace_';
fileNameStateTime = 'state_time_';

for i = 1:param.numTraces
	currTrace = cellAllTraces{i,1};
	fileID = fopen([fileNameTrace num2str(i) '.txt'], 'w');
	fprintf(fileID, '%%t (s)\tIdd (a.u.)\tIda (a.u.)\tIaa (a.u.)\tFRET E\n');
	fprintf(fileID, '%.3f\t%+.3e\t%+.3e\t%+.3e\t%+.3f\n', currTrace.');
	fclose(fileID);
	
	currStateTime = cellAllTraces{i,2};
	fileID = fopen([fileNameStateTime num2str(i) '.txt'], 'w');
	fprintf(fileID, '%%state\tdeg. state\tt_start (s)\tt_dwell (s)\n');
	fprintf(fileID, '%d\t%d\t%+.3e\t%+.3e\n', currStateTime.');
	fclose(fileID);
end


%loop over all states and save cellDwellTimes
fileNameDT = 'dwellTimes_state_';
for i = 1:param.numStates
	currDT = cellDwellTimes{i,1};
	fileID = fopen([fileNameDT num2str(i) '.txt'], 'w');
	fprintf(fileID, ['%%dwell times of state ' num2str(i) ' (s)\n']);
	fprintf(fileID, '%.6f\n', currDT.');
	fclose(fileID);
end


%save config file
[~,name,ext] = fileparts(param.configFileName);
copyfile(param.configFileName, ['Copy_of_' name ext]);


%save struct param?


cd(saveCDF);

assignin('base','cellAllTraces',cellAllTraces);
assignin('base','cellDwellTimes',cellDwellTimes);


toc

mStationaryDist = getStationaryDist(param.transitionRateMatrix);
disp(mStationaryDist)

end





%% Simulating One Trace
%dimensions=3 (I^d_d, I^d_a, I^a_a)

%mu x=1, y=dimensions, z=states
%sigma x=y=dimensions, z=states

%intensity x=state, y=dimensions
%mu+sigma*randn


function [trace, mStateTime] = getTrace(param, R, verbose)
	transitionRateMatrix = param.transitionRateMatrix;
	bleachRate = param.bleachRate;
	minTraceLen = param.minTraceLen;
	lookupCollapse = param.lookupCollapse;
	mu = param.mu;
	numDims = param.numDims;
	numDegStates = param.numDegStates;
	frameRate = param.frameRate;
	rmLastFrame = param.rmLastFrame;
	trVar = param.perTraceVariability;
	
	%get continues time trace
	mStateTime = doSimContTime(transitionRateMatrix, bleachRate, minTraceLen, verbose);
	
	%collapse states that share the same observation
	mStateTime = collapseStates(transitionRateMatrix, lookupCollapse, mStateTime);
	
	%discretize
	bleachTime = mStateTime(end,3);
	traceFrameMax = ceil(bleachTime*frameRate);
	if (rmLastFrame)
		traceFrameMax = traceFrameMax -1;
	end
	
	trace = zeros(traceFrameMax, numDims);
	
	for i = 1:traceFrameMax
		%get ratio of time each state is populated during this frame
		ratioState = getStateRatio(i, 1/frameRate, mStateTime, numDegStates, verbose);
		
		%get random intensities for all states in this frame
		Z = zeros(numDegStates,numDims);
		for j = 1:numDegStates
			Z(j,:) = mu(:,:,j) + randn(1,numDims)*R(:,:,j);
		end
		
		%sum intensities according to the ratio of the corresponding state
		ZZ = ratioState' * Z;
		
		%write intensities to trace
		trace(i,2:4) = ZZ;
	end
	
	%multiply intensities by random factor
	if (trVar)
		%generate uniformly distributed random number in the interval (-1,1)
		r = -1 + 2*rand;
		multiplier = trVar^r;
		trace(:,2:4) = trace(:,2:4) * multiplier;
		%fprintf('multiplier %f\n', multiplier);
	end
	
	%calculate timestamp of frame and FRET efficiency
	trace(:,1) = linspace(0,1/frameRate*(traceFrameMax-1),traceFrameMax);
	trace(:,5) = trace(:,3) ./ (trace(:,2)+trace(:,3));
	
	%trace(:,1)	timestamp
	%trace(:,2)	Idd
	%trace(:,3)	Ida
	%trace(:,4)	Iaa
	%trace(:,5)	FRET eff
end





%% Worker Functions

function mStateTime = doSimContTime(mTransitionRate, bleachRate, minTraceLen, verbose)
	%mTransitionRate	NumStates x NumStates; row: from col: to
	%					in unit of 1/time
	%					cols sum to zero
	%bleachRate			in unit of 1/time
	%					for Hsp90 bleachRate=0.03125 for tau=34sec (160frames in 5Hz)
	%					or k=ln(2)/<tr lenght in sec> // k=1/tau
	%verbose			default (0): shut up
	%
	%mStateTime			[state, collapsed state, start time dwell, dwell time]
	
	
	%in matlab matching of arguments to variables is done positionally, so there is no way to skip an argument.
	%pass [] to indicate a skipped argument
	
	if ~exist('verbose','var')
		verbose = 0;
	end
	
	
	%get stationary distribution
	mStationaryDist = getStationaryDist(mTransitionRate);
	
	
	%get random start state according to stationary distribution
	mStationaryDist_cum = cumsum(mStationaryDist);
	
	r = rand;	%rand returns random number in ]0;1[
	
	for i = 1:size(mStationaryDist, 1)
		if (mStationaryDist_cum(i) > r)
			currS = i;
			if (verbose)
				fprintf('%s : r=%f, currS=%d\n',getBaseName(),r,currS);
			end
			break
		end
	end
	
	
	%get trace length from exponential distribution
	while 1
		r = rand;
		traceLength = -1/bleachRate*log(1-r);	%log(X) returns the natural logarithm ln(x) of each element in array X
		if (traceLength > minTraceLen)
			break
		end
	end
	
	if (verbose)
		fprintf('%s : traceLength=%f\n',getBaseName(),traceLength);
	end
	
	j = 1;
	mStateTime(1,1) = currS;
	mStateTime(1,3) = 0;
	
	while 1
		shortestDwell = +inf;
		%get dwell time for each depopulating transition
		for i =1:size(mStationaryDist, 1)
			if (i == currS)		%skip state that is currently populated
				continue
			end
			transRate = mTransitionRate(currS,i);
			r = rand;
			dwell = -1/transRate*log(1-r);	%log(X) returns the natural logarithm ln(x) of each element in array X
			
			if (dwell<shortestDwell)		%shortest dwell wins
				GoToState = i;
				shortestDwell = dwell;
			end
		end
		
		mStateTime(j,4) = shortestDwell;
		mStateTime(j+1,1) = GoToState;
		mStateTime(j+1,3) = mStateTime(j,3) + shortestDwell;
		
		if (verbose)
			fprintf('%s : GoToState=%d, shortestDwell=%f\n',getBaseName(),GoToState,shortestDwell);
		end
		
		if (mStateTime(j+1,3) > traceLength)	%bleached
			%fprintf('j=%d, traceLength=%f, shortestDwell=%f\n', j, traceLength, shortestDwell);
			mStateTime(j,4) = traceLength - mStateTime(j,3);
			mStateTime(j+1,1) = 0;
			mStateTime(j+1,3) = traceLength;
			mStateTime(j+1,4) = +inf;
			break
		else
			j = j + 1;
			currS = GoToState;
		end
	end
	
	if (verbose)
		%hFig = figure;
		stairs(mStateTime(:,3), mStateTime(:,1:2));
	end
end



%returns the stationary distribution of a transition RATE matrix
%i.e. elements of matrix are rates and rows of the matrix sum to zero
function mStationaryDist = getStationaryDist(mTransitionRate)
	%not sure why the matrix has to be transformed here, but otherwise
	%eigenvalues and eigenvectors are wrong
	[V,D] = eig(mTransitionRate');
	vEigenvalues = diag(D);
	
	%decomposition using floating-point computations, thus result is close to, but not exactly, 0
	vEigenvaluesReal = round(1e10*real(vEigenvalues))/1e10;
	
	%test for == 0 since transition RATE matrix, if transition matrix, this would be == 1
	I = find(vEigenvaluesReal == 0);
	
	if (size(I,1) == 0)
		msg = [getBaseName() ' : Found no stationary distribution. Check input transition rate matrix.'];
		disp(msg);
		mStationaryDist = NaN;
		return
	elseif (size(I,1) > 1)
		msg = [getBaseName() ' : Found more than one stationary distribution. Check input transition rate matrix.'];
		disp(msg);
		mStationaryDist = NaN;
		return
	end
	
	%check if imag(vEigenvalues) == 0
	if (isreal(vEigenvalues(I(1))) == 0)
		msg = [getBaseName() ' : Found no stationary distribution. Check input transition rate matrix.'];
		disp(msg);
		mStationaryDist = NaN;
		return
	end
	
	mStationaryDist = V(:,I(1));
	%normalize to sum=1
	s = sum(mStationaryDist);
	mStationaryDist = mStationaryDist/s;
end



function mStateTime = collapseStates(mTransitionRate, lookupCollapse, mStateTime)
	if (size(lookupCollapse,1) ~= size(mTransitionRate,1))
		msg = [getBaseName() ' : Dimension mismatch.'];
		disp(msg);
	end
	
	
	mState = mStateTime(1:end-1,1);		%last entry is 0
	mStateClps = lookupCollapse(mState(:));
	
	mStateTime(1:end-1,2) = mStateClps;
	mStateTime(end,2) = 0;
end



function ratioState = getStateRatio(frame, frameTime, mStateTime, numDegStates, verbose)
	if (frame <= 0)
		fprintf('%s : frame=%d has to be > 0\n',getBaseName(),frame);
	end
	
	if ~exist('verbose','var')
		verbose = 0;
	end
	
	%mStartEnd: [Obs, t_start, t_end]
	mStartEnd(:,1:2) = mStateTime(:,2:3);
	mStartEnd(:,3) = mStateTime(:,3) + mStateTime(:,4);
	
	start = (frame-1)*frameTime;
	stop = start + frameTime;
	
	k = find(mStartEnd(:,2)>start,1);	%empty if start > traceLength
	j = k-1;
	
	if isempty(k)
		disp('start not found');
	end
	
	ratioState = zeros(numDegStates,1);
	
	while (1)
		if (j == size(mStateTime,1)-1)	%last row is only bleached
			currStart = max([start, mStartEnd(j,2)]);
			currStop = min([stop, mStartEnd(j,3)]);
			ratioState(mStartEnd(j,1)) = ratioState(mStartEnd(j,1)) + currStop - currStart;
			break
		elseif (stop<mStartEnd(j,3))
			currStart = max([start, mStartEnd(j,2)]);
			ratioState(mStartEnd(j,1)) = ratioState(mStartEnd(j,1)) + stop - currStart;
			break
		else
			currStart = max([start, mStartEnd(j,2)]);
			%this is just the dwell time in most cases, but not sure if it
			%would be easier to use the dwell time and have another if in
			%case the currStart is not the start of the dwell
			ratioState(mStartEnd(j,1)) = ratioState(mStartEnd(j,1)) + mStartEnd(j,3) - currStart;
			j = j + 1;
		end
	end
	
	if (verbose)
		fprintf('%s : frameTime=%d, ratioState (not normalized)\n',getBaseName(),frameTime);
		disp(ratioState);
	end
	
	ratioState = ratioState / frameTime;
end





%% Importer

function [param, success] = importSettings(param)
	%param			struct with field numDims
	%				(number of recorded fluorescence channels)
	%				the current function expects numDims==3 with
	%				dimLabels = {'Idd'; 'Ida'; 'Iaa'}
	
	success = 1;
	numDims = param.numDims;
	
	[file,path] = uigetfile('*','Select a simulation config...');
	
	if isequal(file,0)		%user selected Cancel
		disp('User pressed cancel');
		success = -1;
		return
	else
		fileName = fullfile(path,file);
		param.configFileName = fileName;
	end
	
	strFileSetting = fileread(fileName);
	TextAsCells = regexp(strFileSetting, '\n', 'split');
	
	idx = find(strcmp(TextAsCells, '%transition rate matrix'), 1, 'first');
	
	A = importdata(fileName,'\t',idx);
	if (~isfield(A, 'data'))
		A = importdata(fileName,'\t',idx+1);
		if (~isfield(A, 'data'))
			disp('No numeric data found after header ''%transition rate matrix''')
			success = -1;
			return
		end
	end
	
	numStates = size(A.data,2);
	param.numStates = numStates;
	
	param.transitionRateMatrix = A.data(1:numStates,1:numStates);
	
	
	
	Key   = '%frameRate';
	Index = strfind(strFileSetting, Key);
	if (isempty(Index))
		disp('Could not find header ''%frameRate''')
		success = -1;
		return
	else
		Value = sscanf(strFileSetting(Index(1) + length(Key):end), '%f');
		param.frameRate = Value;
	end
	
	
	
	Key   = '%bleachRate';
	Index = strfind(strFileSetting, Key);
	if (isempty(Index))
		disp('Could not find header ''%bleachRate''')
		success = -1;
		return
	else
		Value = sscanf(strFileSetting(Index(1) + length(Key):end), '%f');
		param.bleachRate = Value;
	end
	
	
	
	Key   = '%minimum trace length';
	Index = strfind(strFileSetting, Key);
	if (isempty(Index))
		disp('Could not find header ''%minimum trace length''')
		success = -1;
		return
	else
		Value = sscanf(strFileSetting(Index(1) + length(Key):end), '%f');
		param.minTraceLen = Value;
	end
	
	
	
	Key   = '%lookupCollapse';
	Index = strfind(strFileSetting, Key);
	if (isempty(Index))
		disp('Could not find header ''%lookupCollapse''')
		success = -1;
		return
	else
		Value = sscanf(strFileSetting(Index(1) + length(Key):end), '%f');
		lookupCollapse = Value;
	end
	
	if (numStates ~= size(lookupCollapse, 1))
		disp('lookupCollapse does not match the number of states given in the transition rate matrix.')
		success = -1;
		return
	else
		param.lookupCollapse = lookupCollapse;
		numDegStates = size(unique(lookupCollapse),1);	%is used further down
		param.numDegStates = numDegStates;
	end
	
	
	
	Key   = '%fluorescence intensity, mean';
	Index = find(strcmp(TextAsCells, Key), 1, 'first');
	if (isempty(Index))
		disp('Could not find header ''%fluorescence intensity, mean''')
		success = -1;
		return
	end
	
	param.mu = zeros(1,numDims,numDegStates);	%not really needed, for readability, if removed, change "mu(:,:,1)=" to "mu="
	for i = 1:numDegStates
		Value = sscanf(TextAsCells{1,Index+i}, '%f', [1 3]);
		param.mu(:,:,i) = Value;
	end
	
	
	
	Key   = '%fluorescence intensity, covariance matrix';
	Index = find(strcmp(TextAsCells, Key), 1, 'first');
	if (isempty(Index))
		disp('Could not find header ''%fluorescence intensity, covariance matrix''')
		success = -1;
		return
	end
	
	param.covar = zeros(numDims,numDims,numDegStates);	%not really needed, for readability, if removed, change "mu(:,:,1)=" to "mu="
	for i = 1:numDegStates
		curr = Index + (i-1)*(numDims+1);
		for j = 1:numDims
			Value = sscanf(TextAsCells{1,curr+j}, '%f', [1 3]);
			param.covar(j,:,i) = Value;
		end
	end
	
	
	
	Key   = '%per trace variability';
	Index = strfind(strFileSetting, Key);
	if (isempty(Index))
		disp('Could not find header ''%per trace variability''')
		success = -1;
		return
	else
		Value = sscanf(strFileSetting(Index(1) + length(Key):end), '%f');
		param.perTraceVariability = Value;
	end
	
	
	
	Key   = '%number of traces';
	Index = strfind(strFileSetting, Key);
	if (isempty(Index))
		disp('Could not find header ''%number of traces''')
		success = -1;
		return
	else
		Value = sscanf(strFileSetting(Index(1) + length(Key):end), '%f');
		param.numTraces = Value;
	end

end





%% Little Helper Functions

function basename = getBaseName()
	stack = dbstack;
	basename = stack(2).name;
end





%% Tests for differents parts of the code

function testMatlabRand()
	r = rand(100);	%100x100 matrix of uniformly distributed random numbers
	r_range = [min(r(:)) max(r(:))]
end



function testGetStationaryDist(mTransitionRate)
	len = 1e7;	%1e10 to get close to input
	verbose = 0;
	mStationaryDist = getStationaryDist(mTransitionRate);
	
	occState = zeros(size(mStationaryDist,1),1);
	
	%get random start state according to stationary distribution
	mStationaryDist_cum = cumsum(mStationaryDist);
	
	for i = 1:len
	
		r = rand;	%rand returns random number in ]0;1[
		
		for i = 1:size(mStationaryDist, 1)
			if (mStationaryDist_cum(i) > r)
				currS = i;
				if (verbose)
					fprintf('%s : r=%f, currS=%d\n',getBaseName(),r,currS);
				end
				break
			end
		end
		
		occState(currS) = occState(currS) + 1;

	end
	
	mStationaryDist
	occState = occState / len
	
end



function [counts,centers] = testDistTraceLength()
	len = 1e5;
	traceLen = zeros(len, 1);
	bleachRate = 0.02;	%in unit of 1/time, roughly 0.02 for Hsp90
	
	
	for i = 1:len
		%get trace length from exponential distribution
		r = rand;
		traceLength = -1/bleachRate*log(1-r);	%log(X) returns the natural logarithm ln(x) of each element in array X
		
		traceLen(i) = traceLength;
	end
	
	[counts,centers] = hist(traceLen);
	bar(centers,counts);
	
	%[counts,centers] = testDistTraceLength()
end



function [z, counts,centers] = testMultivariatGauss()

	%from help on randn
	%but watch out, sigma should actually be covariance matrix
	%mu = [1 2];
	%sigma = [1 0.5; 0.5 2];
	%R = chol(sigma);
	%z = repmat(mu,10,1) + randn(10,2)*R
	
	mu = [2 10];
	%sigma = [1 0.5; 0.5 2];
	sigma = [1 0; 0 2];
	R = chol(sigma);
	z = repmat(mu,1e4,1) + randn(1e4,2)*R;
	
	[counts,centers] = hist(z(:,2),20);
	
	%matlab gaussian function returns c1 that is sqrt(2) too large
	%difference between normal dist and matlab funtion
	
	
	
%[z, counts,centers] = testMultivariatGauss();

%%figure;
%%scatter(mR(:,1), mR(:,2));
%%linkdata on
end


