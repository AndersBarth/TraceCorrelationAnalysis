function [Raw, Step] = AnalyzeTraces_v3_Maria(AnalysisAlgorithm,ExpOrSim)
% This program only extracts steps out of data then saves results
% There are two analysis algorithms "AA"
% if AA is 1 then a mixture of a L1 regularized lasso plus the Kalafut
% Vischer with a Penalty Factor of 5 (or to be specified) algorithm is used
% if AA is 2 then the Salapaka Step Finding Algorithm is used
%
% ExpOrSim refers to whether the data is from an experiment or from
% simulation. The difference is because data is stored as a matrix, whereas
% simulations are a structure array.
% 
% Written by Alvaro Crevenna
% LMU Munich, October 2014


% [FileName, PathName] = uigetfile('*.mat', 'Select Data Files to Analyze', 'MultiSelect', 'on');

[FileName, PathName] = uigetfile('*.mat', 'Select Trace to Analyze');
CurrentFileName = [PathName FileName];

% [Raw,Nc, FileName, PathName] = ShortenTRACE(3);
% PF = ; % Penalty Factor to avoid overfitting
if ExpOrSim == 1
    load(CurrentFileName, 'accc','frames'); % Tracy DATA
    Raw = accc;
    
    Nc = max(size(Raw(frames))); % this is the number of traces to analyse
elseif ExpOrSim == 2
    load(CurrentFileName, 'simZMW'); % Simulated DATA
    Raw = simZMW.noise;
    Nc = simZMW.cycles; % this is the number of traces to analyse
end

AA = AnalysisAlgorithm;
h = waitbar(0,'Please wait...');

if AA == 1
    
    PF = 5;
    lambda = 1;
    
    for i = 1:Nc
        waitbar(i / Nc)
        
        if ExpOrSim == 2
            [x, ~, ~, ~] = l1pwc(Raw{i}, lambda, 0); % simulations
            [N, ~] = kvsteps_al(x,PF);
            Step{i} = N;
        elseif ExpOrSim == 1
            [x, ~, ~, ~] = l1pwc(Raw(i,:), lambda, 0); % Exp DATA
            [N, ~] = kvsteps_al(x,PF);
            Step(i,:) = N;
        end
        %     N = simZMW.monomers{i};
        %       N = simZMW.kv_trace{i};
        
        % KV Step Finding Algorithm Alone
        %       [N, ~] = kvsteps_al(accc(i,:),PF);
        
        % KV step finding after L1 regularization from Max Little
        %     [x, ~, ~, ~] = l1pwc(simZMW.noise{i}, lambda, 0); % simulations
%         [x, ~, ~, ~] = l1pwc(Raw(i,:), lambda, 0); % DATA
%         [N, ~] = kvsteps_al(x,PF);
%         Step(i,:) = N;
    end
end

if AA == 2
    for i = 1:Nc
        waitbar(i / Nc)
        if ExpOrSim == 2
            N = stepfit1_alvaro(Raw{i});
            Step{i} = N;
        elseif ExpOrSim == 1
            N = stepfit1_alvaro(Raw(i,:)); % Salapaka Step Finding Algorithm
            % Aggarwal et al Cell Mol Bioen. 2012; 5(1):14-31 Detection of Steps in Single Molecule Data
            Step(i,:) = N;
        end
        
        
    end
end
close(h)

% save('Gel_S', '-struct','D')

if AA == 1
    cd(PathName)
    uisave({'Raw','Step'},[FileName(1:end-4) '_L1KV']);
elseif AA == 2
    cd(PathName)
    uisave({'Raw','Step'},[FileName(1:end-4) '_Salapaka50']);
end