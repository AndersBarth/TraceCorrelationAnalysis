function est=vitpass(y,estin)
global userargs outputnoise measnoise pass wt noise 
global verbose resolution lowsave highsave b a
% partition data into blocks and run viterbi on each block
% global trackcost
%% windowsize
N=length(y);
imp=zeros(1,1000);
imp(1)=1;
impf=filter(b,a,imp);
bw=find(impf>0.2*max(impf),1,'last');
window=max(10,min(length(y),bw));
%% noiseamp
if(pass>=1)%noiseamp
    [h,w]=freqz(b,a,10000);
    noiseamp=sqrt(sum((abs(h)).^2)/length(h));
    noise=(outputnoise-measnoise)/noiseamp;
end
%% ncores
try
    ncores=matlabpool('size');
catch
    ncores=1;
end%ncores
%% Blocksize

steps=diff(estin);
ind=find(steps);
if(pass==1)
    meandwell=500;
else
    meandwell=max(500,round(mean(diff(ind))));
end
blocksize=max(10*meandwell,ceil(N/ncores/ceil(N/10/meandwell)));
r=floor((N-1)/blocksize);
overlap=max(.1*blocksize,2*meandwell);
est=zeros(1,length(y));
%% noise reestimation
if(pass>1 && isempty(userargs{3}))%outputnoise 
    outputnoise=noise_std(y-filter(b,a,estin)); %re-estimate noise
end
%% envelope
if(pass==1)%envelope
    [low high maxstep minstep miny maxy]=findenvelope2(estin,window);
    lowsave=low;
    highsave=high;
    resolution=mean(high-low)/20; % refine resolution 
else
[low high maxstep minstep miny maxy resolution ]=findenvelope3(estin,window,resolution,lowsave,highsave);
lowsave=low;
highsave=high;
end
low=round((low-miny)/resolution)+1;
high=round((high-miny)/resolution)+1;
message=['Pass:' num2str(pass) ', Resolution=' sprintf('%4.3f',resolution) ];
% if(verbose)
    fprintf(message);
% end
%% compute weight 

if(pass<=1 )
    range=minstep-1:resolution:maxstep;
    [val in]=min(abs(range));
    range=range-range(in);
    wt.zeroindex=in;
    wt.pdf=abs(sign(range));
else
   wt=computeweight(estin,minstep,maxstep);
end%compute weight
start=zeros(1,r+1);
stop=zeros(1,r+1);
eststart=zeros(1,r+1);
eststop=zeros(1,r+1);

%% parallel operation
% if(verbose)
    msg=[char(zeros(1,r+1)+'|')  ' '];
    fprintf(msg);
% end
for j=1:r+1
    i=j-1;
    start(j)=max(blocksize*i+1-overlap,1);
    stop(j)=min(blocksize*(i+1)+overlap,N);
    eststart(j)=blocksize*i+1;
    eststop(j)=min(N,blocksize*(i+1));
    yme{j}=y(start(j):stop(j));
    lowtmp{j}=low(start(j):stop(j));
    hightmp{j}=high(start(j):stop(j));
end %%% slice data for parallel operation
param.noise=noise;
param.measnoise=measnoise;
param.outputnoise=outputnoise;
param.wt=wt;
param.resolution=resolution;
param.b=b;
param.a=a;
param.pass=pass;
parfor j=1:r+1  %*****************parfor**************
    [tmp costtmp]=viterbistepdetector(yme{j},miny,lowtmp{j},hightmp{j},param);
    myest{j}=tmp;
%     if(verbose)
        disp([char(8) char(8)  ]);
%     end
end
disp([char(8)]);
for j=1:r+1
    tempest=myest{j};
    est(eststart(j):eststop(j))=tempest(eststart(j)-start(j)+1:eststop(j)-start(j)+1);
end  %% merge sliced data 

%% display progress
if(verbose)
    figure(4);clf;
    %         nr=floor(sqrt(passes));
    %         nc=ceil(passes/nr);
    subplot(2,1,1);
    v=1:N;
    rv=N:-1:1;
    polyx=[v rv];
    polyenv=[(low-1)*resolution+miny (fliplr(high)-1)*resolution+miny];
    polyy=[y rv*0];
    py=fill(polyx,polyy,'k');hold on;set(py,'facealpha',0,'edgealpha',0.1);
    pf=fill(polyx,polyenv,'r');set(pf,'facealpha',0.2,'edgealpha',0);hold on;
    plot(est,'r','linewidth',2);
    xlabel('Samples');
    ylabel('Position');
    title(['Iteration ' num2str(pass)]);
%     legend('Data','Envelope','Fit','location','best');
    axis tight;
    
    subplot(2,1,2);
    cla
    x=resolution*((1:length(wt.pdf))-wt.zeroindex);
    smoothedhistogram=exp(-wt.pdf);
    smoothedhistogram=smoothedhistogram-min(smoothedhistogram);
    smoothedhistogram(wt.zeroindex)=0;        
    [bins nhist]=stephist(est,linspace(min(x),max(x),25),gcf,'normalized','bar');    
    if(max(nhist))
        bar(bins,nhist*max(smoothedhistogram)/max(nhist),'b');
    end        
    hold on;
    plot(x-resolution/2,smoothedhistogram,'r','linewidth',1);
    xlabel('Step-size');       
    ylabel('Probability');
    axis tight;
    drawnow();
    if(verbose==2)
        str=sprintf('StepFitIterationProgress %d .fig',pass);
        hgsave(str);
    end
end