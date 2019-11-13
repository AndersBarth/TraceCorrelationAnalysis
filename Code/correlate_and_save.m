%% load dataset
filename = 'sim_level2_final_publish';
n_states = 3;
time_res = 0.1;
step_finding = true;
linear = false;

addpath(genpath('.'));
load([filename '.mat']);
folder = [filename '_results'];
mkdir(folder);

%% correlate
if linear
    filename = [filename '_lin'];
end
[c,t,e,ca] = tracewise_correlation(time_res,Idd,Idd,linear);
save_to_fcsfit('Idd','Idd',[pwd filesep folder],filename,t,c,e,ca);
[c,t,e,ca] = tracewise_correlation(time_res,Ida,Ida,linear);
save_to_fcsfit('Ida','Ida',[pwd filesep folder],filename,t,c,e,ca);
[c,t,e,ca] = tracewise_correlation(time_res,Idd,Ida,linear);
save_to_fcsfit('Idd','Ida',[pwd filesep folder],filename,t,c,e,ca);
[c,t,e,ca] = tracewise_correlation(time_res,Ida,Idd,linear);
save_to_fcsfit('Ida','Idd',[pwd filesep folder],filename,t,c,e,ca);
%% calculate efficiency histogram
eff = vertcat(E{:});
[hE,xE] = histcounts(eff,linspace(0,1,101));
xE = xE(1:end-1)+min(diff(xE))/2;
if step_finding
    %%% find steps
    E_steps = cell(size(E));
    for i = 1:numel(E)
        %disp(sprintf('Processing trace %d\n',i));
        try % stepfinding sometimes failes
            E_steps{i} = stepfit1_alvaro(E{i},'outputnoise',0.05)';
        end
    end
    E_steps = E_steps(~cellfun(@isempty,E_steps));
    %E_steps = cellfun(@(x) stepfit1_alvaro(x)',E,'UniformOutput',false); % run program without any parameters
    eff_steps = vertcat(E_steps{:});
    hE_steps = histcounts(eff_steps,linspace(0,1,101));
end
%% threshold E trace (quasi filtered FCS)
switch n_states
    case 2
        threshold = 0.5;
        lf = cellfun(@(x) double(x < threshold), E, 'UniformOutput',false);
        hf = cellfun(@(x) double(x >= threshold), E, 'UniformOutput',false);
        [c,t,e,ca] = tracewise_correlation(time_res,lf,lf,linear);
        save_to_fcsfit('LF','LF',[pwd filesep folder],filename,t,c,e,ca);
        [c,t,e,ca] = tracewise_correlation(time_res,hf,hf,linear);
        save_to_fcsfit('HF','HF',[pwd filesep folder],filename,t,c,e,ca);
        [c,t,e,ca] = tracewise_correlation(time_res,lf,hf,linear);
        save_to_fcsfit('LF','HF',[pwd filesep folder],filename,t,c,e,ca);
        [c,t,e,ca] = tracewise_correlation(time_res,hf,lf,linear);
        save_to_fcsfit('HF','LF',[pwd filesep folder],filename,t,c,e,ca);
        if step_finding
            threshold = 0.5;
            lf = cellfun(@(x) double(x < threshold), E_steps, 'UniformOutput',false);
            hf = cellfun(@(x) double(x >= threshold), E_steps, 'UniformOutput',false);
            [c,t,e,ca] = tracewise_correlation(time_res,lf,lf,linear);
            save_to_fcsfit('LF','LF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);
            [c,t,e,ca] = tracewise_correlation(time_res,hf,hf,linear);
            save_to_fcsfit('HF','HF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);
            [c,t,e,ca] = tracewise_correlation(time_res,lf,hf,linear);
            save_to_fcsfit('LF','HF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);
            [c,t,e,ca] = tracewise_correlation(time_res,hf,lf,linear);
            save_to_fcsfit('HF','LF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);
        end
    case 3
        threshold1 = 0.3; threshold2 = 0.72;
        lf = cellfun(@(x) double(x < threshold1), E, 'UniformOutput',false);
        mf = cellfun(@(x) double(x >= threshold1 & x < threshold2), E, 'UniformOutput',false);
        hf = cellfun(@(x) double(x >= threshold2), E, 'UniformOutput',false);
        [c,t,e,ca] = tracewise_correlation(time_res,lf,lf,linear);
        save_to_fcsfit('LF','LF',[pwd filesep folder],filename,t,c,e,ca);
        [c,t,e,ca] = tracewise_correlation(time_res,mf,mf,linear);
        save_to_fcsfit('MF','MF',[pwd filesep folder],filename,t,c,e,ca);
        [c,t,e,ca] = tracewise_correlation(time_res,hf,hf,linear);
        save_to_fcsfit('HF','HF',[pwd filesep folder],filename,t,c,e,ca);
        [c,t,e,ca] = tracewise_correlation(time_res,lf,mf,linear);
        save_to_fcsfit('LF','MF',[pwd filesep folder],filename,t,c,e,ca);     
        [c,t,e,ca] = tracewise_correlation(time_res,mf,lf,linear);
        save_to_fcsfit('MF','LF',[pwd filesep folder],filename,t,c,e,ca);
        [c,t,e,ca] = tracewise_correlation(time_res,lf,hf,linear);
        save_to_fcsfit('LF','HF',[pwd filesep folder],filename,t,c,e,ca);
        [c,t,e,ca] = tracewise_correlation(time_res,hf,lf,linear);
        save_to_fcsfit('HF','LF',[pwd filesep folder],filename,t,c,e,ca);
        [c,t,e,ca] = tracewise_correlation(time_res,mf,hf,linear);
        save_to_fcsfit('MF','HF',[pwd filesep folder],filename,t,c,e,ca);
        [c,t,e,ca] = tracewise_correlation(time_res,hf,mf,linear);
        save_to_fcsfit('HF','MF',[pwd filesep folder],filename,t,c,e,ca);
        
        if step_finding
            threshold1 = .3; threshold2 = .72;            
            lf = cellfun(@(x) double(x < threshold1), E_steps, 'UniformOutput',false);
            mf = cellfun(@(x) double(x >= threshold1 & x < threshold2), E_steps, 'UniformOutput',false);
            hf = cellfun(@(x) double(x >= threshold2), E_steps, 'UniformOutput',false);
            [c,t,e,ca] = tracewise_correlation(time_res,lf,lf,linear);
            save_to_fcsfit('LF','LF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);
            [c,t,e,ca] = tracewise_correlation(time_res,mf,mf,linear);
            save_to_fcsfit('MF','MF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);
            [c,t,e,ca] = tracewise_correlation(time_res,hf,hf,linear);
            save_to_fcsfit('HF','HF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);
            [c,t,e,ca] = tracewise_correlation(time_res,lf,mf,linear);
            save_to_fcsfit('LF','MF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);     
            [c,t,e,ca] = tracewise_correlation(time_res,mf,lf,linear);
            save_to_fcsfit('MF','LF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);
            [c,t,e,ca] = tracewise_correlation(time_res,lf,hf,linear);
            save_to_fcsfit('LF','HF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);
            [c,t,e,ca] = tracewise_correlation(time_res,hf,lf,linear);
            save_to_fcsfit('HF','LF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);
            [c,t,e,ca] = tracewise_correlation(time_res,mf,hf,linear);
            save_to_fcsfit('MF','HF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);
            [c,t,e,ca] = tracewise_correlation(time_res,hf,mf,linear);
            save_to_fcsfit('HF','MF',[pwd filesep folder],[filename '_stepfinding'],t,c,e,ca);
        end
end