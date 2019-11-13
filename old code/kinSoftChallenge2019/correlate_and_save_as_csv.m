%% load dataset
folder = 'sim_190212_202530_level2';
load([folder '.mat']);
folder = [folder '_results'];
mkdir(folder);
time_res = 0.1;

linear = true;
%% correlate linear
[cdd,tdd,edd] = tracewise_correlation(time_res,Idd,Idd,linear);
[cda,tda,eda] = tracewise_correlation(time_res,Ida,Ida,linear);
[cddda,tddda,eddda] = tracewise_correlation(time_res,Idd,Ida,linear);
[cdadd,tdadd,edadd] = tracewise_correlation(time_res,Ida,Idd,linear);
% calculate efficiency histogram
eff = vertcat(E{:});
[hE,xE] = histcounts(eff,linspace(0,1,101));
xE = xE(1:end-1)+min(diff(xE))/2;
%save results
csvwrite([folder '/corr_Idd_lin.txt'],[tdd,cdd,edd]);
csvwrite([folder '/corr_Ida_lin.txt'],[tda,cda,eda]);
csvwrite([folder '/corr_IddIda_lin.txt'],[tddda,cddda,eddda]);
csvwrite([folder '/corr_IdaIdd_lin.txt'],[tdadd,cdadd,edadd]);
csvwrite([folder '/FRETeffciencies.txt'],eff);
csvwrite([folder '/FRETefficiency_histogram.txt'],[xE',hE']);

linear = false;
% correlate log
[cdd,tdd,edd] = tracewise_correlation(time_res,Idd,Idd,linear);
[cda,tda,eda] = tracewise_correlation(time_res,Ida,Ida,linear);
[cddda,tddda,eddda] = tracewise_correlation(time_res,Idd,Ida,linear);
[cdadd,tdadd,edadd] = tracewise_correlation(time_res,Ida,Idd,linear);
csvwrite([folder '/corr_Idd_log.txt'],[tdd,cdd,edd]);
csvwrite([folder '/corr_Ida_log.txt'],[tda,cda,eda]);
csvwrite([folder '/corr_IddIda_log.txt'],[tddda,cddda,eddda]);
csvwrite([folder '/corr_IdaIdd_log.txt'],[tdadd,cdadd,edadd]);

%% threshold E trace (quasi filtered FCS)
n_states = 3;
switch n_states
    case 2
        threshold = 0.5;
        a = cellfun(@(x) double(x < threshold), E, 'UniformOutput',false);
        b = cellfun(@(x) double(x >= threshold), E, 'UniformOutput',false);
        linear = false;
        [cdd,tdd,edd] = tracewise_correlation(time_res,a,a,linear);
        [cda,tda,eda] = tracewise_correlation(time_res,b,b,linear);
        [cddda,tddda,eddda] = tracewise_correlation(time_res,a,b,linear);
        [cdadd,tdadd,edadd] = tracewise_correlation(time_res,b,a,linear);
    case 3
        threshold1 = 0.36; threshold2 = 0.66;
        a = cellfun(@(x) double(x < threshold1), E, 'UniformOutput',false);
        b = cellfun(@(x) double(x >= threshold1 & x < threshold2), E, 'UniformOutput',false);
        c = cellfun(@(x) double(x >= threshold2), E, 'UniformOutput',false);
        linear = false;
        [caa,taa,eaa] = tracewise_correlation(time_res,a,a,linear);
        [cbb,tbb,ebb] = tracewise_correlation(time_res,b,b,linear);
        [ccc,tcc,ecc] = tracewise_correlation(time_res,c,c,linear);
        [cab,tab,eab] = tracewise_correlation(time_res,a,b,linear);
        [cba,tba,eba] = tracewise_correlation(time_res,b,a,linear);
        [cac,tac,eac] = tracewise_correlation(time_res,a,c,linear);
        [cca,tca,eca] = tracewise_correlation(time_res,c,a,linear);
        [cbc,tbc,ebc] = tracewise_correlation(time_res,b,c,linear);
        [ccb,tcb,ecb] = tracewise_correlation(time_res,c,b,linear);
end