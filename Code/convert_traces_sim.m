function convert_traces_sim()

folder = 'sim_level3_final_publish';
% load all trace_X.txt files
try
    for i = 1:1000
        dat = dlmread([folder filesep sprintf('trace_%d.txt',i)],'\t',1,0);
        time{i,1} = dat(:,1);
        Idd{i,1} = dat(:,2);
        Ida{i,1} = dat(:,3);
        Iaa{i,1} = dat(:,4);
        E{i,1} = dat(:,5);
    end
end
%ssave data
save([folder '.mat'],'time','Idd','Ida','Iaa','E');