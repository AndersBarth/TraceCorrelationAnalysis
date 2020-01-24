function convert_traces_exp()

folder = 'expSet3';
fn_root = 'Mov017-2016-07-07.traces_tr';
% load all trace_X.txt files
time = {};
Idd = {};
Ida = {};
E = {};
try
    for i = 1:1000
        try
            dat = dlmread([folder filesep sprintf([fn_root '%.3d.dat'],i)]);
            time{end+1,1} = dat(:,1)-dat(1,1);
            Idd{end+1,1} = dat(:,2);
            Ida{end+1,1} = dat(:,3);
            E{end+1,1} = dat(:,3)./(dat(:,2)+ dat(:,3));
        end
    end
end
% save data
save([folder '.mat'],'time','Idd','Ida','E');