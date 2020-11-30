function convert_traces_exp()

folder = 'Traces EG 1ms';
fn_root = '';
% load all trace_X.txt files
time = {};
Idd = {};
Ida = {};
E = {};
try
    for i = 1:1000
        try
            try
                dat = dlmread([folder filesep sprintf([fn_root '%.3d.dat'],i)]);
            catch
                dat = dlmread([folder filesep sprintf([fn_root '%.3d.txt'],i)]);
            end
            time{end+1,1} = dat(:,1)-dat(1,1);
            Idd{end+1,1} = dat(:,2);
            Ida{end+1,1} = dat(:,3);
            E{end+1,1} = dat(:,3)./(dat(:,2)+ dat(:,3));
        end
    end
end
% save data
save([folder '.mat'],'time','Idd','Ida','E');