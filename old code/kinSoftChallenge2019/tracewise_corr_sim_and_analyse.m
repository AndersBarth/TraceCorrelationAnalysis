%%% simulate traces
Sig = [];
Data = [];
n_traces = 1000;
a1 = 0.2;
a2 = 1-a1;
tau = [20, 100];
I = [300, 100];

for i = 1:n_traces
    dwell = 0;
    s = 0;
    dur(i) = randi([300,1000]);
    %%% simulate 2state kinetics
    s(1) = binornd(1,a2)+1;
    dwell(1) = 1;
    while sum(dwell) < dur(i)
        dwell(end+1) = round(exprnd(tau(s(end))));
        if s(end) == 1
            s(end+1) = 2;
        elseif s(end) == 2
            s(end+1) = 1;
        end
    end
    dwell = cumsum(dwell);
    dwell(end) = dur(i); % truncate last dwell time
    for j = 1:numel(dwell)-1
        Sig{i}(dwell(j):(dwell(j+1)-1),1) = poissrnd(I(s(j)),dwell(j+1)-dwell(j),1);
    end
end

maxtime = cellfun(@numel,Sig);
for i = 1:n_traces
    Data{i} = (1:maxtime(i))';
end

[Cor_Array,Timeaxis] = CrossCorrelation(Data,Data,maxtime,Sig,Sig,2);
Cor_error = std(Cor_Array,0,2);
Cor = mean(Cor_Array,2);
Timeaxis = Timeaxis;

%% plot and fit
c = Cor(1:end-10);
t = Timeaxis(1:end-10);
fitfun = @(a,b,c,x) a*exp(-x./b)+c;
cf = fit(t,c,fitfun,'StartPoint',[c(1),Timeaxis(end)/10,0]);

figure;hold on;
scatter(t,c,'o');
p = plot(cf);
p.LineWidth = 2;
set(gca,'XScale','log');

% calculate fitted parameters
k1 = (cf.b*(1+cf.a))^(-1);
k2 = (cf.a*(k1));
tau1 = 1/k1;
tau2 = 1/k2;