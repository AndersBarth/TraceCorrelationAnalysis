function [c,t,e] = tracewise_correlation(time_res,sig1,sig2,linear)
if nargin < 3
    sig2 = sig1;
end
maxtime = cellfun(@numel,sig1);
for i = 1:numel(maxtime)
    time{i} = (1:maxtime(i))';
end
[Cor_Array,Timeaxis] = CrossCorrelation(time,time,maxtime,sig1,sig2,2,linear);
Cor_error = std(Cor_Array,0,2);
Cor = mean(Cor_Array,2);

%% plot and fit
offset = 0;
c = Cor(1:end-offset);
e = Cor_error(1:end-offset);
t = time_res*Timeaxis(1:end-offset);
fitfun = @(f,xdata) (f(1)*exp(-xdata./f(2))+f(3));
cf = lsqcurvefit(@(x,xdata) fitfun(x,xdata)./e,[20*c(2)./e(2),Timeaxis(end)/10,-Inf],t,c./e,[-Inf,0,0],[Inf,Inf,Inf]);
fitres = fitfun(cf,0:0.01:max(t));
figure;hold on;
errorbar(t,c,e);
p = plot(0:0.01:max(t),fitres);
p.LineWidth = 2;
set(gca,'XScale','log');
if cf(1) < 0
    cf(1) = -cf(1);
end
% calculate fitted parameters
k1 = (cf(2)*(1+cf(1)))^(-1);
k2 = (cf(1)*(k1));
tau1 = 1/k1;
tau2 = 1/k2;
disp(sprintf('tau = %f s',cf(2)));
disp(sprintf('k1 = %f s^(-1)',k1));
disp(sprintf('k2 = %f s^(-1)',k2));
disp(sprintf('offset = %f',cf(3)));