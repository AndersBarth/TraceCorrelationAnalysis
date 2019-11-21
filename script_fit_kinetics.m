addpath(genpath(fileparts(matlab.desktop.editor.getActiveFilename)));
% load data
% 1 = LF, 2 = MF, 3 = HF
% order:
% 1x1, 2x2, 3x3, 1x2, 1x3, 2x1, 2x3, 3x1, 3x2
fn = 'sim_level3_final_publish';
linear = false;
stepfinding = true;
use_weights = true;
n_states = 2;
degenerate = 4; % degenerate, i.e. 2 FRET efficiencies, but 4 states. Specifiy as the number of exponentials for empirical model
if linear
    ext = '_lin';
else
    ext = '';
end
if stepfinding
    sf = '_stepfinding';
else
    sf = '';
end
switch n_states
    case 3
        filenames = {[fn ext sf '_LF_x_LF.mcor'],...
                     [fn ext sf '_MF_x_MF.mcor'],...
                     [fn ext sf '_HF_x_HF.mcor'],...
                     [fn ext sf '_LF_x_MF.mcor'],...
                     [fn ext sf '_LF_x_HF.mcor'],...
                     [fn ext sf '_MF_x_LF.mcor'],...
                     [fn ext sf '_MF_x_HF.mcor'],...
                     [fn ext sf '_HF_x_LF.mcor'],...
                     [fn ext sf '_HF_x_MF.mcor']...
                     };
    case 2
        if stepfinding
            filenames = {[fn ext sf '_LF_x_LF.mcor'],...                    
                         [fn ext sf '_HF_x_HF.mcor'],...
                         [fn ext sf '_LF_x_HF.mcor'],...
                         [fn ext sf '_HF_x_LF.mcor'],...
                         };
        else
            %default to color-FCS
            filenames = {[fn ext '_Idd_x_Idd.mcor'],...                    
                         [fn ext '_Ida_x_Ida.mcor'],...
                         [fn ext '_Idd_x_Ida.mcor'],...
                         [fn ext '_Ida_x_Idd.mcor'],...
                         };
        end
end
 filenames = cellfun(@(x) [fn '_results' filesep x],filenames,'UniformOutput',false);
 
 Cor_Times = [];
 Cor_Average = [];
 Cor_SEM = [];
 for i = 1:numel(filenames)
    d = load(filenames{i},'-mat');
    if i == 1 || numel(d.Cor_Times) == size(Cor_Times,1)
        Cor_Times(:,end+1) = d.Cor_Times;
        Cor_Average(:,end+1) = d.Cor_Average;
        Cor_SEM(:,end+1) = d.Cor_SEM;
    else
        Cor_Times(:,end+1) = Cor_Times(:,end);
        Cor_Average(:,end+1) = zeros(size(Cor_Times,1),1);
        Cor_SEM(:,end+1) = ones(size(Cor_Times,1),1);
        % probably shorter due to short trace
        Cor_Average(1:numel(d.Cor_Average),end) = d.Cor_Average;
        Cor_SEM(1:numel(d.Cor_SEM),end) = d.Cor_SEM;
    end
 end
 t = Cor_Times(:,1);
 valid = t < 150;
 t = t(valid);
 Cor_Average = Cor_Average(valid,:);
 Cor_SEM = Cor_SEM(valid,:);
 if ~use_weights
     Cor_SEM = ones(size(Cor_SEM));
 end
 
 discard_first_point = false;
 if discard_first_point
     t = t(2:end);
     Cor_Average = Cor_Average(2:end,:);
     Cor_SEM = Cor_SEM(2:end,:);
 end
 
 switch n_states
     case 3
         % initial rates
         k0 = rand(1,6)/10;
         lb = zeros(size(k0)); ub = 100*ones(size(k0));
         use_offset = false;
         if use_offset
             k0 = [k0, zeros(1,9)];
             lb = [lb (-1)*ones(1,9)];
             ub = [ub ones(1,9)];
         end
     case 2
         % initial rates
         k0 = rand(1,2);
         lb = zeros(size(k0)); ub = 100*ones(1,2);
         use_offset = false;
         if use_offset
             k0 = [k0, zeros(1,4)];
             lb = [lb (-1)*ones(1,4)];
             ub = [ub ones(1,4)];
         end         
 end

switch n_states
     case 3
        fun = @(x,method) FCS_three_state_kinetics(x,t,Cor_Average,Cor_SEM,method);
     case 2
         if stepfinding
             if ~degenerate
                fun = @(x,method) FCS_two_state_kinetics_fFCS(x,t,Cor_Average,Cor_SEM,method);
             else % consider that states are degenerate
                 switch degenerate
                     case 1
                         % try 1 exp fit
                         k0 = ones(1,3);
                         lb = zeros(1,3);
                         ub = inf(1,3);
                         fun = @(x,method) FCS_two_state_kinetics_fFCS_oneexp(x,t,Cor_Average,Cor_SEM,method);
                     case 2
                         % try 2 exp fit
                         k0 = ones(1,7);
                         lb = zeros(1,7);
                         ub = inf(1,7);
                         fun = @(x,method) FCS_two_state_kinetics_fFCS_biexp(x,t,Cor_Average,Cor_SEM,method);
                     case 3
                         % try 3 exp fit
                         k0 = ones(1,11);
                         lb = zeros(1,11);
                         ub = inf(1,11);
                         fun = @(x,method) FCS_two_state_kinetics_fFCS_threeexp(x,t,Cor_Average,Cor_SEM,method);
                     case 4
                         % try 4 exp fit
                         k0 = 0.1*ones(1,15);
                         lb = zeros(1,15);
                         ub = inf(1,15);
                         fun = @(x,method) FCS_two_state_kinetics_fFCS_fourexp(x,t,Cor_Average,Cor_SEM,method);
                 end
             end
         else %default to colorFCS
            fitE = false;
            if fitE
                fun = @(x,method) FCS_two_state_kinetics_colorFCS(x,0,0,t,Cor_Average,Cor_SEM,method);
                k0 = [k0,rand(1,2)]; lb = [lb,0,0]; ub = [ub,1,1];
            else
                % E1 = 0.31; E2 = 0.70;
                % E1 = 0.24; E2 = 0.72;
                E1 = 0.18; E2 = 0.73;
                fun = @(x,method) FCS_two_state_kinetics_colorFCS(x,E1,E2,t,Cor_Average,Cor_SEM,method);
            end
         end
end
method = 2;
switch method
     case 1 % lsqnonlin
        % first fit simplex
        options = optimset('MaxFunEvals',1E7,'MaxIter',1E7);
        k = fminsearchbnd(@(x) fun(x,2),k0,lb,ub,options);
        % then refine
        [k,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@(x) fun(x,method),k0,lb,ub);
        % get error from nlparci
        ci = nlparci(k,residual,'jacobian',jacobian);
        ci = (ci(:,2)-ci(:,1))/2;
    case 2
        options = optimset('MaxFunEvals',1E7,'MaxIter',1E7);
        k = fminsearchbnd(@(x) fun(x,method),k0,lb,ub,options);
        % use lsqnonlin to estimate errors
        options = optimset('MaxIter',1);
        [k,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@(x) fun(x,1),k,lb,ub);
        % get error from nlparci
        ci = nlparci(k,residual,'jacobian',jacobian);   
        ci = (ci(:,2)-ci(:,1))/2;
end

[w_res,C] = fun(k,1);

mcmc = true;
if mcmc % do Markov Chain Monte Carlo
    logpdf = @(x) (-1)*fun(x,2);
    [smpl, accept] = mhsample(k,1E4,'logpdf',logpdf,'proprnd',@(x) normrnd(x,ci'/10),'symmetric',true,'thin',100);
    k_mcmc = mean(smpl,1);
    ci_mcmc = 1.96*std(smpl,1);
    for i = 1:size(smpl,1);
        p(i) = logpdf(smpl(i,:));
    end
    plot_mcmc = true;
    if plot_mcmc
        if ~(stepfinding && degenerate)
            %%% plot the result of the sampling
            figure('Color',[1,1,1]);
            [~,ax] = plotmatrix_kde(smpl);
            %%% add labels
            if n_states == 2
                ax(1,1).YLabel.String = 'k_{12}';
                ax(2,1).YLabel.String = 'k_{21}';
                ax(2,1).XLabel.String = 'k_{12}';
                ax(2,2).XLabel.String = 'k_{21}';
            elseif n_states == 3
                ax(1,1).YLabel.String = 'k_{12}';
                ax(2,1).YLabel.String = 'k_{31}';
                ax(3,1).YLabel.String = 'k_{21}';
                ax(4,1).YLabel.String = 'k_{23}';
                ax(5,1).YLabel.String = 'k_{31}';
                ax(6,1).YLabel.String = 'k_{32}';
                ax(6,1).XLabel.String = 'k_{12}';
                ax(6,2).XLabel.String = 'k_{31}';
                ax(6,3).XLabel.String = 'k_{21}';
                ax(6,4).XLabel.String = 'k_{23}';
                ax(6,5).XLabel.String = 'k_{31}';
                ax(6,6).XLabel.String = 'k_{32}';
            end
            set(get(gcf,'Children'),'LineWidth',1,'FontSize',14,'Color',[1,1,1]);
        else % fitting descriptive model, just show relaxation times
            %%% plot the result of the sampling
            figure('Color',[1,1,1]);
            [~,ax] = plotmatrix_kde(smpl(:,1:degenerate));
            %%% add labels
            if degenerate == 2
                ax(1,1).YLabel.String = '\tau_{1}';
                ax(2,1).YLabel.String = '\tau_{2}';
                ax(2,1).XLabel.String = '\tau_{1}';
                ax(2,2).XLabel.String = '\tau_{2}';
            elseif degenerate == 3
                ax(1,1).YLabel.String = '\tau_{1}';
                ax(2,1).YLabel.String = '\tau_{2}';
                ax(3,1).YLabel.String = '\tau_{3}';
                ax(3,1).XLabel.String = '\tau_{1}';
                ax(3,2).XLabel.String = '\tau_{2}';
                ax(3,3).XLabel.String = '\tau_{3}';
            elseif degenerate == 4
                ax(1,1).YLabel.String = '\tau_{1}';
                ax(2,1).YLabel.String = '\tau_{2}';
                ax(3,1).YLabel.String = '\tau_{3}';
                ax(4,1).YLabel.String = '\tau_{4}';
                ax(4,1).XLabel.String = '\tau_{1}';
                ax(4,2).XLabel.String = '\tau_{2}';
                ax(4,3).XLabel.String = '\tau_{2}';
                ax(4,4).XLabel.String = '\tau_{4}';
            end
            set(get(gcf,'Children'),'LineWidth',1,'FontSize',14,'Color',[1,1,1]);
        end
    end
end

real_solution = false;
if real_solution
    % instead, plot the real solution
    k(1:6) = [0.1, 0.05, 0.25, 0, 0.025, 0.075];
    if numel(k) > 6
        k(7:end) = 0;
    end
    [w_res,C] = FCS_three_state_kinetics(k,t,Cor_Average,Cor_SEM,1);
end

w_res = reshape(w_res,size(C));

colors = hsv(n_states^2);
figure('Color',[1,1,1]);
ax = axes('Color',[1,1,1],'LineWidth',2,'FontSize',20,'Box','on');
hold on;
error_bar_plot = true;
for i = 1:(n_states^2)
    if error_bar_plot
        errorbar(t,Cor_Average(:,i),Cor_SEM(:,i),'.','MarkerSize',10,'Color',colors(i,:));
    else
        semilogx(t,Cor_Average(:,i),'.','MarkerSize',10,'Color',colors(i,:));
    end
    semilogx(t,C(:,i),'LineWidth',2,'Color',colors(i,:));
end
ax.XScale = 'log';


switch n_states
    case 3
        legend(flipud(ax.Children(1:2:end)),{'LFxLF','MFxMF','HFxHF','LFxMF','LFxHF','MFxLF','MFxHF','HFxLF','HFxMF'},'FontSize',12);
    case 2
        legend(flipud(ax.Children(1:2:end)),{'LFxLF','HFxHF','LFxHF','HFxLF'},'FontSize',12);
end
ax.Position(4) = ax.Position(4)*0.8;
ax_res = axes('ColorOrder',jet(9),'Color',[1,1,1],'LineWidth',2,'FontSize',20,'XTickLabel',[],'Box','on',...
    'Units','normalized','Position',[ax.Position(1),ax.Position(2)+ax.Position(4)+0.05,ax.Position(3),0.15]);
hold on;

for i = 1:(n_states^2)
    semilogx(t,w_res(:,i),'LineWidth',2,'Color',colors(i,:));
end
ax_res.XScale = 'log';
linkaxes([ax,ax_res],'x');

ax.XLim(2) = t(end);
if stepfinding
    ax.YLim(1) = -1.1;
end

rates = k;
if n_states == 3
    K = [-(rates(1)+rates(2)), rates(1), rates(2);...
        rates(3),-(rates(3)+rates(4)),rates(4);...
        rates(5),rates(6),-(rates(5)+rates(6))];
    K
elseif n_states == 2
    k
end
% compute BIC
BIC = chi2_to_bic(-0.5*fun(k,2),numel(k),numel(C));
disp(sprintf('The BIC is: %.2f',BIC));
disp(sprintf('Red. Chi2 is: %.2f',sum(w_res(:).^2)./(numel(C)-numel(k))));
%confint = nlparci(k,residual,'jacobian',jacobian);