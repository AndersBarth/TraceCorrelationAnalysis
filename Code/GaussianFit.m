%%% perform Gaussian fitting of FRET histogram
stepfinding = false;
save_figures = true;
fn = 'expSet3';
load([fn '.mat']);
if stepfinding
    load([fn '_results' filesep fn '_eff_steps.mat']);
end
n = 5;
if stepfinding
    E = eff_steps; % assume that we just analyzed the dataset
else
    E = vertcat(E{:});
end

n_comp = 6;
fits = cell(n_comp,1);
BIC = zeros(n_comp,1); AIC = BIC;
opt = statset('MaxIter',1E5);
for i = 1:n_comp % consider up to four states

    S.mu = linspace(0,1,i)';
    S.Sigma = repmat(0.1^2,1,1,i);
    
    fits{i} = fitgmdist(E,i,'Start',S,'Options',opt);
    BIC(i) = fits{i}.BIC;
    AIC(i) = fits{i}.AIC;
end
%%
bins = linspace(-0.25,1.25,151);
f = figure; hold on;
[h,bins_edges] = histcounts(E,bins); bins = bins_edges(1:end-1) + min(diff(bins_edges))/2;
bar(bins,h,'EdgeColor','none','BarWidth',1,'FaceColor',[0.5,0.5,0.5]);
if ~stepfinding
    p = pdf(fits{n},bins'); p = numel(E).*p./sum(p);
    plot(bins,p,'LineWidth',2,'Color',[0,0,0]);
    %%% add individual populations
    if n > 1
        mu = fits{n}.mu;
        sigma = fits{n}.Sigma;
        amp = fits{n}.ComponentProportion';
        for i = 1:n
            p_ind = pdf(gmdistribution(mu(i),sigma(i)),bins');
            p_ind = numel(E).*amp(i).*p_ind./sum(p_ind);
            plot(bins,p_ind,'--','LineWidth',2,'Color',[0.25,0.25,0.25]);
        end
        sigma = sqrt(squeeze(sigma));
    end
else
    color = lines(3);
    if n == 2
        plot([1/2,1/2],[0,max(h)],'--','LineWidth',2,'Color',color(1,:));
    elseif n==3
        plot([1/3,1/3],[0,max(h)],'--','LineWidth',2,'Color',color(1,:));
        plot([2/3,2/3],[0,max(h)],'--','LineWidth',2,'Color',color(2,:));
    end
    p = pdf(fits{n},bins'); p = numel(E).*p./sum(p);
    plot(bins,p,'LineWidth',2,'Color',[0,0,0]);
    %%% add individual populations
    if n > 1
        mu = fits{n}.mu;
        sigma = fits{n}.Sigma;
        amp = fits{n}.ComponentProportion';
        for i = 1:n
            p_ind = pdf(gmdistribution(mu(i),sigma(i)),bins');
            p_ind = numel(E).*amp(i).*p_ind./sum(p_ind);
            plot(bins,p_ind,'--','LineWidth',2,'Color',[0.25,0.25,0.25]);
        end
        sigma = sqrt(squeeze(sigma));
    end
end
xlabel('FRET efficiency, E');
ylabel('Occurrence');
set(gca,'Box','on','FontSize',20,'LineWidth',2,'Layer','top');
axis('tight');

if ~stepfinding
    ax = gca; ax.Position(4) = .7;
    ax_res = axes('Position',ax.Position,'Box','on','FontSize',20,'LineWidth',2,'Layer','top');


    linkaxes([ax,ax_res],'x');

    w_res = (p-h')./sqrt(h');
    valid = isfinite(w_res);
    chi2 = sum(w_res(valid).^2)./(numel(valid)-3*numel(mu));
    plot(ax_res,bins,w_res,'LineWidth',2,'Color',[0,0,0]);
    set(gca,'Box','on','FontSize',20,'LineWidth',2,'Layer','top');
    ylabel('w. res.');
    ax_res.XTickLabel = [];   
    ax_res.Position(2) = 0.835;%(ax.Position(2)+ax.Position(4));
    ax_res.Position(4) = 0.1;
    axis('tight');
    ax.YTickMode = 'auto'; ax.YTickLabelMode = 'auto';
    switch n
        case 2
            text(ax,ax.XLim(1)+0.1,ax.YLim(2)*0.85,sprintf('E_1 = %.2f\nE_2 = %.2f',fits{n}.mu),'FontSize',20);
        case 3
            text(ax,ax.XLim(1)+0.1,ax.YLim(2)*0.8,sprintf('E_1 = %.2f\nE_2 = %.2f\nE_3 = %.2f',fits{n}.mu),'FontSize',20);
    end
end
if save_figures
    if ~stepfinding
        print(f,[fn '_E.png'],'-dpng');
        savefig(f,[fn '_E.fig']);
    else
        print(f,[fn '_E_sf.png'],'-dpng');
        savefig(f,[fn '_E_sf.fig']);
    end
end