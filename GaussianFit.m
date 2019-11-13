%%% perform Gaussian fitting of FRET histogram
stepfinding = false;
save = false;
fn = 'sim_level1_final_publish';
load([fn '.mat']);
n = 2;
if stepfinding
    E = eff_steps; % assume that we just analyzed the dataset
else
    E = vertcat(E{:});
end

n_comp = 4;
fits = cell(n_comp,1);
BIC = zeros(n_comp,1); AIC = BIC;
opt = statset('MaxIter',1E5);
for i = 1:n_comp % consider up to four states
    fits{i} = fitgmdist(E,i,'Options',opt);
    BIC(i) = fits{i}.BIC;
    AIC(i) = fits{i}.AIC;
end
%%
bins = linspace(-0.1,1.1,101);
f = figure; hold on;
[h,bins_edges] = histcounts(E,bins); bins = bins_edges(1:end-1) + min(diff(bins_edges))/2;
bar(bins,h,'EdgeColor','none','BarWidth',1,'FaceColor',[0.5,0.5,0.5]);
if ~stepfinding
    p = pdf(fits{n},bins'); p = numel(E).*p./sum(p);
    plot(bins,p,'LineWidth',2,'Color',[0,0,0]);  
else
    color = lines(3);
    if n == 2
        plot([1/2,1/2],[0,max(h)],'--','LineWidth',2,'Color',color(1,:));
    elseif n==3
        plot([1/3,1/3],[0,max(h)],'--','LineWidth',2,'Color',color(1,:));
        plot([2/3,2/3],[0,max(h)],'--','LineWidth',2,'Color',color(2,:));
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
    plot(ax_res,bins,w_res,'LineWidth',2);
    set(gca,'Box','on','FontSize',20,'LineWidth',2,'Layer','top');
    ylabel('w. res.');
    ax_res.XTickLabel = [];
    ax_res.Position(2) = ax.Position(2)+ax.Position(4);
    ax_res.Position(4) = 0.1;
    axis('tight');

    switch n
        case 2
            text(ax,ax.XLim(1)+0.1,ax.YLim(2)*0.85,sprintf('E_1 = %.2f\nE_2 = %.2f',fits{n}.mu),'FontSize',20);
        case 3
            text(ax,ax.XLim(1)+0.1,ax.YLim(2)*0.8,sprintf('E_1 = %.2f\nE_2 = %.2f\nE_3=%.2f',fits{n}.mu),'FontSize',20);
    end
end
if save
    if ~stepfinding
        print(f,[fn '_E.png'],'-dpng');
    else
        print(f,[fn '_E_sf.png'],'-dpng');
    end
end