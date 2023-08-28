function plot_chains_level(results)

% plots N first modes for both stars;
M = size(results.XR,2);
Mmodes = 6;
if ~isfield(results,'XR_burnthin')
    error('Burnthin before you plot chains!')
end

vec = {'00','01','0-1','10','-10','11','1-1','-11','-1-1','02','0-2','20'};
ha = tight_subplot(6,2,[.08 .1],[.05 .05],[.07 .03]);

for i = 1:(2*Mmodes)
    XR_star = results.XR_burnthin(:,i);
    
    %subplot(6,2,i)
    axes(ha(i));
    plot(XR_star,'k','linewidth',1.5);
    title(['$\hat{\theta}_{ ',vec{i},'}$'],'interpreter','latex','fontsize',18);
    ax=gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;

end


        
