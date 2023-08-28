function plot_autocorr_level(results)
% plots N first modes for both stars;
M = size(results.XR,2);
Mmodes = 6;
if ~isfield(results,'XR_burnthin')
    error('Burnthin before you plot chains!')
end

vec = {'00','01','0-1','10','-10','11','1-1','-11','-1-1','02','0-2','20'};
ESS = zeros(2*Mmodes,1);
for i = 1:(2*Mmodes)
    XR_star = results.XR_burnthin(:,i);
    
    subplot(6,2,i)
    [acf,lags,~] = autocorr(XR_star,size(XR_star,1)-1);
    ESS(i) = computeESS(acf);
    plot(lags,acf,'k','linewidth',1.5);
    title(['$\hat{\theta}_{ ',vec{i},'}, \mathrm{ESS}=',sprintf('%.0f',ESS(i)),'$'],'interpreter','latex','fontsize',18);
    ylim([-1 1]);
    ax=gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;

end

        
