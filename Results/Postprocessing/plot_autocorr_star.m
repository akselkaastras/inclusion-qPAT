function plot_autocorr_star(results,N)
% plots N first modes for both stars;
nstars = 2;
M = size(results.XR,2);
Mmodes = M/nstars;

if ~isfield(results,'XR_burnthin')
    error('Burnthin before you plot chains!')
end

vec = [0,-1,1,-2,2,-3,3,-4,4];
ESS = zeros(N,nstars);

for i = 1:nstars
    XR_star = results.XR_burnthin(:,(i-1)*Mmodes+(1:Mmodes));
    for j = 1:N
        subplot(N,nstars,nstars*(j-1)+i)
        [acf,lags,~] = autocorr(XR_star(:,j),size(XR_star,1)-1);
        Kess = computeESS(acf);
        ESS(j,i) = Kess;
        plot(lags,acf,'k','linewidth',1.5);
        title(['$\hat{\theta}_{ ',num2str(vec(j)),'}, \mathrm{ESS}=',sprintf('%.0f',Kess),'$'],'interpreter','latex','fontsize',18);
        ylim([-1 1]);
            ax=gca;
        ax.XAxis.FontSize = 12;
        ax.YAxis.FontSize = 12;
    end
end

        
