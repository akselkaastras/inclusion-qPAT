function plot_chains_star(results,N)
% plots N first modes for both stars;
nstars = 2;
M = size(results.XR,2);
Mmodes = M/nstars;

if ~isfield(results,'XR_burnthin')
    error('Burnthin before you plot chains!')
end

vec = [0,-1,1,-2,2,-3,3,-4,4];
ha = tight_subplot(N,nstars,[.08 .1],[.05 .05],[.07 .03]);

for i = 1:nstars
    XR_star = results.XR_burnthin(:,(i-1)*Mmodes+(1:Mmodes));
    for j = 1:N
        axes(ha(nstars*(j-1)+i))
        %subplot(N,nstars,nstars*(j-1)+i)
        plot(XR_star(:,j),'k','linewidth',1.5);
        title(['$\hat{\theta}_{ ',num2str(vec(j)),'}$'],'interpreter','latex','fontsize',18);
            ax=gca;
        ax.XAxis.FontSize = 12;
        ax.YAxis.FontSize = 12;
    end
end

        
