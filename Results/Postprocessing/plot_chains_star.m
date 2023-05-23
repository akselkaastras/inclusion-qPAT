function plot_chains_star(results,N)
% plots N first modes for both stars;
nstars = 2;
M = size(results.XR,2);
Mmodes = M/nstars;

if ~isfield(results,'XR_burnthin')
    error('Burnthin before you plot chains!')
end



for i = 1:nstars
    XR_star = results.XR_burnthin(:,(i-1)*Mmodes+(1:Mmodes));
    for j = 1:N
        subplot(N,nstars,nstars*(j-1)+i)
        plot(XR_star(:,j),'k','linewidth',1.5);
        title(['Coefficient ',num2str(j)]);
    end
end

        
