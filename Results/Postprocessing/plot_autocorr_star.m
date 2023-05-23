function plot_autocorr_star(results,N)
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
        [acf,lags,~] = autocorr(XR_star(:,j),size(XR_star,1)-1);
        plot(lags,acf,'k','linewidth',1.5);
        title(['Autocorrelation coefficient ',num2str(j)]);
        ylim([-1 1]);
    end
end

        
