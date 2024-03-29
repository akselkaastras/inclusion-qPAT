%% This code is for specific figuremaking given data run on the DTU gbar
%% WILL return errors if these files are not available

%% Make mean plot star-shape (5*1e5)

p{1} = 'Results/StarDG/xr_x0seed_1_noiseseed_1_0.01_0.0175_starDG_0.02_2206.mat';
p{2} = 'Results/StarDG/xr_x0seed_1_noiseseed_1_0.01_0.0175_starDG_0.04_2206.mat';
p{3} = 'Results/StarDG/xr_x0seed_1_noiseseed_1_0.01_0.0175_starDG_0.08_2206.mat';
p{4} = 'Results/StarDG/xr_x0seed_1_noiseseed_1_0.01_0.0175_starDG_0.16_2206.mat';
% Make curves
t = linspace(0,2*pi,1000);
kitecurve = [cos(t)+0.65*cos(2*t)-0.65; 1.5*sin(t)];
r = sqrt(0.8+0.8*(cos(4*t)-1).^2);
cushioncurve = [r.*cos(t); r.*sin(t)]; 
curves = [0.18*kitecurve+[0.4;-0.4]; 0.12*cushioncurve+[-0.4;0.4]];
values = [0.2,0.4,0.1];
%%


res = 1000;
figure(1);
ha = tight_subplot(2,2,[.07 .05],[.03 .04],[.03 .03]);
noise = [2,4,8,16];
for i=1:4
    axes(ha(5-i))
    load(p{i})
    results = burnthin(results,5e5,10000);
    error = plot_from_coef_star_mean(results,results.priorpar,res);
    hold on 
    plot3(curves(1,:),curves(2,:),0*curves(1,:)+2,'r:','linewidth',2)
    plot3(curves(3,:),curves(4,:),0*curves(1,:)+2,'r:','linewidth',2)
    hold off
    set(gcf, 'Position',  [100, 100, 1200, 1000])
    set(gcf,'color','w');
    title(['$',num2str(noise(i)) ,'\% $ relative noise'],'interpreter','latex','fontsize',20)
    ax=gca;
    ax.XAxis.FontSize = 14;
    ax.YAxis.FontSize = 14;
end
set(gcf, 'Position',  [100, 100, 980, 1000])
set(gcf,'color','w');
%%

%%
%export_fig 'Figures/recon_plot/mean_plot_star_red.pdf' -opengl
print('Figures/recon_plot/mean_plot_star_red1.eps','-depsc2')
%% Make mean plot level-set

p{1} = 'Results/Level/xr_x0seed_1_noiseseed_1_0.01_0.0175_level_0.02_2206.mat';
p{2} = 'Results/Level/xr_x0seed_1_noiseseed_1_0.01_0.0175_level_0.04_2206.mat';
p{3} = 'Results/Level/xr_x0seed_1_noiseseed_1_0.01_0.0175_level_0.08_2206.mat';
p{4} = 'Results/Level/xr_x0seed_1_noiseseed_1_0.01_0.0175_level_0.16_2206.mat';


res = 1000;
figure(1);
ha = tight_subplot(2,2,[.07 .05],[.03 .04],[.03 .03]);

for i=1:4
    axes(ha(5-i))
    load(p{i})
    results = burnthin(results,1.2e6,100000);
    error = plot_from_coef_level_mean(results,results.priorpar,res);
    hold on 
    plot3(curves(1,:),curves(2,:),0*curves(1,:)+2,'r:','linewidth',2)
    plot3(curves(3,:),curves(4,:),0*curves(1,:)+2,'r:','linewidth',2)
    hold off
    set(gcf, 'Position',  [100, 100, 1200, 1000])
    set(gcf,'color','w');
    title(['$',num2str(noise(i)) ,'\% $ relative noise'],'interpreter','latex','fontsize',20)
    ax=gca;
    ax.XAxis.FontSize = 14;
    ax.YAxis.FontSize = 14;
end
set(gcf, 'Position',  [100, 100, 980, 1000])
set(gcf,'color','w');
%%
%export_fig 'Figures/recon_plot/mean_plot_level_red.pdf' -opengl
print('Figures/recon_plot/mean_plot_level_red1.eps','-depsc2')
%% Autocorrelation plot for 4 % star (6 first coefficients)
p{2} = 'Results/StarDG/xr_x0seed_1_noiseseed_1_0.01_0.0175_starDG_0.04_2206.mat';
load(p{2});
results = burnthin(results,5e5,1);
figure(1);

plot_autocorr_star(results,7);
set(gcf, 'Position',  [100, 100, 600, 1000])
    set(gcf,'color','w');

%%
export_fig 'Figures/recon_plot/acf_star.pdf' -opengl

%% Autocorrelation plot for 4 % level (6 first coefficients)
p{2} = 'Results/Level/xr_x0seed_1_noiseseed_1_0.01_0.0175_level_0.04_2206.mat';
load(p{2});
results = burnthin(results,1.2e6,1);
plot_autocorr_level(results)
set(gcf, 'Position',  [100, 100, 600, 1000])
set(gcf,'color','w');

%%
export_fig 'Figures/recon_plot/acf_level.pdf' -opengl

%% Star chains
p{2} = 'Results/StarDG/xr_x0seed_1_noiseseed_1_0.01_0.0175_starDG_0.04_2206.mat';
load(p{2});
results = burnthin(results,5e5,1);
plot_chains_star(results,6)
set(gcf, 'Position',  [100, 100, 600, 1000])
set(gcf,'color','w');
%%
%export_fig 'Figures/recon_plot/chains_star.pdf' -opengl
print('Figures/recon_plot/chains_star1.eps','-depsc2')
%% Level chains
p{2} = 'Results/Level/xr_x0seed_1_noiseseed_1_0.01_0.0175_level_0.04_2206.mat';
load(p{2});
results = burnthin(results,1.2e6,1);
plot_chains_level(results)
set(gcf, 'Position',  [100, 100, 600, 1000])
set(gcf,'color','w');
%%
%export_fig 'Figures/recon_plot/chains_level.pdf' -opengl
print('Figures/recon_plot/chains_level1.eps','-depsc2')

%% Now we come to error plot
%% Start with star-shaped and then level-set
noiseseed = [1,2,3,4,5];
N = length(noiseseed);
error_star = zeros(5,N);
eps_star = zeros(1,5);
for i = 1:N
    p{1} = ['Results/StarDG/xr_x0seed_1_noiseseed_',num2str(noiseseed(i)),'_0.01_0.0175_starDG_0.01_2206.mat'];
    p{2} = ['Results/StarDG/xr_x0seed_1_noiseseed_',num2str(noiseseed(i)),'_0.01_0.0175_starDG_0.02_2206.mat'];
    p{3} = ['Results/StarDG/xr_x0seed_1_noiseseed_',num2str(noiseseed(i)),'_0.01_0.0175_starDG_0.04_2206.mat'];
    p{4} = ['Results/StarDG/xr_x0seed_1_noiseseed_',num2str(noiseseed(i)),'_0.01_0.0175_starDG_0.08_2206.mat'];
    p{5} = ['Results/StarDG/xr_x0seed_1_noiseseed_',num2str(noiseseed(i)),'_0.01_0.0175_starDG_0.16_2206.mat'];


    res = 1000;
    figure(1);
    
    for j=1:5
        load(p{j})
        results = burnthin(results,5e5,10000);
        eps_star(j) = results.datapar.epssq_approx;
        error_star(j,i) = computeL2error_star(results,results.priorpar,1000);
        j
    end
end
error_level = zeros(5,N);
eps_level = zeros(1,5);
for i = 1:N
    p{1} = ['/work3/akara/qPAT-level/Results/Level/xr_x0seed_1_noiseseed_',num2str(noiseseed(i)),'_0.01_0.0175_level_0.01_2206.mat'];
    p{2} = ['/work3/akara/qPAT-level/Results/Level/xr_x0seed_1_noiseseed_',num2str(noiseseed(i)),'_0.01_0.0175_level_0.02_2206.mat'];
    p{3} = ['/work3/akara/qPAT-level/Results/Level/xr_x0seed_1_noiseseed_',num2str(noiseseed(i)),'_0.01_0.0175_level_0.04_2206.mat'];
    p{4} = ['/work3/akara/qPAT-level/Results/Level/xr_x0seed_1_noiseseed_',num2str(noiseseed(i)),'_0.01_0.0175_level_0.08_2206.mat'];
    p{5} = ['/work3/akara/qPAT-level/Results/Level/xr_x0seed_1_noiseseed_',num2str(noiseseed(i)),'_0.01_0.0175_level_0.16_2206.mat'];


    res = 1000;
    figure(1);
    
    for j=1:5
        load(p{j})
        results = burnthin(results,1.2e6,10000);
        eps_level(j) = results.datapar.epssq_approx;
        error_level(j,i) = computeL2error_level(results,results.priorpar,1000);
        j
    end
end

%%
cmap = colororder();
eps_star_vec = repmat(eps_star,1,5);
eps_level_vec = repmat(eps_level,1,5);

semilogx(eps_level,mean(error_level,2),'-square','color',cmap(1,:),'MarkerFaceColor',cmap(1,:),'linewidth',1.5)
hold on
scatter(eps_level_vec,error_level(:),'square','MarkerEdgeColor',cmap(1,:))
semilogx(eps_star,mean(error_star,2),'-o','color',cmap(2,:),'MarkerFaceColor',cmap(2,:),'linewidth',1.5)
scatter(eps_star_vec,error_star(:),'MarkerEdgeColor',cmap(2,:))
hleg = legend('show');

legend({'Mean, level set','Error, level set','Mean, star set','Error, star set'},'location','northwest','fontsize',14,'interpreter','latex')
xlabel('$\varepsilon_{\mathrm{app}}$','fontsize',20,'interpreter','latex')
ylabel('$\|E[\gamma|y]-\gamma_0\|_{L^2(D)}$','fontsize',20,'interpreter','latex')
%set(gca,'YTickLabel',0.02:0.1:0.4,'fontsize',15)
yticks([0.02,0.03,0.04,0.05,0.06,0.07])
%set(gca,'fontsize',13)
%yticklabels({'A','B','C','D'})
set(gcf, 'Position',  [100, 100, 600, 450])
%hleg.String = []; % delete the last legend entry of the very last plot
%%
set(gcf,'color','w');
%export_fig 'Figures/recon_plot/error_plot.pdf' -opengl
print('Figures/recon_plot/error_plot1.eps','-depsc2')