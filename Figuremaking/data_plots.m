hmax = 0.0175;
fine_hmax = 0.01;
meshpar = mesh_comp(hmax);

% Make curves
t = linspace(0,2*pi,1000);
kitecurve = [cos(t)+0.65*cos(2*t)-0.65; 1.5*sin(t)];
r = sqrt(0.8+0.8*(cos(4*t)-1).^2);
cushioncurve = [r.*cos(t); r.*sin(t)]; 
curves = [0.18*kitecurve+[0.4;-0.4]; 0.12*cushioncurve+[-0.4;0.4]];
values = [0.2,0.4,0.1];


noiseed = 1;
noiselevel = 0.02;
%% Read user specified curves
ninclusions = size(curves,1)/2;
ncurve = size(curves,2);

%% Initialize mesh
meshpar_fine = mesh_comp(fine_hmax);

% Find mesh on boundary
npoints = length(meshpar_fine.p);
npoints_coarse = length(meshpar.p);
pN = length(meshpar_fine.p);
pNN = pN-size(meshpar_fine.e(1,:),2);
meshpar.NZ = setdiff(1:length(meshpar.p),meshpar.e(1,:));
meshpar_fine.NZ = setdiff(1:pN,meshpar_fine.e(1,:));

%% Make gamma from curves
Gamma = zeros(ninclusions,meshpar_fine.HN);
for i = 1:ninclusions
    index = 2*(i-1);
    nodes = [curves(index+1,:)',curves(index+2,:)'];
    Gamma(i,:) = Gamma(i,:) + mean(reshape(inpoly2([meshpar_fine.xq, meshpar_fine.yq],nodes),[meshpar_fine.HN meshpar_fine.n]),2)';
end
gamma = values(ninclusions+1)+zeros(1,meshpar_fine.HN);
for i = 1:ninclusions
    gamma = gamma + values(i)*Gamma(i,:);
end


Mua = zeros(ninclusions,length(meshpar_fine.p));
for i = 1:ninclusions
    index = 2*(i-1);
    nodes = [curves(index+1,:)',curves(index+2,:)'];
    Mua(i,:) = Mua(i,:) + inpoly2([meshpar_fine.p(1,:)', meshpar_fine.p(2,:)'],nodes)';
end
mua = values(ninclusions+1)+zeros(1,length(meshpar_fine.p));
for i = 1:ninclusions
    mua = mua + values(i)*Mua(i,:);
end

%% Make smoothened scattering mus with Gaussian smoothing
nimage = 300;
[X,Y] = meshgrid(linspace(-1,1,nimage));
x = X(:);
y = Y(:);
Musimage = zeros(ninclusions, nimage*nimage);

% Make image
for i = 1:ninclusions
    index = 2*(i-1);
    Musimage(i,:) = inpolygon(x,y,curves(index+1,:),curves(index+2,:));
end

% Include factor 100
musimage = 100*values(ninclusions+1)+zeros(1,nimage*nimage);
for i = 1:ninclusions
    musimage = musimage + 100*values(i)*Musimage(i,:);
end
musimage = reshape(musimage,nimage,nimage);

% Gaussian filtering
musimage = imgaussfilt(musimage,15);

% Interpolate onto mesh
mus = interp2(X,Y,musimage,meshpar_fine.p(1,:),meshpar_fine.p(2,:));


%% Make diffusion from q and mus
constantg = 0.8;   % Anisotropy factor of the Heneyey-Greenstein scattering function
constantkind = 2;  % Dimension

murs = mus.*(1-constantg);
D = 1./constantkind * (1./(mua + murs));

%% plot parameters
figure(1);
ha = tight_subplot(1,2,[.05 .05],[.1 .1],[.1 .1]);
axes(ha(1))
pdesurf(meshpar_fine.p,meshpar_fine.t,gamma)
colormap default
view(2)
%title('Absorption coefficient, $\gamma$','interpreter','latex','fontsize',200)
xticks([-1,0,1])
yticks([-1,0,1])
%xlim([0,2*pi])

set(gca,'XTickLabel',{'-1','0','1'},'fontsize',15)
colorbar
cbh = colorbar('XTick', [0.1+eps,0.3,0.5],'fontsize',11);

box on
axis square
ax=gca;
ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 11;
%subplot(132)
%plot_from_gamma(full(mus),meshpar_fine)
%title('Scattering coefficient')
%colorbar
axes(ha(2))
plot_from_gamma(full(D),meshpar_fine)
%hT = title('Diffusion coefficient, a','interpreter','latex','fontsize',200);
xticks([-1,0,1])
yticks([-1,0,1])
%xlim([0,2*pi])
set(gca,'XTickLabel',{'-1','0','1'},'fontsize',15)
grid off
colorbar
set(gcf, 'Position',  [100, 100, 1000, 360])
set(gcf,'color','w');
box on
axis square
ax=gca;
ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 11;
cbh = colorbar('XTick', 0.01*(6:4:22),'fontsize',11);
%%


%% Export figure
%export_fig 'Figures/data_plot/coefficients1.eps'
print('Figures/data_plot/coefficients1.eps','-depsc2')
%% evaluate forward model
% define source from wfun
sigma = 0.5;
m1 = 0.5*[sqrt(2),sqrt(2)];
m2 = 0.5*[-sqrt(2),sqrt(2)];
scale1 = 10;
scale2 = 5;
scale3 = 2;
figure(2);

% Source is Gaussian
wfunfun = @(x1,x2,m,scale) scale*exp(-1/(2*sigma^2)*((x1-m(1)).^2+(x2-m(2)).^2));
wfungradfun = @(x1,x2,m,scale) 1/(sigma^2)*wfunfun(x1,x2,m,scale).*[m(1)-x1 m(2)-x2]; 
wfun = @(x1,x2) wfunfun(x1,x2,m1,scale1) + wfunfun(x1,x2,-m1,scale2) + wfunfun(x1,x2,m2,scale3);  %wfunfun(x1,x2,-m2,scale2);
wfungrad = @(x1,x2) wfungradfun(x1,x2,m1,scale1) + wfungradfun(x1,x2,-m1,scale2) + wfungradfun(x1,x2,m2,scale3);  %wfungradfun(x1,x2,-m2,scale2);

fmdl = precomputeFEM_DG(meshpar_fine);
fmdl = precomputeRHS_DG(meshpar_fine,fmdl,wfun,wfungrad);

% Precomputing stiffness
fmdl = fixingD(meshpar_fine,fmdl,D);

% Evaluating forward model from precomputed matrices
N = 13;
trunc = N*(2*N+1);
fmdl = computeProjectionMatrices_fine(fmdl,meshpar_fine,trunc);
u = evalFowardModel(fmdl,meshpar_fine,gamma');

%% Add source back again
U = zeros(pN,1);
U(meshpar_fine.NZ) = u;

u = U + wfun(meshpar_fine.p(1,:)',meshpar_fine.p(2,:)');


%% Data
u = u(meshpar_fine.t(1:3,:));
u = fmdl.G*u(:);
u = reshape(u,[3 meshpar_fine.HN]);
        
        
data = gamma.*u;
dataq = data(1,:)*fmdl.U_proj1 + data(2,:)*fmdl.U_proj2 + data(3,:)*fmdl.U_proj3;

%% make noise realization
rng(1);
noiselevel = 0.04;
xi = randn(trunc,1);
%rel_noise_level = 0.01;
eps = noiselevel * norm(dataq,2) / norm(xi,2);
epssq = eps^2;
eps = eps*xi;

load('/work3/akara/qPAT-level/Data/noise_model/eigenv/E_coarse_0.01_0.0175.mat')
N = 13;
trunc = (2*N+1)*N;
datav = datapar.bq'*E_coarse(:,1:trunc)';
epsv = eps'*E_coarse(:,1:trunc)';

%% Figure number 2
figure(1);
ha = tight_subplot(1,2,[.05 .05],[.1 .1],[.1 .1]);
axes(ha(1))
plot_from_gamma(datav,meshpar);
colormap default
view(2)
%title('Absorption coefficient, $\gamma$','interpreter','latex','fontsize',200)
xticks([-1,0,1])
yticks([-1,0,1])
%xlim([0,2*pi])

set(gca,'XTickLabel',{'-1','0','1'},'fontsize',15)
colorbar
cbh = colorbar('XTick', [0,0.5,1,1.5,2],'fontsize',11);

box on
axis square
ax=gca;
ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 11;
%subplot(132)
%plot_from_gamma(full(mus),meshpar_fine)
%title('Scattering coefficient')
%colorbar
axes(ha(2))
plot_from_gamma(epsv,meshpar)
%hT = title('Diffusion coefficient, a','interpreter','latex','fontsize',200);
xticks([-1,0,1])
yticks([-1,0,1])
%xlim([0,2*pi])
set(gca,'XTickLabel',{'-1','0','1'},'fontsize',15)
grid off
colorbar
set(gcf, 'Position',  [100, 100, 1000, 360])
set(gcf,'color','w');
box on
axis square
ax=gca;
ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 11;
cbh = colorbar('XTick', 0.1*(-1:0.25:1),'fontsize',11);
%%


%% Export figure
%export_fig 'Figures/data_plot/data.pdf' -opengl
print('Figures/data_plot/data1.eps','-depsc2')