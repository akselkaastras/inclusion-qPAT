function plot_from_gamma(gamma,meshpar)

trisurf(meshpar.t(1:3,:)', meshpar.p(1, :), meshpar.p(2, :), gamma,'EdgeColor','none','FaceColor','interp')
view(2)