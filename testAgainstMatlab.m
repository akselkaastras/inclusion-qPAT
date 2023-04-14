%%% Solve pde in MATLAB
model = createpde();
geometryFromEdges(model,@circleg);

figure 
pdegplot(model,"EdgeLabels","on"); 
axis equal

myFun = @(location,state) 3+2*location.x; 

applyBoundaryCondition(model,"dirichlet", ...
                             "Edge",1:model.Geometry.NumEdges, ...
                             "u",myFun);

specifyCoefficients(model,"m",0,"d",0,"c",0.001,"a",5,"f",0);

hmax = 0.05;
generateMesh(model,"Hmax",hmax);
figure
pdemesh(model); 
axis equal

results = solvepde(model);

%%
u = results.NodalSolution;
pdeplot(model,"XYData",u)
title("Numerical Solution");
xlabel("x")
ylabel("y")
hold on
plot(model.Mesh.Nodes(1,I),model.Mesh.Nodes(2,I),'rx')
