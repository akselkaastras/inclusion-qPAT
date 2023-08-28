function G = computePlanes(meshpar)

G = sparse(3*meshpar.HN,3*meshpar.HN);
for i = 1:meshpar.HN
    if rem(i,round(meshpar.HN/20)) == 0
        disp(['- - ', num2str(round((i/meshpar.HN)*100)),' %'])
    end
    x = meshpar.p(1,meshpar.t(1:3,i))';
    y = meshpar.p(2,meshpar.t(1:3,i))';
    G_elem = computePlaneMatrix(x,y);
    G((i-1)*3+(1:3),(i-1)*3+(1:3)) = G_elem;
end