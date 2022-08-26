function Sol = f2Da3D(MP1, MP2, pi1, pi2)

% MP son las matrices de proyección
% pi1 son las coordenadas XY en la cámara 1
% pi2 son las coordenadas XY en la cámara 2


% Resolvemos el sistema 
Resol = [MP1(1,1:3)-pi1(1)*MP1(3,1:3); MP1(2,1:3)-pi1(2)*MP1(3,1:3);...
    MP2(1,1:3)-pi2(1)*MP2(3,1:3); MP2(2,1:3)-pi2(2)*MP2(3,1:3)];

Indep = -[MP1(1,4)-pi1(1)*MP1(3,4); MP1(2,4)-pi1(2)*MP1(3,4);...
    MP2(1,4)-pi2(1)*MP2(3,4); MP2(2,4)-pi2(2)*MP2(3,4)];

Sol(1,:) = Resol([1 2 3],:)\Indep([1 2 3]);
Sol(2,:) = Resol([1 2 4],:)\Indep([1 2 4]);
Sol(3,:) = Resol([1 3 4],:)\Indep([1 3 4]);
Sol(4,:) = Resol([2 3 4],:)\Indep([2 3 4]);

end

