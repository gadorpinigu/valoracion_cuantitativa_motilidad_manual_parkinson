% % % SCRIPT ORIGINAL: PROYECTO TTKT
% ADAPTADO PARA TEST PARKINSON: GADOR PIÑEYRO, UNAV

function [MP1,MP2, Dif] = calibrado(p1,p2)

% Función para calcular las matrices de proyección para la conversión de
% coordenadas 2D de las cámaras a la 3D real, resolviendo
% el sistema por descomposición en valores singulares (svd)

% p1 y p2 son las coordenadas de los puntos del maniquí, en el sistema de
% coordenadas 2D de las cámaras. p1(punto,xy)

%Coordenadas 21jul22
p(1,:)=[3.5,1.0,15.0];
p(2,:)=[28.4,2.7,15.0];
p(3,:)=[15.8,17.2,15.0];
p(4,:)=[16.6,22.7,4.9];
p(5,:)=[2.9,38.0,4.9];
p(6,:)=[29.2,37.3,4.9];

% MR es la matriz que utilizaremos para poder realizar la descomposición 
% en valores singulares
% MR = [p(1,:),1; p(2,:),1; p(3,:),1; p(4,:),1; p(5,:),1; p(6,:),1];
MR = [p, ones(6,1)];

% Z es una matriz con ceros para poder contruir la matriz del sistema homogeneo
Z = zeros(6,4);
I = ones(1,4);

% Matriz del sistema homogeneo que se resuelve con la descomposición en 
% valores singulares (svd)
MF1 = [MR, Z, -MR.*(p1(:,1)*I); Z, MR, MR.*(p1(:,2)*I)];

% Resolvemos por la descomposición en valores singulares (svd)
[~,~,V1] = svd(MF1);

EMP1 = V1(:, size(MF1, 2))';
MP1 = [-EMP1(1:4); EMP1(5:8); -EMP1(9:12)];

% Lo mismo para los puntos de la cámara 2

MF2 = [MR, Z, -MR.*(p2(:,1)*I); Z, MR, MR.*(p2(:,2)*I)];

[~,~,V2] = svd(MF2);
EMP2 = V2(:, size(MF2, 2))';
MP2 = [EMP2(1:4); -EMP2(5:8); EMP2(9:12)];

Dif = zeros(4,3,6);
    for i = 1:6
        % para poder calcular sus coordenadas en 3
        % Utilizamos una funcíon de 2D a 3D
        Sol = f2Da3D(MP1,MP2,[p1(i,1) p1(i,2)],[p2(i,1) p2(i,2)]);
        Sol([1 3],:)
        for j = 1:4
            Sol(j,:) = Sol(j,:)-p(i,:);
        end
        Dif(:,:,i) = Sol(:,:);
    end % del for i
end
