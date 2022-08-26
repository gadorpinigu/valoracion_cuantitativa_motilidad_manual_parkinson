function [MP1,MP2, Dif] = calibrado(p1,p2)

% Función para calcular las matrices de proyección para la conversión de
% coordenadas 2D de las cámaras a la 3D real, resolviendo
% el sistema por descomposición en valores singulares (svd)

% p1 y p2 son las coordenadas de los puntos del maniquí, en el sistema de
% coordenadas 2D de las cámaras. p1(punto,xy)


% Coordenadas reales conocidas de los puntos de referencia.

%Coordenadas 21jul22
p(1,:)=[3.5,1.0,15.0];
p(2,:)=[28.4,2.7,15.0];
p(3,:)=[15.8,17.2,15.0];
p(4,:)=[16.6,22.7,4.9];
p(5,:)=[2.9,38.0,4.9];
p(6,:)=[29.2,37.3,4.9];

%Coordenadas 23/jun/2022
% p(1,:)=[3.5,38.0,15.0];
% p(2,:)=[28.3,36.8,15.0];
% p(3,:)=[15.8,22.6,15.0];
% p(4,:)=[16.7,17.2,4.9];
% p(5,:)=[3.0,2.0,4.9];
% p(6,:)=[29.2,2.6,4.9];

% maniquí metálico tumbado
% p(1,:)=[35,380,150];
% p(2,:)=[283.0,368.0,150];
% p(3,:)=[158,226,150];
% p(4,:)=[167,172,50];
% p(5,:)=[30,20,50]; 
% p(6,:)=[292,26,50];
% 

% %Coordenada de caja metalica alta
% p(1,:)=[0,0,15.10];
% p(2,:)=[24.80,0,15.10];
% p(3,:)=[12.4,-17.20,15.10];
% p(4,:)=[12.40,-20.20,5.20];
% p(5,:)=[0,-35.20,5.20];
% p(6,:)=[24.80,-35.20,5.20];



% % maniquí metálico de pie
% p(6,:)=[29,20,50];
% p(5,:)=[35,380,150];
% p(4,:)=[170,128.5,50]; %el de arriba
% p(3,:)=[157,227,150]; 
% p(2,:)=[290,27,50];
% p(1,:)=[283.5,372.5,150];

% maniquí de madera a 3/12/2016 (S. Javier)
% p(1,:)=[38.5,47.8,50.9]; 
% p(2,:)=[38.8,354.7,101.4];
% p(3,:)=[176.8,151.6,51.8];
% p(4,:)=[147.4,247.8,101.3];
% p(5,:)=[281.6,53.5,51.9];
% p(6,:)=[286.9,353.0,100.8];

% % maniquí de madera versión 2, a 14/1/17
% p(1,:)=[57.1,51,51.2]; 
% p(2,:)=[67.3,346.8,149.5];
% p(3,:)=[196.2,153.9,52.2];
% p(4,:)=[196.8,258.3,150.4];
% p(5,:)=[300.1,55.9,50.6];
% p(6,:)=[302.6,350.2,150];
% load('maniqui_madera_2.mat');

% Calibracion en modelo papel en el labo de Fisica

% load('calibrado_papel.mat');

% MR es la matriz que utilizaremos para poder realizar la descomposición 
% en valores singulares
% MR = [p(1,:),1; p(2,:),1; p(3,:),1; p(4,:),1; p(5,:),1; p(6,:),1];
MR = [p, ones(6,1)];

% Z es una matriz con ceros para poder contruir la matriz del sistema homogeneo
Z = zeros(6,4);
I = ones(1,4);
% Matriz del sistema homogeneo que se resuelve con la descomposición en 
% valores singulares (svd)
% MF1 = [MR, Z, -MR.*(p1(:,1)*[1, 1, 1, 1]); Z, MR, MR.*(p1(:,2)*[1, 1, 1, 1])];
MF1 = [MR, Z, -MR.*(p1(:,1)*I); Z, MR, MR.*(p1(:,2)*I)];
% Resolvemos por la descomposición en valores singulares (svd)
[~,~,V1] = svd(MF1);

% EMP1 = V1(:,size(MF1,2));  % Montamos la matriz de proyección de la imagen 1
% MP1(1,:) = EMP1(1:4)';
% MP1(2,:) = -EMP1(5:8)';
% MP1(3,:) = EMP1(9:12)';

EMP1 = V1(:, size(MF1, 2))';
MP1 = [-EMP1(1:4); EMP1(5:8); -EMP1(9:12)];

% Lo mismo para los puntos de la cámara 2

%MF2 = [MR, Z, -MR.*(p2(:,1)*[1, 1, 1, 1]); Z, MR, MR.*(p2(:,2)*[1, 1, 1, 1])];
MF2 = [MR, Z, -MR.*(p2(:,1)*I); Z, MR, MR.*(p2(:,2)*I)];

[~,~,V2] = svd(MF2);
EMP2 = V2(:, size(MF2, 2))';
MP2 = [EMP2(1:4); -EMP2(5:8); EMP2(9:12)];

% EMP2 = V2(:,size(MF2,2));  % Montamos la matriz de proyección de la imagen 2
% 
% MP2(1,:) = EMP2(1:4)';  
% MP2(2,:) = -EMP2(5:8)';
% MP2(3,:) = EMP2(9:12)';

% % Guardamos las matrices de proyección
% save(('Calib.mat'),'MP1', 'MP2')
% 
% a ver qué tal ha quedado...
% parece que las soluciones buenas son la 1 y la 3
Dif = zeros(4,3,6);
    for i = 1:6
%         pi1 = [p1(i,1), p1(i,2)];  % de la imagen 1
%         pi2 = [p2(i,1), p2(i,2)];  % de la imagen 2
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