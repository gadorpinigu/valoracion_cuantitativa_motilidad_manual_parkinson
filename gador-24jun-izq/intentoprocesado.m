clc
clf
close all
clear all
%%
archivos_mat=dir('*.mat');
numfiles=length(archivos_mat);
nombres_mat=[];

for i=1:numfiles
    nombres_mat{i}=archivos_mat(i).name;
    if length(nombres_mat{i})==13
        matrices=load(nombres_mat{i});
        MP1=matrices.MP1*10;
        MP2=matrices.MP2*10;

    else
        tiempos_mat=load(nombres_mat{i})
    end
end
%%
tiempos=tiempos_mat.time_frame;
for i=1:length(tiempos)
    for j=1:2
t(i,j)=str2num(tiempos{j,i}(17:end)); %columna 1 tiempo camara 1 y columna 2 tiempo camara 2
    end
end
t1=t(2:end,1);
t2=t(2:end,2);
%% 
files=dir('*.csv');
numfiles=length(files);
nombres_archivos=[];
p=zeros(1000,4);
iter=1;
for i=1:numfiles
    nombres_archivos{i}=files(i).name
    datos=open(nombres_archivos{i});
    datos=datos.data;
    p(1:length(datos),iter:(iter+1))=datos(1:length(datos),2:3);
    iter=3;
end
p=p(1:length(datos),:); %col1ycol2=puntos cam1 (x1,y1)
%col3ycol4=puntos cam2 (x2,y2)
x1=p(:,1);
y1=p(:,2);
x2=p(:,3);
y2=p(:,4);

%% Supongamos que escogemos t1
x21=interp1(t2,x2,t1);
y21=interp1(t2,y2,t1);
r1=[x1 y1];
r2=[x21 y21];
r2(1,:)=r2(2,:);
t1=t1(:)-min(t1);
%Ya tengo los valores de la  camara 2 para tiempos 1
%%
for i=1:length(x1)
    Sol=f2Da3D(MP1,MP2,r1(i,:),r2(i,:));
    Sol_def(i,:)=Sol(1,:);
end
%%
x3d=Sol_def(:,1)-min(Sol_def(:,1));
y3d=Sol_def(:,2)-min(Sol_def(:,2));
z3d=Sol_def(:,3)-min(Sol_def(:,3));

figure(1)
% subplot(1,3,1)
% plot(t1,x3d,'r')
% subplot(1,3,2)
% plot(t1,y3d,'g')
% subplot(1,3,3)
% plot(t1,z3d,'b')
plot(t1,x3d,'-.r')

%%
b=1;
i=1;
amplitudes=[];
periodo=[];
while b==1
    [x,y,b]=ginput(2)
    amplitudes(i)=abs(y(2)-y(1));
    periodo(i)=abs(x(2)-x(1))*2;
    i=i+1;
end
%%
Ax=amplitudes(1:end-1);
Tx=periodo(1:end-1);

parametros_x=[transpose(Ax) transpose(Tx)];
%% Estadística sobre los parametros

