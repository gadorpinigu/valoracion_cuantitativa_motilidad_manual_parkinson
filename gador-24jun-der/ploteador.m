close all
clear all
clc



csvs=dir('*.csv');
mp4s=dir('*_corregido.mp4');
numfiles=length(csvs);
nombrescsv=[];
nombresmp4=[];
for i=1:numfiles
    nombrescsv{i}=csvs(i).name;
    nombresmp4{i}=mp4s(i).name;
end
%Hasta aqui ya tengo los nombres de los archivos que voy a analizar para
%poer representar la trayectoria obtenida por DLC sobre los vídeos
%originales
%Están ordenados (el i-esimo elemento de ada una de las carpetas coincide
%con el csv correpsondiente al mp4 del que estamos hablando

for j=1:numfiles
    csv=nombrescsv{j}
    datos=open(nombrescsv{j});
    datos=datos.data;
    indice=datos(:,2:3);
    V=VideoReader(nombresmp4{j});
    nombre_nuevo=[strcat(nombresmp4{i}(1:(end-4)),'_plot')]
    props=get(V);
    num_frames=props.NumFrames;
    v=VideoWriter(nombre_nuevo,'MPEG-4');
    v.FrameRate=props.FrameRate;
    open(v)
    for l=1:num_frames
        I=read(V,l);
        figure(1)
        imshow(I)
        hold on
        plot(indice(l,1),indice(l,2),'r*')
        frame=getframe(gcf);
        writeVideo(v,frame)
    end
close(v)
end




