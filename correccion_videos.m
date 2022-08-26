% % % SCRIPT ORIGINAL: GADOR PIÑEYRO, TFM PARKINSON, UNAV

% % % % % % % % % CORRECCIÓN DE LÍNEAS
close all
clear all
clc

path=pwd;
carpeta=input('Introduce el número de historia clínica: ','s');
%cd(strcat(path,'\',carpeta)) % WINDOWS
cd(strcat(path,'/',carpeta)) % MAC

comprobacion=strcat(carpeta,'*_corregido.mp4');
if exist(comprobacion)==0
% for i=1:length(comprobacion)
%     nom{i}=comprobacion(i).name;
% end
% 
% if isempty(nom)~=0
    
    mp4files = dir('*.mp4');
    numfiles = length(mp4files);
    nombresarchivos=[];
    
    for i=1:numfiles
        nombresarchivos{i}=mp4files(i).name;
        nombre_nuevo=[strcat(nombresarchivos{i}(1:(end-4)),'_corregido')];
        V=VideoReader(nombresarchivos{i});
        I=read(V,1);
        [m,n,p]=size(I);
        props_cam1=get(V);
             num_frames=props_cam1.NumFrames; % EN ESTE MATLAB NO ESTA ACTUALIZADO
        %     Y NO EXISTE LA PROPIEDAD NUMFRAMES
        %num_frames=floor(props_cam1.Duration*props_cam1.FrameRate);
        v=VideoWriter(nombre_nuevo,'MPEG-4');
        v.FrameRate=props_cam1.FrameRate;
        open(v)
        for j=1:num_frames
            Ii=read(V,j);
            I2=Ii;
            for k=1:3
                for M=2:2:m-1
                    I2(M,:,k)=mean([I2(M-1,:,k); I2(M+1,:,k)]);
                end
            end
            writeVideo(v,I2)
        end
        close(v)
    end
else
    
end
