% % % SCRIPT ORIGINAL: GADOR PIÑEYRO, TFM PARKINSON, UNAV

% % % % % % % % % % POST PROCESAMIENTO % % % % % % % % % %

close all
clear all
clc
%%
path=pwd;
carpeta=input('Introduce el número de historia clínica: ','s');
cd(strcat(path,'/',carpeta)) %PARA MAC
%cd(strcat(path,'\',carpeta)) %PARA WINDOWS

if exist('calibrado.mat')==0;
    sprintf('Hace falta realizar la calibración.')
    sprintf('Por favor, marque en orden los 6 puntos para el calibrado en cada imagen.')
    tiff_cam1=dir(strcat(carpeta,'*Cam1.tiff')).name;
    tiff_cam2=dir(strcat(carpeta,'*Cam2.tiff')).name;
    Icam1=imread(tiff_cam1);
    Icam2=imread(tiff_cam2);

    p1 = correlacion(Icam1,1);
    p2 = correlacion(Icam2,2);

    sb = [1 3];
    [MP1,MP2, Dif,p] = calibrado(p1,p2);
    % Revisa la máxima diferencia, para aceptar o no el calibrado
    calidad = max(max(max(Dif(sb,:,:))));
    folder_name=pwd;
    path_calibrado=strcat(folder_name,'/calibrado.mat');

    if calidad > 0.3
        mensaje = sprintf('Hay diferencias de hasta %0.1f en la posición de los puntos. No se acepta el calibrado',calidad);
    else
        mensaje = sprintf('Diferencias máximas de %0.2f en la posición de los puntos. Calibrado aceptado',calidad);
        save(path_calibrado,'MP1','MP2')
    end

    load('calibrado.mat')
    opc=input('¿Quieres comprobar que calcula bien los puntos de calibración? responde "si" o "no": ','s');

    if opc=="si"

        for i=1:6
            Sol=f2Da3D(MP1,MP2,p1(i,:),p2(i,:));
            Sol1(i,:)=Sol(1,:);
        end

        close all
        figure(1)
        subplot(3,1,1)
        plot(Sol1(:,1),'*r')
        hold on
        yline(p(:,1))
        title('Coordenada x')
        legend('Coordenadas calculadas','Coordenadas reales')
        subplot(3,1,2)
        plot(Sol1(:,2),'*r')
        hold on
        yline(p(:,2))
        title('Coordenada y')
        legend('Coordenadas calculadas','Coordenadas reales')
        subplot(3,1,3)
        plot(Sol1(:,3),'*r')
        hold on
        yline(p(:,3))
        title('Coordenada z')
        legend('Coordenadas calculadas','Coordenadas reales')
        savefig('Comprobación_calibrado')
    else
    end
    % % OPCIONAL: COMPROBACIÓN DE QUE SE HAN CALCULADO BIEN LOS PUNTOS DE CALIBRACIÓN POST CALIBRACION

else

    load('calibrado.mat')


end


tarea=input('Introduce el número de tarea que desea evaluar: ','s');
num_tarea=strcat('Tarea',tarea);

ext=strcat(carpeta,'*',num_tarea,'*times.mat');
tiempos=dir(ext);
load(tiempos.name)

for i=1:length(time_frame)
    for j=1:2
        t(i,j)=str2num(time_frame{1,i}(17:end));
    end
end

t1=t(1:end-1,1);
t2=t(1:end-1,2);

% % % % REALMENTE QUEREMOS QUE SEAN PARA LOS MISMOS TIEMPOS ASI QUE
% ESCOGEMOS T1 E INTERPOLAMOS PARA OBTENER X2 Y2
ext2=strcat(carpeta,'*',num_tarea,'*.csv');
csv=dir(ext2);
csvs=[];
mp4=[];
if tarea=='1'|tarea=='2'
    repes=1;
else
    repes=4;
end

list_dedos={'indice','corazon','anular','menique'};
bp=2;
iter=1;
for bodypart=1:repes
    t1=t(1:end-1,1);
    dedo=zeros(length(t1),4);
    for i=1:length(csv)
        csvs{i}=csv(i).name;
        datos=open(csvs{i});
        if i==1
            dedo(:,1:2)=datos.data(:,bp:(bp+1));
        else
            dedo(:,3:4)=datos.data(:,bp:(bp+1));
        end
    end
    x1=dedo(:,1);
    y1=dedo(:,2);
    x2=dedo(:,3);
    y2=dedo(:,4);
    x2_t1=interp1(t2,x2,t1);
    y2_t1=interp1(t2,y2,t1);

    t1=t1-min(t1);

    r1=[x1,y1];
    r2=[x2_t1,y2_t1];

    for i=1:length(x1)
        Sol=f2Da3D(MP1,MP2,r1(i,:),r2(i,:));
        Sol1(i,:)=Sol(1,:);
        Sol2(i,:)=Sol(2,:);
        Sol3(i,:)=Sol(3,:);
        Sol4(i,:)=Sol(4,:);
    end
    close all
    figure(1)
    subplot(3,1,1)
    plot(t1,Sol1(:,1),'r.-')
    hold on
    plot(t1,Sol3(:,1),'b.-')
    title('Coordenada x')
    xlabel('Tiempo')
    ylabel('Posición')
    subplot(3,1,2)
    plot(t1,Sol1(:,2),'r.-')
    hold on
    plot(t1,Sol3(:,2),'b.-')
    title('Coordenada y')
    xlabel('Tiempo')
    ylabel('Posición')
    subplot(3,1,3)
    plot(t1,Sol1(:,3),'r.-')
    hold on
    plot(t1,Sol3(:,3),'b.-')
    title('Coordenada z')
    xlabel('Tiempo')
    ylabel('Posición')
    linkaxes

    savefig(strcat(num_tarea,list_dedos{iter},'_3d'))
    bp=bp+3;

    

    iter=iter+1;
    pause
end



visualizacion_grafica=input('¿Desea graficar sobre los vídeos las posiciones obtenidas? Responda "si" o "no": ','s');

if visualizacion_grafica=="si"
    close all
    path=strcat(carpeta,'*',num_tarea,'*_corregido.mp4');
    mp4s=dir(path);
    numfiles=length(csvs);
    nombrescsv=[];
    nombresmp4=[];

    for i=1:numfiles
        nombrescsv{i}=csvs{i};
        nombresmp4{i}=mp4s(i).name;
    end

    for j=1:numfiles
        csv=nombrescsv{j}
        datos=open(nombrescsv{j});
        datos=datos.data;
        indice=datos(:,2:3);
        corazon=datos(:,5:6);
        anular=datos(:,8:9);
        menique=datos(:,11:12);
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
            if repes==4
                plot(corazon(l,1),corazon(l,2),'b*')
                plot(anular(l,1),anular(l,2),'g*')
                plot(menique(l,1),menique(l,2),'m*')
            else
            end
            frame=getframe(gcf);
            writeVideo(v,frame)
        end

        close(v)
    end

else

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % %   FUNCIONES    % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function p = correlacion(A, n_fig)

%muestra la imagen problema
figure(n_fig);
imshow(uint8(A));
r=4;
% pide que se introduzcan manualmente los
% que faltan

    for i = 1:6
        % pide la posición aproximada del punto
        [x0,y0] = ginput(1);
        % busca la posición del centroide del punto
        p(i,:) = contrasteJ(x0,y0,A,r,r);
    
    % muestra la posición calculada
    figure(2)
    figure(n_fig)
    hold on
        plot(p(i,1),p(i,2),'xr')
    end

    figure(n_fig);
    for i = 1:6
        text(p(i,1),p(i,2)+2,int2str(i));
     end

end

function Puntos = contrasteJ(x, y, imagen, semiDistP, semiDistG)

%%%%%%%%%%%%%%%%%%% Empieza Contraste...
 %semiDistP=6;
 %semiDistG=dist;
% if (dist==0)
%     semiDistG=12;
%    
% else
%     semiDistG=dist;
%     
% end

%%%%%%%%%%%% PARA ROTACION CON VELOCIDAD DE HASTA UNOS 200 mm / s
% semiDistG=25;
% semiDistP=10;
[XG,YG] = meshgrid(-semiDistG:semiDistG, -semiDistG:semiDistG);
[X,Y] = meshgrid(-semiDistP:semiDistP, -semiDistP:semiDistP);
[X2,Y2] = meshgrid(-semiDistP:.2:semiDistP, -semiDistP:.2:semiDistP);
% DEFINO MASCARAS PARA TOMAR SOLO LO QUE TIENE RADIO INF A R
MascaraRadioG = sqrt(XG.^2+YG.^2);
MascaraRadioG = round((sign(semiDistG-MascaraRadioG)+1)/2);
% MascaraRadio = sqrt(X.^2+Y.^2);
% MascaraRadio = round((sign(semiDistP-MascaraRadio)+1)/2);
MascaraRadio2 = sqrt(X2.^2+Y2.^2);
MascaraRadio2 = round((sign(semiDistP-MascaraRadio2)+1)/2);

% Obtenemos recuadros donde se encuentra el punto y obtenemos una imagen de
% ese recuadro. Busco el mínimo....
RangoXG = round(x)-semiDistG:round(x)+semiDistG;
RangoYG = round(y)-semiDistG:round(y)+semiDistG;
%disp(['round(x) ' num2str(round(x)) '; round(y) ' num2str(round(y)) '; semiDistG ' num2str(semiDistG)]);
R = double(imagen(RangoYG,RangoXG,:));
R = (squeeze((R(:,:,1)+R(:,:,2)+R(:,:,3))/3));
R = (255-R).*MascaraRadioG;
[i,j] = find(R==max(R(:)));
% Puntos = [mean(mean(XG(i,j)))+round(x), mean(mean(YG(i,j)))+round(y)];
% x = Puntos(1);
% y = Puntos(2);
x = mean(mean(XG(i,j)))+round(x);
y = mean(mean(YG(i,j)))+round(y);
%
% Ahora repito el cï¿½lculo centrado en el punto del mï¿½nimo...
%
RangoX = round(x)-semiDistP:round(x)+semiDistP;
RangoY = round(y)-semiDistP:round(y)+semiDistP;
R = double(imagen(RangoY,RangoX,:));
R = squeeze((R(:,:,1)+R(:,:,2)+R(:,:,3))/3);

% Conseguimos hacer un ajuste del contraste para cada punto seleccionado.
ValorMax = mean(R(:))+.333*std(R(:));
R(R>ValorMax) = ValorMax;
F = 1.001-(R-min(R(:)))./(ValorMax-min(R(:)));
F2 = interp2(X,Y,F,X2,Y2,'cubic').*MascaraRadio2;

% Defino Mascara...
M = round((sign(F2-0.25)+1)/2);
F2 = F2.*M;


Puntos = [sum(X2(:).*F2(:))./sum(F2(:))+round(x), sum(Y2(:).*F2(:))./sum(F2(:))+round(y)];

end
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
function [MP1,MP2, Dif,p] = calibrado(p1,p2)

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


MR = [p, ones(6,1)];
Z = zeros(6,4);
I = ones(1,4);
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
