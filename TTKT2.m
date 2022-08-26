% Begin initialization code - DO NOT EDIT
function varargout = TTKT2(varargin)
    
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @TTKT2_OpeningFcn, ...
                       'gui_OutputFcn',  @TTKT2_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
                   
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1}); %03824225
                                                        %197001
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
     
end
% End initialization code - DO NOT EDIT

% --- Executes just before TTKT2 is made visible.
function TTKT2_OpeningFcn(hObject, ~, handles, varargin)

% cada vez que se ejecuta en un nuevo usuario, hay que dar de alta en
% MatLab la dll. Escribiendo esto:
% imaqregister('C:\Program Files (x86)\TIS IMAQ for MATLAB R2013b\x64\TISImaq_R2013.dll')
% ver http://www.theimagingsource.com/support/documentation/ic-matlab-2013-extension/

% PARTE GENERADA AUTOMÁTICAMENTE
% Choose default command line output for TTKT2
    handles.output = hObject;
% Update handles structure
    guidata(hObject, handles);
    
    %set(handles.figure1, 'Position', get(0, 'Screensize'));

% PARTE MODIFICADA
    global vid sb dir_res mi_axes;
    
    %set(handles.figure1, 'Position', get(0, 'Screensize')); % metido por
    %Alejandro para aumentar la ventana. Lo quito porque fastidia.
    dir_res = 'C:\Users\imed\Desktop\ttkt_neuro\RESULTADOS\';
    sb = [1 3]; %de las cuatro soluciones, las dos más consistentes
    imaqreset;
    mi_axes = {handles.axes1 handles.axes2};


% define los objetos de video

    vid{1} = TTK_videoinput(1);
    vid{2} = TTK_videoinput(2);
    
% Inicia las cámaras, que quedan pendientes de disparo con trigger

    start(vid{1});
    start(vid{2});


% muestra los vídeos en la interfaz
%cb_live_Callback(hObject, eventdata, handles);
    
    for i = 1:2       
        axes(mi_axes{i})     
        vidRes = vid{i}.VideoResolution;
        nBands = vid{i}.NumberOfBands;
        hImage = image(zeros(vidRes(1),vidRes(2),nBands)); 
        setappdata(hImage,'UpdatePreviewWindowFcn',@gira_ima);
        preview(vid{i}, hImage);   
    end

end

function gira_ima(~, event, himage)

    Img = rot90(event.Data,2);
    set(himage, 'cdata', Img);
    
end

% --- Outputs from this function are returned to the command line.
function varargout = TTKT2_OutputFcn(~, ~, handles)

    varargout{1} = handles.output;
    
end

% --- Executes on button press in PB_calibrar.
function PB_calibrar_Callback(hObject, ~, handles)

    global vid Dif sb;
% saca dos imágenes de las cámaras

    if isempty(get(handles.txt_nhc,'String')) == 1
        warndlg('Introduzca el número de historia clínica');

        return;
    end

    for i=1:2
        I = getsnapshot(vid{i});
%         IM{i}=I;
        IM{i} = imrotate(I,180);
    end
    
    % control de NHC y carpeta
    WD = cd; %Current directory
    folder_name=strcat(WD,filesep,get(handles.txt_nhc,'String'));
    
    if exist(folder_name,'dir')~=7% si la carpeta no existe se crea
        mkdir(folder_name)
    end
    
    % Si se corren las tereas 3 y 4, se añade la repeticion al nombre 
    % Control de filename para evitar overwriting añadiendo _02,03,04...
    outFile_cam1 = avoidOverwrite(strcat(get(handles.txt_nhc,'String'),'_Calibracion_Cam1.tiff'),folder_name,2,2);
    outFile_cam2 = avoidOverwrite(strcat(get(handles.txt_nhc,'String'),'_Calibracion_Cam2.tiff'),folder_name,2,2);
    name_cam1=strcat(folder_name,filesep,outFile_cam1);
    name_cam2=strcat(folder_name,filesep,outFile_cam2);    
    
    % camara 1
    imwrite(IM{1},name_cam1)
    imwrite(IM{2},name_cam2)
    
% % emplea la función correlacion.m para buscar los puntos. Saca la posición
% % de los puntos, ordenados, en un sistema de coordenadas centrado en la
% % imagen y con coordenadas X e Y en el sentido convencional.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ESTO ES LO QUE ESTABA COMENTADO PARA PROBAR CON LOS SENSORES INFRARROJOS%
% ES LA **CALIBRACIÓN** INICIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p1 = correlacion(IM{1},1);
    p2 = correlacion(IM{2},2);

    % p1 = correlacion(IM{1},'Patron_D.png',1);
    % p2 = correlacion(IM{2},'Patron_I.png',2);
    % p1=IM{1};
    % p2=IM{2}; 
    
%save('centros.mat','p1','p2');
%load('centros.mat');

% la función calibrado.m devuelve las matrices de calibrado (MP1 y MP2), y
% las diferencias entre las posiciones reales de los seis puntos y las
% calculadas con el algoritmo
    [MP1,MP2, Dif] = calibrado(p1,p2);

% Revisa la máxima diferencia, para aceptar o no el calibrado
    calidad = max(max(max(Dif(sb,:,:))));
    path_calibrado=strcat(folder_name,'/calibrado.mat');
    if calidad > 0.3 %estaba en 0.3 pero lo he puesto asi para mirar una cosa
        mensaje = sprintf('Hay diferencias de hasta %0.1f en la posición de los puntos. No se acepta el calibrado',calidad);
    else
        mensaje = sprintf('Diferencias máximas de %0.2f en la posición de los puntos. Calibrado aceptado',calidad);
        save(path_calibrado,'MP1','MP2')
    end
          
%set(handles.text1,'String',mensaje);
    uiwait(msgbox(mensaje,'Resultados del calibrado','modal'));
  
%msgbox(mensaje,'Resultados del calibrado');
    close(1);
    close(2);

% En caso de que haya que cambiar los puntos, PARA CALIBRAR EL MANIQUI
% comentar las lineas desde [MP1,... hasta aquí, y descomentar las
% siguientes
load('calibrado.mat');
for i=1:6
    % guardamos las coordenadas del primer punto seleccionada
    pi1=[p1(i,1) p1(i,2)];  % de la imagen 1
    pi2=[p2(i,1) p2(i,2)];  % de la imagen 2
    % para poder calcular sus coordenadas en 3
    % Utilizamos una funcíon de 2D a 3D
    [Sol]=f2Da3D(MP1,MP2,pi1,pi2);
    Dif(:,:,i)=Sol(:,:);
end % del for i


end

% --- Executes on button press in PB_start.
function PB_start_Callback(~, ~, handles)

    global vid dir_res mi_axes sigo sb lim ref referencia 
    clear time_frame name_cam1 V_cam1 name_cam2 V_cam2;
    lim = [0, 0]; 
    ref = 0;
    referencia = 0;
    
    set(handles.PB_stop,'Visible','on');
    set(handles.PB_start,'Visible','off');

    %% MIGUEL -Inicializar Video
    WD = cd; %Current directory
    folder_name=strcat(WD,filesep,get(handles.txt_nhc,'String'));
    
    if exist(folder_name,'dir')~=7% si la carpeta no existe se crea
        mkdir(folder_name)
    end
    
    % Si se corren las tareas 3 y 4, se añade la repeticion al nombre 
    % Control de filename para evitar overwriting añadiendo _02,03,04...

    if strcmp(get(findobj(handles.uibuttongroup1,'Value',1),'tag'),'rb3')
    
    outFile_cam1 = avoidOverwrite(strcat(get(handles.txt_nhc,'String'),'_',get(findobj(handles.uibuttongroup1,'Value',1),'String'),'_rep_',get(handles.txt_rep,'String'),'_Cam1','.mp4'),folder_name,2,2);
    outFile_cam2 = avoidOverwrite(strcat(get(handles.txt_nhc,'String'),'_',get(findobj(handles.uibuttongroup1,'Value',1),'String'),'_rep_',get(handles.txt_rep,'String'),'_Cam2','.mp4'),folder_name,2,2);
   
    elseif strcmp(get(findobj(handles.uibuttongroup1,'Value',1),'tag'),'rb4')
    outFile_cam1 = avoidOverwrite(strcat(get(handles.txt_nhc,'String'),'_',get(findobj(handles.uibuttongroup1,'Value',1),'String'),'_rep_',get(handles.txt_rep,'String'),'_Cam1','.mp4'),folder_name,2,2);
    outFile_cam2 = avoidOverwrite(strcat(get(handles.txt_nhc,'String'),'_',get(findobj(handles.uibuttongroup1,'Value',1),'String'),'_rep_',get(handles.txt_rep,'String'),'_Cam2','.mp4'),folder_name,2,2);

    else 
    outFile_cam1 = avoidOverwrite(strcat(get(handles.txt_nhc,'String'),'_',get(findobj(handles.uibuttongroup1,'Value',1),'String'),'_Cam1','.mp4'),folder_name,2,2);
    outFile_cam2 = avoidOverwrite(strcat(get(handles.txt_nhc,'String'),'_',get(findobj(handles.uibuttongroup1,'Value',1),'String'),'_Cam2','.mp4'),folder_name,2,2);
      
    end

    name_cam1=strcat(folder_name,filesep,outFile_cam1);
    name_cam2=strcat(folder_name,filesep,outFile_cam2);    
    
    V_cam1= VideoWriter(name_cam1,'MPEG-4');
    V_cam2= VideoWriter(name_cam2,'MPEG-4');

    set(V_cam1,'Quality',100)
    set(V_cam2,'Quality',100)
    open(V_cam1)
    open(V_cam2)
%**************************************************************************
%****************** 1. COMPROBACIONES PREVIAS******************************
%**************************************************************************

% desconecta la presentación de video en vivo
% set(handles.cb_live,'value',0);

% Si encuentra datos de calibración, sigue adelante. Si no, pide calibrar
    try
        load('calibrado.mat')   
    catch 
        %set(handles.text1,'String','Faltan datos de calibrado: pulse CALIBRAR');
        errores(1);
        return;
    end

% comprueba que se han rellenado los campos de número de puntos y nhc
% si faltan datos, da un mensaje y se corta
    if isempty(get(handles.txt_np,'String')) == 1
        warndlg('Introduzca el número de puntos');
        set(handles.PB_stop,'Visible','off');
        set(handles.PB_start,'Visible','on');
        return;
    else
        np = str2double(get(handles.txt_np,'String'));
    end

    if isempty(get(handles.txt_nhc,'String')) == 1
        warndlg('Introduzca el número de historia clínica');
        set(handles.PB_stop,'Visible','off');
        set(handles.PB_start,'Visible','on');
        return;
    else
        nhc = get(handles.txt_nhc,'String');
    end


%**************************************************************************
%****************** 2. INICIALIZA VARIABLES *******************************
%**************************************************************************
    sigo = 1;% variable para poder parar el proceso. Cuando se pulsa stop, sigo=0
%Al no conocer el número de fotogramas que se van a capturar, se
%sobredimensionan las matrices 
    P3D = zeros(50000,12,np); %posiciones de los puntos en el sistema 3D
%tasa=zeros(1,50000); % frame rate
    c = zeros(2,2,50000,np);%(XY,cámara,fotograma,punto) Posiciones en las imágenes (2D)
    c0 = zeros(2,2,np);%(XY,cámara,punto) Posiciones corregidas en las imágenes (2D) en cada fotograma
    p_ant = zeros(2,2,np); % posiciones del punto en el último fotograma (x/y,cámara,punto)
%tasa=zeros(1,50000);
    t = zeros(3,50000);
    tnf=zeros(5,50000); %variable analiza los tiempos de ejecución por cada fotograma
    rP = 6; % radio pequeño de búsqueda
    rG = 6; % radio grande de búsqueda
% tamaño de las capturas (lo necesita después para la conversión 2D->3D)
    I = getsnapshot(vid{1});
    Ima = imrotate(I,180);
%     Ima=I;
    [s_y,s_x,~] = size(Ima);
%     cla(handles.axes3);
%     axis(handles.axes3,[5,35,-15,15]); %tamaño inicial
% gráficas animadas para representar desplazamientos en la dirección Z.
% Si el maniquí está tumbado y con la parte alta hacia el scanner, X es
% izda-dcha, Y es cabeza-pies, Z es arriba-abajo
    color = 'rgbkymc';
    set(handles.axes3,'Visible','off')
%     for i=1:np
%        hx{i} = animatedline(handles.axes3,'Color',color(i),'LineWidth',3); 
%     end


% carpeta en la que se escribirán los resultados
    cl = clock;
    ANYO = num2str(cl(1,1));
    MES = num2str(cl(1,2));
    DIA = num2str(cl(1,3));
    f_salida = strcat(nhc,'_',ANYO,'_',MES,'_',DIA,'_');

%Busca en la carpeta RESULTADOS si existe otro archivo de la misma fecha y
%mismo nhc e incrementa en 1 el numero
    Files1 = dir(strcat(dir_res,f_salida,'*.png'));
    nn = int2str(size(Files1,1)+1);
    f_salida = strcat(f_salida,nn);


%**************************************************************************
%****************** 3. POSICIONES INICIALES ******************************
%**************************************************************************

% hace sendas capturas con las cámaras

    for i=1:2 % número de cámara
        I = getsnapshot(vid{i});
        IM{i} = imrotate(I,180);
%         IM{i}=I;
        axes(mi_axes{i}); 
        image(IM{i});%muestra la imagen
    
            for j=1:np % número de punto
                
                [x0,y0] = ginput(1);
                hold on;
                plot(x0,y0,'+g');
                c(:,i,1,j) = contrasteJ(x0,y0,IM{i},rP,rG);
                p_ant(:,i,j) = c(:,i,1,j);
                plot(c(1,i,1,j),c(2,i,1,j),'*r');
                text(c(1,i,1,j)+8,c(2,i,1,j),int2str(j),'color','b');
                hold off;
            end
    end
% guarda una imagen, para poder identificar los puntos en el futuro
    name = strcat(dir_res,f_salida,'_Im.png');
    F = getframe;
    Ima = frame2im(F);
    imwrite(Ima,name);


%**************************************************************************
%****************** 4. SEGUIMIENTO DE PUNTOS ******************************
%**************************************************************************
profile clear 
profile on

    nf = 0; % contador de fotogramas
    t0 = crono();
    
% recordemos: c(XY,cámara,fotograma,punto) Posiciones en las imágenes (2D)
%             p_ant(x/y,cámara,punto)  
    while sigo == 1 % empieza una iteración hasta que se presione stop.
        tic;
        nf = nf+1;
        
    % a. Localiza las posiciones 2D de los puntos en cada cámara
    
   
    % hace un disparo con cada cámara
        trigger([vid{1} vid{2}]); %mirar si quitándolo se soluciona
    
    % lee el disparo de cada cámara
        for i = 1:2 % para cada cámara
            [I, ~, d] = getdata(vid{i});
%              I=permute(I,[2 1 3]);
%              (size(I,1):-1:1,:,:); % SOLUCION 3
%              IM{i}=flip(I); %SOLUCION 2
%              IM{i}=rot90(I); %SOLUCION 1
               IM{i} = imrotate(I,180); % rota la imagen 90º
%                     IM{i}=I;
t_a = d.AbsTime; % tiempo absoluto en que se hizo la adquisición
            t(i,nf) = t_a(4)*3600+t_a(5)*60+t_a(6)-t0;
            
            %% MIGUEL -GrabarVideo
            if i==1
            writeVideo(V_cam1,IM{i})
            time_frame{i,nf} = strcat(num2str(t_a(1)),'_',num2str(t_a(2)), '_',num2str(t_a(3)), '_',num2str(t_a(4)), '_',num2str(t_a(5)), '_',num2str(t_a(6)));
            elseif i==2
            writeVideo(V_cam2,IM{i})
            time_frame{i,nf} = strcat(num2str(t_a(1)),'_',num2str(t_a(2)), '_',num2str(t_a(3)), '_',num2str(t_a(4)), '_',num2str(t_a(5)), '_',num2str(t_a(6)));
            end
        
            
        end

% % % % % %         for i = 1:2 % para cada cámara 
% % % % % %         
% % % % % %             for j = 1:np %itera para cada punto
% % % % % %             % busca la posición del punto poniendo como aproximacion
% % % % % %             % inicial la del último punto encontrado
% % % % % %                 c(:,i,nf,j) = contrasteJ(p_ant(1,i,j),p_ant(2,i,j),IM{i},rP,rG);
% % % % % %             % almacena la posición, para emplearla en la siguiente iteración
% % % % % %                 p_ant(:,i,j) = c(:,i,nf,j);
% % % % % %             end % del j en cada punto
% % % % % %         end % del i en cada cámara
        tnf(1,nf)=toc; 
        tic;
    % b. muestra uno de cada cinco fotogramas
        if mod(nf,5) == 0
            for i = 1:2 % para cada cámara
                axes(mi_axes{i});
                image(IM{i});
                set(gca, 'visible', 'off');
                hold on;         
% % %                 for j = 1:np % para cada punto
% % %                     plot(c(1,i,nf,j),c(2,i,nf,j),'o','Color',color(j));
% % %                 end % del j en puntos
                hold off;
            end % del i en cámaras
        
        end % del if mod(nf,5)==0;
    tnf(2,nf)=toc; 
    tic;
    % corregimos el desfase entre cámaras
        t(3,nf) = min([t(1,nf) t(2,nf)]); %el menor de los dos tiempos t(1,nf) o t(2,nf)
    % c0(XY,cámara,punto): posiciones en el tiempo t(3,nf)
    %delta=1;
        if nf > 1 % sólo corregimos a partir del segundo fotograma
            if t(2,nf) > t(1,nf)
                delta = (t(1,nf)-t(2,nf-1))/(t(2,nf)-t(2,nf-1));
                c0(:,1,:) = c(:,1,nf,:);
                c0(:,2,:) = (1-delta)*c(:,2,nf-1,:)+delta*c(:,2,nf,:);
            else
                delta = (t(2,nf)-t(1,nf-1))/(t(1,nf)-t(1,nf-1));
                c0(:,1,:) = (1-delta)*c(:,1,nf-1,:)+delta*c(:,1,nf,:);
                c0(:,2,:) = c(:,2,nf,:);
            end
        else
            c0(:,1,:) = c(:,1,nf,:);
            c0(:,2,:) = c(:,2,nf,:);
        end

    % c. calcula las posiciones 3D
    % recordemos: c(XY,cámara,fotograma,punto) Posiciones en las imágenes (2D)
    % antes, tiene que convertir de coord matlab a normales
        for i = 1:np % por cada punto
        % invertimos las coordenadas y
            c0(2,:,i) = s_y-c0(2,:,i);
        % trasladamos los ejes de coordenadas al centro de la imagen.
            c0(1,:,i) = c0(1,:,i)-s_x/2;
            c0(2,:,i) = c0(2,:,i)-s_y/2;
            Sol = f2Da3D(MP1,MP2,c0(:,1,i),c0(:,2,i));
            P3D(nf,:,i) = Sol(:);
        end % de la iteración en puntos
        tnf(3,nf)=toc;
        tic;
    
    % d. representamos las coordenadas.
        n_promedio = 100;
        if nf == n_promedio % cuando lleva cien fotogramas, calcula el movimiento medio
            Xmed = zeros(1,np);
            for i = 1:np % por cada punto
                X0(:) = (P3D(1:n_promedio,sb(1)+8,i)+P3D(1:n_promedio,sb(2)+8,i))/2;
                Xmed(i) = mean(X0(:));
            end
%        Y(:,:)=(P3D(1:n_promedio,sb(1)+4,:)+P3D(1:n_promedio,sb(2)+4,:))/2;
%        Ymed=mean(Y,1);
%        Z(:,:)=(P3D(1:n_promedio,sb(1)+8,:)+P3D(1:n_promedio,sb(2)+8,:))/2;
%        Zmed=mean(Z,1);
        elseif nf > n_promedio
            T = t(3,nf);
            X = zeros(1,np);
            for i = 1:np
                X(i) = (P3D(nf,sb(1)+8,i)+P3D(nf,sb(2)+8,i))/2-Xmed(i);
% % % % %                 if lim(1) ~= lim(2) && X(i) >= lim(2) && X(i) <= lim(1) %comprobar limtes
% % % % % %                 bombilla = 1;                     %primero el inferior
% % % % %                     set(handles.text5, 'BackgroundColor', [1 1 0])
% % % % %                 else
% % % % % %                 bombilla = 0;
% % % % %                     set(handles.text5, 'BackgroundColor', [1 0 0])
% % % % %                 end
% % % % %                 if referencia == 1
% % % % %                     %poner handles.axes3
% % % % %                     h1 = refline([0 lim(1)]); %dibuja lineas de referencia de  e el eje Y hasta
% % % % %                     h1.Color = 'b';         %que el plot crece
% % % % %                     h1.LineWidth = 1;
% % % % %                     h2 = refline([0 lim(2)]); 
% % % % %                     h2.Color = 'g';
% % % % %                     h2.LineWidth = 1;
% % % % %                 end
% % % % %                 if ref == 1  %hay que cambiar esto para que cree la referencia en la media de los puntos
% % % % %                     hline = refline([0 X(i)]); %revisar donde va a dibujar la linea
% % % % %                     hline.Color = 'g';
% % % % %                     hline.LineWidth = 3;
% % % % %                     delete(hline);
% % % % %                 end
            end
%        Y=(P3D(nf,sb(1)+4,1)+P3D(nf,sb(2)+4,1))/2-Ymed(1);
%        Z=(P3D(nf,sb(1)+8,1)+P3D(nf,sb(2)+8,1))/2-Zmed(1);
            if T > 20 % cuando pasan 20 segundos, cambia la escala
                minimo = zeros(1,np);
                maximo = zeros(1,np);
                for i = 1:np
                    minimo(i) = min(P3D(nf-100:nf,sb(1)+8,i)-Xmed(i));
                    maximo(i) = max(P3D(nf-100:nf,sb(1)+8,i)-Xmed(i));
                end
                ax_min = min(minimo);
                ax_max = max(maximo);
%                 axis(handles.axes3,[T-15,T+2,ax_min-1,ax_max+1]); % re-escala los ejes
            end
%             for i = 1:np
% %                 addpoints(hx{i},T,X(i));
%             end
        end
   tnf(4,nf)=toc;
    end
profile off
tnf(5,:)=tnf(1,:)+tnf(2,:)+tnf(3,:)+tnf(4,:);
tfotograma=tnf(:,1:nf);
% para las cámaras
%     stop(vid{1});
%     stop(vid{2});


% solución final
    Salida = P3D(1:nf, [sb sb+4 sb+8], :);
    name = strcat(dir_res,f_salida,'.mat');
    tiempo(:) = t(3,1:nf);
    
    for j = 1:np
        for i = 1:3
            pos(i,:,j) = (Salida(:,2*i-1,j)+Salida(:,2*i,j))/2;
            media = mean(pos(i,:,j));
            pos(i,:,j) = pos(i,:,j)-media;
        end
    end
    
    save(name,'Salida','tiempo','tfotograma')

    disp(strcat('fotograma',num2str(nf)));

    %% MIGUEL-Cerrar Video
       close(V_cam1)
       close(V_cam2)
       
       
       name_file=strcat(folder_name,filesep,outFile_cam1(1:end-4),'_frame_times.mat');
       save(name_file,'time_frame')
% figure;
% hold on
% plot(Salida(:,2,1)-Salida(:,1,1),'r');
% plot(Salida(:,4,1)-Salida(:,3,1),'g');
% plot(Salida(:,6,1)-Salida(:,5,1),'b');
% figure;
% plot(t(2,2:nf)-t(1,2:nf));
% figure;
% plot(1./diff(t(2,2:nf)));
end 

% --- Executes on button press in PB_stop.
function PB_stop_Callback(hObject, ~, handles)

    
    global sigo mi_axes vid seguir;
    sigo = 0;
    seguir = 0;
% muestra los vídeos en la interfaz
    for i = 1:2
        axes(mi_axes{i})
        vidRes = vid{i}.VideoResolution;
        nBands = vid{i}.NumberOfBands;
        hImage = image(zeros(vidRes(1),vidRes(2),nBands)); 
        setappdata(hImage,'UpdatePreviewWindowFcn',@gira_ima);
        preview(vid{i}, hImage);  
    end
    % Activar/Desactivar botones
    set(handles.PB_stop,'Visible','off');
    set(handles.PB_start,'Visible','on');
    
    if strcmp(get(findobj(handles.uibuttongroup1,'Value',1),'tag'),'rb3')
        push_ford_Callback(hObject, [], handles)
    elseif strcmp(get(findobj(handles.uibuttongroup1,'Value',1),'tag'),'rb4')
        push_ford_Callback(hObject, [], handles)
    end
end

function txt_np_Callback(~, ~, ~)
end

% --- Executes during object creation, after setting all properties.
function txt_np_CreateFcn(hObject, ~, ~)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

end

function txt_nhc_Callback(~, ~, ~)
end
    
function txt_nhc_CreateFcn(hObject, ~, ~)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% --- Executes on button press in salir.
function salir_Callback(~, ~, ~)

    global vid;
% para las cámaras
    stop(vid{1});
    stop(vid{2});
 
% resetea la tarjeta capturadora (por si acaso)
    imaqreset;
% cierra la GUI
    close();
    clear;
    
end

% % --- Executes on button press in togglebutton1.
% function togglebutton1_Callback(~, ~, ~)
% end

% --- Executes on button press in limites.
function pushbutton6_Callback(~, ~, ~)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global lim lim2 referencia;
    %axes(mi_axes{3}); 
    [x1,y1] = ginput(1); %limite superior
    %plot(x1,y1,'+g');
    [x2,y2] = ginput(1); %limite inferior
    %plot(x2,y2,'+g');
    lim = [y1, y2];
    lim2 = [x1, x2]; %para crear cuadrado
    h1 = refline([0 y1]); %dibuja lineas de referencia desde el eje Y hasta
    h1.Color = 'b';         %que el plot crece
    h1.LineWidth = 1;
    h2 = refline([0 y2]); 
    h2.Color = 'g';
    h2.LineWidth = 1;
    referencia = 1;
end

function errores(valor)

    switch valor
        case 1
            errordlg('Faltan datos de calibrado: pulse CALIBRAR','Error datos calibración');
        case 2
            errordlg('Error al borrar los datos','Error límites');
        case 3
            errordlg('File not found','File Error');
        otherwise
            errordlg('File not found','File Error');
    end
    
end

% --- Executes on button press in borrar limites.
function pushbutton7_Callback(~, ~, ~)
    global lim;
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    try
        lim = [0, 0];
        uiwait(msgbox('Se han borrado los límites','Límites'));
    catch
        errores(2);
    end
    
end


% --- Executes on button press in calibrado nuevo.
function pushbutton8_Callback(~, ~, ~)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global vid Dif sb;
% saca dos imágenes de las cámaras
    for i=1:2
        I = double(getsnapshot(vid{i}));
%         IM{I}=I;
        IM{i} = imrotate(I,180);
    end
% % % % %         I = double(getsnapshot(vid{1}));
% % % % %         IM{1} = imrotate(I,90);
% % % % %         I = double(getsnapshot(vid{2}));
% % % % %         IM{2} = imrotate(I,-90);
% % % % %         
% Mirar como saca las matrices p1 y p2 en el calibrado normal
    p1 = calibrado_nuevo(IM{1}); 
    p2 = calibrado_nuevo(IM{2});

%save('centros.mat','p1','p2');
%load('centros.mat');

% la función calibrado.m devuelve las matrices de calibrado (MP1 y MP2), y
% las diferencias entre las posiciones reales de los seis puntos y las
% calculadas con el algoritmo
    [~,~, Dif] = calibrado(p1,p2);

% Revisa la máxima diferencia, para aceptar o no el calibrado
    calidad = max(max(max(Dif(sb,:,:))));
    if calidad > 0.3
        mensaje = sprintf('Hay diferencias de hasta %0.1f en la posición de los puntos. No se acepta el calibrado',calidad);
    else
        mensaje = sprintf('Diferencias máximas de %0.2f en la posición de los puntos. Calibrado aceptado',calidad);
        save('calibrado.mat','MP1','MP2')
    end
%set(handles.text1,'String',mensaje);
    uiwait(msgbox(mensaje,'Resultados del calibrado','modal'));
%msgbox(mensaje,'Resultados del calibrado');
    close(1);
    close(2);
end


% --- Executes on button press in grafica paciente.
function pushbutton9_Callback(~, ~, ~)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global lim lim2 ref;
    x1 = lim2(1);
    x2 = lim(1)-lim(2);
    y1 = lim(1);
    y2 = lim(2);
    x = [x1, x2, x2, x1, x1];
    y = [y1, y1, y2, y2, y1];
    plot(x, y, 'b-', 'LineWidth', 3);
    hold on;
    ref = 1;

end



function txt_rep_Callback(hObject, eventdata, handles)
% hObject    handle to txt_rep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_rep as text
%        str2double(get(hObject,'String')) returns contents of txt_rep as a double
end

% --- Executes during object creation, after setting all properties.
function txt_rep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_rep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in push_back.
function push_back_Callback(hObject, eventdata, handles)
% hObject    handle to push_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.push_back,'Value',0)
num = str2num(get(handles.txt_rep,'String'));

if num > 0
num = num - 1;
set(handles.txt_rep,'String',num2str(num));
end
end
% --- Executes on button press in push_ford.
function push_ford_Callback(hObject, eventdata, handles)
% hObject    handle to push_ford (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.push_ford,'Value',0)
num = str2num(get(handles.txt_rep,'String'));
num = num + 1;
set(handles.txt_rep,'String',num2str(num));

end


% --- Executes on button press in rb3.
function rb3_Callback(hObject, eventdata, handles)
% hObject    handle to rb3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb3
set(handles.txt_rep,'String','1');

end


% --- Executes on button press in rb4.
function rb4_Callback(hObject, eventdata, handles)
% hObject    handle to rb4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb4
set(handles.txt_rep,'String','1');
end


% --------------------------------------------------------------------

% --------------------------------------------------------------------

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PB_start_Callback(hObject, [], handles)
end

% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PB_stop_Callback(hObject,[], handles)
end


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end
