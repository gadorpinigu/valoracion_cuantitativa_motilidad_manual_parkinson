function p = correlacion(A, n_fig)
% A es la imagen en la que buscar
% 'patron' es el nombre de la imagen patrón a usar
% n_fig es la figura en la que tiene que mostrar las imágenes
% busca las zonas de 'imagen' que tienen mayor correlación con la imagen
% patrón, por medio de una función de correlación 2D.
% Empieza con un nivel de correlación exigente (un threshold de 50%) y va
% bajando el nivel hasta que encuentra seis o más puntos. Si encuentra
% seis, considera el problema resuelto. Si encuentra más de seis (o sea, ha
% detectado correlaciones falsas) se queda con las encontradas en la
% anterior iteración, y pide al usuario que señale las que faltan a mano.

r = 4; % radio de búsqueda del punto bueno en el entorno del punto aproximado



%--------------------------------------------------------------------------
% 1. LEE LAS IMÁGENES Y HACE LA CORRELACIÓN
%--------------------------------------------------------------------------

%muestra la imagen problema
figure(n_fig);
imshow(uint8(A));

% lee la imagen patrón
% B = double(imread(patron));

% En ambas imágenes, se queda solo con el primer canal, convierte a double y normaliza
% ANri = double(A(:,:,1));
% AN = AN-mean(AN(:));
% BN = double(B(:,:,1));
% BN = BN-mean(BN(:));
% AN = A(:,:,1)-mean2(A(:,:,1));
% BN = B(:,:,1)-mean2(B(:,:,1));

% Calcula la matriz de correlación 2D entre las dos imágenes
% Z = xcorr2(AN,BN)/max(max(xcorr2(AN,BN)));
% Z = Z/max(max(Z));



%--------------------------------------------------------------------------
% 2. LOCALIZA LOS PUNTOS DE MAYOR CORRELACIÓN
%--------------------------------------------------------------------------

% % Itera para sucesivos thresholds. Va bajando el threshold hasta que
% % encuentra seis o más puntos
% ther = 0.5;
% encontrados = 0;
% prop = [];
% prop_ant = [];
%     while encontrados < 6
%         prop_ant = prop;
%         ther = ther-0.05;
%         % Convierte la matriz de correlación en una imagen binaria, con un threshold del 30%
%         BW = imbinarize(Z,ther);
% %figure(2)
% %imshow(BW)
% 
% % busca regiones en la imagen BW
%         [LabeledIm,~] = bwlabel(BW,4);
%         prop = regionprops(LabeledIm,'Eccentricity','Centroid','Area');
%         encontrados = size(prop,1);
%     end
% % Si, en la última iteración, ha encontrado más de seis puntos, se queda
% % con la anterior iteración.
%     if encontrados > 6
%         prop = prop_ant;
%     end
% 
% % Reactiva la imagen inicial
%     figure(n_fig)
%     hold on
%     %p = size(prop,1);
% % muestra los resultados encontrados
%     for i = 1:length(prop) %itera para cada punto
%         %ecc(i) = prop(i).Eccentricity;
%         % corrijo la posición por la distancia al centro de la imagen patrón
%         p(i,:) = prop(i).Centroid-[size(B,2)/2 size(B,1)/2];
%     
%     % muestro la primera aproximación al punto
%         plot(p(i,1),p(i,2),'xg')
%     % dibuja una circunferencia
%         circ = [p(i,1)-r, p(i,2)-r, 2*r, 2*r];
%         rectangle('Position',circ,'Curvature',[1 1],'EdgeColor','g')
%     
%     % calcula la posición del centroide del punto
%         p(i,:) = contrasteJ(p(i,1),p(i,2),A,r,r);
%     
%     % muestra la posición calculada
%         plot(p(i,1),p(i,2),'xr')
%     end


%--------------------------------------------------------------------------
% 3. BUSCA LOS PUNTOS QUE NO HA DETECTADO AUTOMÁTICAMENTE
%--------------------------------------------------------------------------

% si no ha encontrado seis puntos, pide que se introduzcan manualmente los
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

%--------------------------------------------------------------------------
% 4. ORDENA LOS PUNTOS
%--------------------------------------------------------------------------

% ordeno los puntos para que estén 2   6
%                                    4
%                                    3
%                                  1   5

% % primero los ordeno por la coordenada x
%     [~,k] = sort(p(:,1));
%     p = [p(k,1),p(k,2)];
% % ordeno las pareja por la coordenada y, en orden descendente
%     for i = 1:3
%         n = 2*i-1;
%         [~,k] = sort(p(n:n+1,2),'descend');    
%         p(n:n+1,:) = [p(k+2*(i-1),1),p(k+2*(i-1),2)];
%     end

% numero los puntos en la figura
    figure(n_fig);
    for i = 1:6
        text(p(i,1),p(i,2)+2,int2str(i));
     end
% 
% pause
% %--------------------------------------------------------------------------
% % 4. CAMBIA A COORDENADAS XY NORMALES
% %--------------------------------------------------------------------------
% 
% % Los puntos están en el sistema de coordenadas de las imágenes de matlab,
% % en que la coordenada Y crece hacia abajo, y el origen de coordenadas está
% % en la esquina superior izquierda. Paso a un sistema de coordenadas
% % centrado en el centro de la imagen y con las orientaciones convencionales
% 
% % tamaño de la imagen
% [s_y,s_x,~] = size(A);
% 
% % invertimos las coordenadas y
% p(:,2) = s_y-(p(:,2));
% 
% % trasladamos los ejes de coordenadas al centro de la imagen. 
% p(:,1) = p(:,1)-s_x/2;
% p(:,2) = p(:,2)-s_y/2;
