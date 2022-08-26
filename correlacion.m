function p = correlacion(A, n_fig)

r = 4; % radio de búsqueda del punto bueno en el entorno del punto aproximado



%--------------------------------------------------------------------------
% 1. LEE LAS IMÁGENES Y HACE LA CORRELACIÓN
%--------------------------------------------------------------------------

%muestra la imagen problema
figure(n_fig);
imshow(uint8(A));


%--------------------------------------------------------------------------
% 3. BUSCA LOS PUNTOS 
%--------------------------------------------------------------------------

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


% numero los puntos en la figura
    figure(n_fig);
    for i = 1:6
        text(p(i,1),p(i,2)+2,int2str(i));
     end
