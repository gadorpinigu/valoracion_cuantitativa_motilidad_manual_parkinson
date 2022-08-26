  function F = calibrado_nuevo(A)

    %D1 = imread(A);
    % D1 = imread('Cam2.png');
    D2 = rgb2gray(uint8(A));
    D = imbinarize(D2, 0.5);
    % 
    % 
    % [centersDark, radiiDark] = imfindcircles(B,[6 38],'ObjectPolarity','dark');
    [centersDark, radiiDark] = imfindcircles(D,[6 40],'ObjectPolarity','dark');
    % 
    % imshow(B);
    imshow(D2);
    j = 1;
    for i = 1:length(radiiDark)
        if radiiDark(i)<10
            radio(j,1) = radiiDark(i);
            centro(j,1) = centersDark(i,1);
            centro(j,2) = centersDark(i,2);
            j = j+1;
        end
    end
    viscircles(centro, radio,'Color','b');
    centro2 = centro;
    matriz_A = zeros(16,4);
    for i = 1:length(radio)
        for k = 1:length(radio)
            if i == length(radio) && k == length(radio)
                break;
            elseif k == i
                k = k+1;
            end
            comprobar_y1 = centro2(k,2)-10;
            comprobar_y2 = centro2(k,2)+10;
            comprobar_x1 = centro2(k,1)-40;
            comprobar_x2 = centro2(k,1)+40;
            comprobar_d = centro2(i,1)-centro(k,1);
            if centro2(i,2) >= comprobar_y1 && centro2(i,2) <= comprobar_y2...
                    && centro2(i,1) >= comprobar_x1 && centro2(i,1) <= comprobar_x2...
                    && abs(comprobar_d) > 10
                matriz_A(i,1) = centro2(i,1);
                matriz_A(i,2) = centro2(i,2);
                matriz_A(i,3) = centro2(k,1);
                matriz_A(i,4) = centro2(k,2);
                centro2(i,1) = 0;
                centro2(i,2) = 0;
            end
        end
    end
    hold on;
    plot(matriz_A(:,1), matriz_A(:,2), 'r+');
    plot(matriz_A(:,3), matriz_A(:,4), 'g*');

    radio = 3;
    k = 1;
    C = zeros(16,2);
    for i = 1:size(matriz_A,1)
        if sum(matriz_A(i,:)) ~= 0
            media_y = round((matriz_A(i,2) + matriz_A(i,4))/2);
            media_x = round((matriz_A(i,1) + matriz_A(i,3))/2);
            %al hacer coincidencia con la foto las coordenadas están al revés
    %         matriz_B = B((media_y-radio:media_y+radio),(media_x-radio:media_x+radio));
            matriz_B = D((media_y-radio:media_y+radio),(media_x-radio:media_x+radio));
            centro_x = sum(matriz_B,1);
            centro_y = sum(matriz_B,2);
            for j = 1:length(centro_x)
                if centro_x(j) < 7
                    centro3_x(k) = j;
                    k = k+1;
                end
            end
            k = 1;
            for j = 1:length(centro_y)
                if centro_y(j) < 7
                    centro3_y(k) = j;
                    k = k+1;
                end
            end
            k = 1;
            centro4 = [round(length(centro3_y)/2), round(length(centro3_x)/2)];
    %         matriz_C = B((media_y-radio+centro3_y(1,1)-1:media_y+radio-length(centro3_y)+1)...
    %            ,(media_x-radio+centro3_x(1,1)-1:media_x+radio-length(centro3_x)-1));
            punto = [media_y-radio+centro3_y(1,1)-2+centro4(1,1),...
                media_x-radio+centro3_x(1,1)-2+centro4(1,2)];
            plot(punto(1,2), punto(1,1), 'bx');
            clear centro3_y centro3_x;
            C(i,1) = punto(1,2);
            C(i,2) = punto(1,1);
        end
    end
    C(C==0)=[];
    F = zeros(6,2);
    for i = 1:length(C)/2
        F(i,1) = C(i);
        F(i,2) = C(i+6);
    end
end