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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Presentacion grafica para tests...

% size(i)
% figure(3)
% clf
% subplot(131)
% plot(X2(i,:),F2(i,:),'b')
% hold on
% plot(Y2(:,j),F2(:,j),'r')
% subplot(132)
% image(uint8(F2*255))
% colormap(gray(256))
% subplot(133)
% image(uint8(F*255))
% colormap(gray(256))
% pause 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% En dos etapas: una primera para centrar el punto basado en la posicion
% del mï¿½ximo, y luego una vecindad mï¿½s pequeï¿½a para el centroide.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basado en posicion del maximo

% [i,j]=find(F2==max(F2(:)));
% Puntos=[mean(mean(X2(i,j)))+round(x) mean(mean(Y2(i,j)))+round(y)]

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Basado en ajuste minimos cuadrados a un paraboloide...
% k=find(M==1);
% X3=[X2(k)*0.0+1 X2(k) Y2(k) X2(k).*Y2(k) X2(k).^2 Y2(k).^2]; 
% Z3=F2(k);
% 
% % Ajuste a la ecuacion 
% % A(1) + A(2) x + A(3) y + A(4) xy + A(5) x^2 + A(6) y^2
% A=pinv(X3)*Z3;
% % Y ahora encuentro el maximo del chisme...
% B=[2*A(5) A(4); A(4) 2*A(6)]\[-A(2); -A(3)];
% 
% Puntos=[B(1)+round(x) B(2)+round(y)]
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basado en posicion del centro de gravedad...
% Puntos=[sum(X(:).*F(:))./sum(F(:))+round(x), sum(Y(:).*F(:))./sum(F(:))+round(y)];
Puntos = [sum(X2(:).*F2(:))./sum(F2(:))+round(x), sum(Y2(:).*F2(:))./sum(F2(:))+round(y)];

end
