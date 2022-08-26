function fabrica_patrones% para buscar los patrones

vid{1} = TTK_videoinput(1);
vid{2} = TTK_videoinput(2);
Im{1} = getsnapshot(vid{1});
Im{2} = getsnapshot(vid{2});
imwrite(Im{1},'Cam1.png');
imwrite(Im{2},'Cam2.png');
stop(vid{1});
stop(vid{2});

% carga la imagen 1 y la visualiza
IM = imread('C:\ttkt2\Cam1.png');
IM = imrotate(IM,90);
imshow(IM)
% buscamos los limites de un circulo
IM2 = IM(404:428,185:237);
imshow(IM2)
imwrite(IM2,'C:\ttkt2\Patron_D.png');

% carga la imagen 2 y la visualiza
IM = imread('C:\ttkt2\Cam2.png');
IM = imrotate(IM,90);
imshow(IM)
% buscamos los limites de un circulo
IM2 = IM(422:446,382:434);
imshow(IM2)
imwrite(IM2,'C:\ttkt2\Patron_I.png');

end