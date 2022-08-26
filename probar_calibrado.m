
vid{1} = TTK_videoinput(1);
vid{2} = TTK_videoinput(2);
start(vid{1});
start(vid{2});
    
for i=1:2
    I = double(getsnapshot(vid{i}));
    IM{i} = imrotate(I,90);
end
p1 = correlacion(IM{1},'Patron_D.png',1);
p2 = correlacion(IM{2},'Patron_I.png',2);
[MP1,MP2, Dif] = calibrado(p1,p2);
stop(vid{1});
stop(vid{2});
imaqreset;

