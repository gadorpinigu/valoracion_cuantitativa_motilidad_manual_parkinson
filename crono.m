function t = crono()
    c = clock;
%     anyo = c(1,1);
%     mes = c(1,2);
%     dia = c(1,3);
    hora = c(1,4);
    min = c(1,5);
    seg = c(1,6);
%     t = hora*3600+minu*60+seg-t0;
    t = hora*3600+min*60+seg;
end