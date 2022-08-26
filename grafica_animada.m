function grafica_animada

h = animatedline;
x = linspace(0,40*pi,1000);
y = sin(x);
axis([0,4*pi,-1,1])
    for k = 1:length(x)
    addpoints(h,x(k),y(k));
        if x(k) > 3.9*pi
            axis([x(k)-3.9*pi,x(k)+0.1*pi,-1,1])
        end
    drawnow
    pause(0.1)
    end

end