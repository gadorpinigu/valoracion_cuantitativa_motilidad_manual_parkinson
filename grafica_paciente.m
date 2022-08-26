figure
get(0,'MonitorPositions');
scrsz = get(0,'ScreenSize');
imshow(im, 'parent', gca);

x1=0;
x2=1;
y1=0;
y2=1;
x = [x1, x2, x2, x1, x1];
y = [y1, y1, y2, y2, y1];
plot(x, y, 'b-', 'LineWidth', 3);
hold on;

xlim([-1, 2]);
ylim([-1, 2]);

hline = refline([0 y1]);
hline.Color = 'g';
hline.LineWidth = 3;
%h = findobj('Color','g');
for i = 1:1000
    hline = refline([0 i]);
    hline.Color = 'g';
    pause(0.05);
    delete(hline);
end