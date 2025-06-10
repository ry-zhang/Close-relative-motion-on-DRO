function set_axis(xlu)
xlu = [min(xlu(:,1:3))-0.2,max(xlu(:,4:6))+0.2];
xlim([xlu(1) xlu(4)]);
ylim([xlu(2) xlu(5)]);

