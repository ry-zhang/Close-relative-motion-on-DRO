function xyLim(x,y,x_ratio,y2x_ratio)
x_max = max(x);
x_min = min(x);
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(y)+min(y))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 