l = 0.5;
x = z_quasi_dynamic(1,1:304);
y = z_quasi_dynamic(2,1:304);
theta = z_quasi_dynamic(4,1:304);
xB = x - l/2*cos(theta);
yB = y - l/2*sin(theta);
xB_dot = (xB(2:end)-xB(1:end-1))/h; 
yB_dot = (yB(2:end)-yB(1:end-1))/h; 