x0 = 0;
y0 = 0.5;
z0 = 0;
theta0 = pi/6;
t0 = 0;
q0=[x0;y0;z0;theta0];

h = 1/100;
n = 1000;
z_static = zeros(20, n);

t = t0;
z_static(1:4,1)=q0;

profile on
for i = 1:1:n-1
  z_static(:,i+1) = RodTimeStepping5(z_static(1:4,i), t, h);
  t = t+h;
end
profile viewer

%% plot
l = 0.5;
im = cell(n/10,1);
for i=1:10:n
  xc = z_static(1,i);
  yc = z_static(2,i);
  theta = z_static(4,i);
  x = [xc - cos(theta)*l/2; xc + cos(theta)*l/2];
  y = [yc - sin(theta)*l/2; yc + sin(theta)*l/2];
  plot(x, y, '-ro')
  hold on
  
  % original position
%   xl0 = [x0 - cos(theta0)*l/2; x0 + cos(theta0)*l/2];
%   yl0 = [y0 - sin(theta0)*l/2; y0 + sin(theta0)*l/2];
%   plot(xl0, yl0, '--ko')
  
  yf = yG(i*h);
  plot([-0.5, 0.5], [yf, yf])
  hold off
  axis([-0.5 0.5 0 2.5])
  axis equal
  grid on
  drawnow
  frame = getframe(1);
  im{floor(i/10)+1} = frame2im(frame);
  % pause(2*h)
end

%% save to gif
filename = 'quasi_static_vs_quasi_dynamic_vs_dynamic.gif'; % Specify the output file name
for idx = 1:size(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.02);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
    end
end