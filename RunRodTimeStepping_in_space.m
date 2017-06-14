x0 = 0;
y0 = 0.5;
theta_0 = pi/6;
t0 = 0;
q0=[x0;y0;theta_0];

h = 1/100;
n = 400;
z = zeros(11, n);

t = t0;
z(1:3,1)=q0;
for i = 1:1:n-1
  z(:,i+1) = RodTimeStepping3_2(z(1:3,i), t, h);
  t = t+h;
end

%% plot
l = 0.5;
for i=1:1:n
  xc = z(1,i);
  yc = z(2,i);
  theta = z(3,i);
  x = [xc - cos(theta)*l/2; xc + cos(theta)*l/2];
  y = [yc - sin(theta)*l/2; yc + sin(theta)*l/2];
  plot(x, y, '-ro')
  hold on
  % xf = i*h;
  xf = 0;
  plot([-0.5, 0.5], [xf, xf])
  hold off
  axis([-0.5 0.5 0 2.5])
  axis equal
  grid on
  drawnow
  frame = getframe(1);
  im{i} = frame2im(frame);
  % pause(2*h)
end

%% save to gif
filename = 'testAnimated.gif'; % Specify the output file name
for idx = 1:n
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.02);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
    end
end