x0 = 0;
y0 = 0.5;
theta_0 = pi/6;
t0 = 0;
q0=[x0;y0;0;theta_0];
v0=[0;0;0;0];

h = 1/100;
n = 1000;
z = zeros(16, n);
q = zeros(4, n);
v = q;
q(:,1) = q0;
v(:,1) = v0;
t = t0;
nc = 2;
for i = 1:1:n-1
  [z(:,i+1), q(:,i+1), v(:,i+1)] = RodTimeSteppingDynamicOnTable(q(:,i), v(:, i), h, t);
  t = t+h;
end


%% plot
l = 0.5;
for i=1:10:n
  xc = q(1,i);
  yc = q(2,i);
  theta = q(4,i);
  x = [xc - cos(theta)*l/2; xc + cos(theta)*l/2];
  y = [yc - sin(theta)*l/2; yc + sin(theta)*l/2];
  plot(x, y, '-ro')
  hold on
  
  yf = yG(i*h);
  plot([-0.5, 0.5], [yf, yf])
  hold off
  axis([-0.5 0.5 0 2.5])
  axis equal
  grid on
  % legend('quasi-static', 'quasi-dynamic')
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