x0 = 0;
y0 = 0.5;
z0 = 0;
theta0 = pi/6;
t0 = 0;
q0=[x0;y0;z0;theta0];

h = 5/100;
n = 200;
z_static = zeros(25, n);

t = t0;
z_static(1:4,1)=q0;
z_MILP = z_static; % quasi-dynamic
for i = 1:1:n-1
  if i==304
    disp('hello')
  end
  z_static(:,i+1) = RodTimeStepping2(z_static(1:4,i), t, h);
  z_MILP(:,i+1) = RodTimeStepping2_MILP(z_MILP(1:4,i), h, yG(t), yG_dot(t));
  t = t+h;
end

%% plot
l = 0.5;
im = cell(n,1);
for i=1:1:n
  xc = z_static(1,i);
  yc = z_static(2,i);
  theta = z_static(4,i);
  x = [xc - cos(theta)*l/2; xc + cos(theta)*l/2];
  y = [yc - sin(theta)*l/2; yc + sin(theta)*l/2];
  plot(x, y, '-ro')
  hold on
  
  xc = z_MILP(1,i);
  yc = z_MILP(2,i);
  theta = z_MILP(4,i);
  x = [xc - cos(theta)*l/2; xc + cos(theta)*l/2];
  y = [yc - sin(theta)*l/2; yc + sin(theta)*l/2];
  plot(x, y, '--go')
  
  %%%%%%%% plot dynamic trajectory from dynamic rod2d simulation
%   xc = q(1,i);
%   yc = q(2,i);
%   theta = q(4,i);
%   x = [xc - cos(theta)*l/2; xc + cos(theta)*l/2];
%   y = [yc - sin(theta)*l/2; yc + sin(theta)*l/2];
%   plot(x, y, '--bo')
  %%%%%%%%%%%%%%
  
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
  legend('LCP', 'MIQP')
  drawnow
  frame = getframe(1);
  im{i} = frame2im(frame);
  % pause(2*h)
end

%% save to gif
filename = 'LCP_as_MIQP_big_M_cost_on_force_balance.gif'; % Specify the output file name
for idx = 1:size(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.02);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
    end
end