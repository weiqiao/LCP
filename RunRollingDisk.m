clear

x0 = 0;
y0 = 0.5;
theta0 = 0;
t0 = 0;
q0=[y0;theta0;x0];

h = 5/100;
n = 50;
nu = 2;
na = 1;
n1 = nu+na;
z_MIQP = zeros(11, n-1); 
delta_q_u = zeros(nu,n-1);
q = zeros(n1,n);

t = t0;
z_MIQP(1:n1,1) = q0; % qa(t=0) = 0.25;
q(:,1) = q0;

for i = 1:1:n-1
  if i==26
    disp('hello')
  end
  disp(i)
  [z_MIQP(:,i), delta_q_u(:,i)] = RollingDiskTimeStepping(q(1:nu,i), q(nu+1:nu+na,i), 0.1, h);
  q(:,i+1) = q(:,i) + [delta_q_u(:,i);z_MIQP(nu+1:nu+na, i)];
  t = t+h;
end

%% plot
close all
load drake_solution
q_drake = log(:, 2:151);
l = 0.5;
steps_per_frame = 1;
im = cell(n/steps_per_frame,1);
for i=1:steps_per_frame:n
  xc = q(1,i);
  yc = q(2,i);
  theta = q(4,i);
  x = [xc - cos(theta)*l/2; xc + cos(theta)*l/2];
  y = [yc - sin(theta)*l/2; yc + sin(theta)*l/2];
  plot(x, y, '--go')
  
  xc = q_drake(1,i);
  yc = q_drake(2,i);
  theta = q_drake(4,i);
  x = [xc - cos(theta)*l/2; xc + cos(theta)*l/2];
  y = [yc - sin(theta)*l/2; yc + sin(theta)*l/2];
  plot(x, y, '--bo')
  
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
  legend('LCP', 'MIQP', 'drake MIQP')
  drawnow
  frame = getframe(1);
  % im{floor(i/steps_per_frame)+1} = frame2im(frame);
  im{i} = frame2im(frame);
  % pause(2*h)
end

%% save to gif
filename = 'LCP_vs_MIQP+QP.gif'; % Specify the output file name
for idx = 1:size(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.02);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
    end
end