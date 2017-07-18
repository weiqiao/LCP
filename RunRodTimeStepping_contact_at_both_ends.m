clear

x0 = 0;
y0 = 0.5;
z0 = 0;
theta0 = pi/6;
t0 = 0;
qa0 = 0.25;
q0=[x0;y0;z0;theta0;0;qa0];


h = 5/100;
n = 150;
nu = 5;
na = 1;
n1 = nu+na;
z_static = zeros(25, n);
z_MIQP = zeros(46, n-1); 
delta_q_u = zeros(nu,n-1);
q = zeros(n1,n);

t = t0;
z_static(1:nu,1)=q0(1:nu);
z_MIQP(1:6,1) = q0; % qa(t=0) = 0.25;
q(:,1) = q0;

for i = 1:1:n-1
  if i==26
    disp('hello')
  end
  disp(i)
  z_static(:,i+1) = RodTimeStepping2(z_static(1:5,i), t, h);
  [z_MIQP(:,i), delta_q_u(:,i)] = RodTimeStepping2_MIQP(q(1:nu,i), q(nu+1:nu+na,i), 0.1, h);
  q(:,i+1) = q(:,i) + [delta_q_u(:,i);z_MIQP(nu+1:nu+na, i)];
  t = t+h;
end

%% plot
close all

l = 0.5;
steps_per_frame = 1;
im = cell(n/steps_per_frame,1);
for i=1:steps_per_frame:n
  xc = z_static(1,i);
  yc = z_static(2,i);
  theta = z_static(4,i);
  x = [xc - cos(theta)*l/2; xc + cos(theta)*l/2];
  y = [yc - sin(theta)*l/2; yc + sin(theta)*l/2];
  plot(x, y, '-ro')
  hold on
  
  xc = q(1,i);
  yc = q(2,i);
  theta = q(4,i);
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