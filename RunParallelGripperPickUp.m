clear 
clc

r = 0.05;
x0 = 0;
y0 = r;
xl0 = -r-0.02;
xr0 = r+0.02;
yg0 = -r;
q0=[x0;y0;xl0;xr0;yg0];
qa_dot_desired = [0.1;-0.1;0.1];

h = 1/100;
n = 150;
nu = 2;
na = 3;
z_MIQP = zeros(29, n-1); 
delta_q_u = zeros(nu,n-1);
q = zeros(nu+na,n);
z_MIQP(1:nu+na,1)=q0;
q(:,1) = q0;

for i = 1:1:n-1
  if i == 51
    disp(i)
  end
  [z_MIQP(:,i), delta_q_u(:,i)] = ParallelGripperPickUp_MIQP(q(1:nu,i), q(nu+1:nu+na,i), qa_dot_desired, h);
  q(:,i+1) = q(:,i) + [delta_q_u(:,i);z_MIQP(nu+1:nu+na, i)];
end

%% plot
close all
steps_per_frame = 1;
im = cell(n/steps_per_frame,1);
for i=1:steps_per_frame:n  
  xc = q(1,i);
  yc = q(2,i);
  rectangle('Position',[xc-r, yc-r, 2*r, 2*r],'Curvature',[1 1],'FaceColor',[0 .5 .5])
  hold on
  xl = q(3,i);
  xr = q(4,i);
  yg = q(5,i);
  plot([-1, 1], [0, 0]) % ground
  plot([xl, xl], [yg, yg+3*r]);
  plot([xr, xr], [yg, yg+3*r])
  axis([-0.2 0.2 -0.1 0.5])
  axis equal
  grid on
  hold off
  drawnow
  frame = getframe(1);
  im{i} = frame2im(frame);
  clf
  % pause(2*h)
end

%% save to gif
filename = 'BallSqueezeAndPush.gif'; % Specify the output file name
for idx = 1:size(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.02);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
    end
end