clear 

x0 = 0.1;
y0 = 0;
t0 = 0;
q0=[x0;y0;0.2;-0.1];
qa_dot_desired = [-0.1;0.1];

h = 1/100;
n = 150;
z_MIQP = zeros(30, n); 

t = t0;
z_MIQP(1:4,1)=q0;

for i = 1:1:n-1
  disp(i)
  z_MIQP(:,i+1) = BallSqueezAndPush_MILP(z_MIQP(1:2,i), z_MIQP(3:4,i), qa_dot_desired, h);
end

%% plot
close all
r = 0.05;
steps_per_frame = 1;
im = cell(n/steps_per_frame,1);
for i=1:steps_per_frame:n  
  xc = z_MIQP(1,i);
  yc = z_MIQP(2,i);
  rectangle('Position',[xc-r, yc-r, 2*r, 2*r],'Curvature',[1 1],'FaceColor',[0 .5 .5])
  hold on
  xf = z_MIQP(3,i);
  yf = z_MIQP(4,i);
  plot([xf, xf], [-1, 1])
  plot([-1,1], [yf, yf]);
  plot([0, 0], [-1, 1])
  axis([-0.2 0.2 -0.3 0.3])
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