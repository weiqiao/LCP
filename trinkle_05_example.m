q0 = [0.6, 0, 0]';
t0 = 0;
h = 0.01;

n = 200;
q = zeros(3, n);
q(:,1) = q0;
t = t0;
for i = 1:1:n-1
  q(:,i+1) = TimeSteppingLCP(q(:,i), t, h);
  t = t+h;
end

%% plot
for i=1:1:n
  plot(q(1,i), q(2,i), 'ro')
  hold on
  xf = 0.5 + 0.4*sin(t0+(i-1)*h);
  plot([xf, xf], [0, 1])
  hold off
  axis([0 1 0 1])
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