clear 
clc
close all

r = 0.05;
x0 = 0;
y0 = r;
theta0 = 0;
xl0 = x0 - r-0.02;
xr0 = x0 + r+0.02;
yg0 = r;
q0=[x0;y0;theta0;xl0;xr0;yg0];
qa_dot_desired_0 = [0.1;-0.1;0.1];
phi_n0 = zeros(3,1);
phi_n0(1) = x0 - xl0 - r; 
phi_n0(2) = xr0 - x0 - r;
phi_n0(3) = y0 - r;

h = 1/100;
n = 60;
nu = 3;
na = 3;
nc = 3;
z_MIQP = zeros(30, n-1); 
delta_q_u = zeros(nu,n-1);
q = zeros(nu+na,n);
penetration = zeros(nc, n);
penetration(:,1) = phi_n0;
z_MIQP(1:nu+na,1)=q0;
q(:,1) = q0;

is_MIQP = true;
for i = 1:1:n-1
  if i == 20
    disp(i)
  end
  disp(i)
  [z_MIQP(:,i), delta_q_u(:,i), penetration(:,i+1)] = ...
  ParallelGripperPickUp_MIQP(q(1:nu,i), q(nu+1:nu+na,i), desired_qa_dot(i,is_MIQP), h, penetration(:,i), is_MIQP);
  q(:,i+1) = q(:,i) + [delta_q_u(:,i);z_MIQP(nu+1:nu+na, i)];
end
%% animation
close all
steps_per_frame = 1;
im = cell(n/steps_per_frame,1);
for i=1:steps_per_frame:n  
  xc = q(1,i);
  yc = q(2,i);
  rectangle('Position',[xc-r, yc-r, 2*r, 2*r],'Curvature',[1 1],'FaceColor',[0 .5 .5])
  hold on
  xl = q(4,i);
  xr = q(5,i);
  yg = q(6,i);
  plot([-1, 1], [0, 0]) % ground
  plot([xl, xl], [yg-r, yg+2*r]);
  plot([xr, xr], [yg-r, yg+2*r])
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
filename = 'ParallelGripperPickUp_without_QP.gif'; % Specify the output file name
for idx = 1:size(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.02);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
    end
end

%% plots
close all
width = 500;
height = 300;
font_size = 20;
qa_dot_desired = zeros(3,n-1);
for i=1:n-1
  qa_dot_desired(:,i) = desired_qa_dot(i, is_MIQP);
end
% plot 1: left finger velocity command vs. time and actual left finger velocity vs. time
h1 = figure(1);
stairs(1:n-1, qa_dot_desired(1,:)*h, 'LineWidth', 3)
hold on
stairs(1:n-1, z_MIQP(nu+1,:), '--','LineWidth', 3)
% plot 1 formatting
grid on
axis([1, n-1, -1.2e-3, 1.2e-3])
xlabel 'l (time step)'
ylabel '\Delta x_l (meter)'
l=legend('$$\mathrm{commanded} \: \Delta x_l$$', '$$\mathrm{actual} \: \Delta x_l$$',  'Location', 'southwest');
set(l, 'Interpreter','latex', 'FontSize', font_size);
set(h1,'units','points','position',[10,10,width,height])
ax = gca;
set(ax, 'FontSize', font_size)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% plot 2: right finger velocity command vs. time and actual right finger velocity vs. time
h2 = figure(2);
stairs(1:n-1, qa_dot_desired(2,:)*h, 'LineWidth', 3)
hold on
stairs(1:n-1, z_MIQP(nu+2,:), '--','LineWidth', 3)
grid on
axis([1, n-1, -1.2e-3, 1.2e-3])
xlabel 'l (time step)'
ylabel '\Delta x_r (meter)'
l=legend('$$\mathrm{commanded} \: \Delta x_r$$', '$$\mathrm{actual} \: \Delta x_r$$',  'Location', 'northwest');
set(l, 'Interpreter','latex', 'FontSize', font_size);
set(h2,'units','points','position',[510,10,width,height])
ax = gca;
set(ax, 'FontSize', font_size)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% plot 3: vertical velocity of both fingers vs. time, commanded and actual
h3 = figure(3);
stairs(1:n-1, qa_dot_desired(3,:)*h, 'LineWidth', 3)
hold on
stairs(1:n-1, z_MIQP(nu+3,:), '--','LineWidth', 3)
grid on
axis([1, n-1, -1.2e-3, 1.2e-3])
xlabel 'l (time step)'
ylabel '\Delta y_g (meters)'
l=legend('$$\mathrm{commanded} \: \Delta y_g$$', '$$\mathrm{actual} \: \Delta y_g$$',  'Location', 'southwest');
set(l, 'Interpreter','latex', 'FontSize', font_size);
set(h3,'units','points','position',[1010,10,width,height])
ax = gca;
set(ax, 'FontSize', font_size)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% plot 4: CG of cylinder x and y position vs. time
h4 = figure(4);
stairs(1:n, q(1,:), '-r', 'LineWidth', 2)
hold on
stairs(1:n, q(2,:), '-g', 'LineWidth', 2)
stairs(1:n, r*q(3,:), '-m','LineWidth', 2)
stairs(1:n, penetration(1,:), '-b', 'LineWidth', 2)
grid on
axis([1, n, -0.01, 0.07])
xlabel 'l (time step)'
ylabel '(meter)'

l = legend('$$x_c$$', '$$y_c$$','$$\theta r$$','$$\bar{\phi}_{n_1}$$', 'Location', 'best');
set(l, 'Interpreter','latex', 'FontSize', font_size);
set(h4,'units','points','position',[10,400,width,height+50])
ax = gca;
set(ax, 'FontSize', font_size)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% plot 5: contact forces at gripper contacts vs. time
h5 = figure(5);
yyaxis left
stairs(2:n, z_MIQP(6,:), 'LineWidth', 2)
hold on
yyaxis right
ylabel '(meter)'
stairs(2:n, penetration(1,2:end), 'LineWidth', 2)
refline(0,0)
grid on
yyaxis left
axis([2, n, -1, 11])
xlabel 'l (time step)'
ylabel '(Newton)'
l = legend('$$\lambda_{n_{1/2}}$$', '$$\bar{\phi}_{n_{1/2}}$$', 'Location', 'best');
set(l, 'Interpreter','latex', 'FontSize', font_size);
set(h5,'units','points','position',[10,400,width,height+50])
ax = gca;
set(ax, 'FontSize', font_size)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% plot 6: contact forces at ground contacts vs. time
h6 = figure(6);
yyaxis left
stairs(2:n, z_MIQP(8,:), 'LineWidth', 2)
hold on
yyaxis right
ylabel '(meter)'
stairs(2:n, penetration(3,2:end), 'LineWidth', 2)
refline(0,0)
grid on
yyaxis left
axis([2, n, -1, 10])
xlabel 'l (time step)'
ylabel '(Newton)'
l = legend('$$\lambda_{n_{3}}$$', '$$\bar{\phi}_{n_{3}}$$', 'Location', 'best');
set(l, 'Interpreter','latex', 'FontSize', font_size);
set(h6,'units','points','position',[10,400,width,height+50])
ax = gca;
set(ax, 'FontSize', font_size)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%% save plots
if(is_MIQP)
  print(h1,'~\Dropbox (MIT)\QuasistaticSimulatorPaper\delta_xl_vs_time','-depsc')
  print(h2,'~\Dropbox (MIT)\QuasistaticSimulatorPaper\delta_xr_vs_time','-depsc')
  print(h3,'~\Dropbox (MIT)\QuasistaticSimulatorPaper\delta_yg_vs_time','-depsc')
  print(h4,'~\Dropbox (MIT)\QuasistaticSimulatorPaper\xy_object_vs_time','-depsc')
  print(h5,'~\Dropbox (MIT)\QuasistaticSimulatorPaper\gripper_normal_force_vs_time','-depsc')
  print(h6,'~\Dropbox (MIT)\QuasistaticSimulatorPaper\ground_normal_force_vs_time','-depsc')
end

if(~is_MIQP) 
  print(h4,'~\Dropbox (MIT)\QuasistaticSimulatorPaper\xy_object_vs_time_lcp','-depsc')
  print(h5,'~\Dropbox (MIT)\QuasistaticSimulatorPaper\gripper_normal_force_vs_time_lcp','-depsc')
  print(h6,'~\Dropbox (MIT)\QuasistaticSimulatorPaper\ground_normal_force_vs_time_lcp','-depsc')
end




























