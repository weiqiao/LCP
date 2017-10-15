function [z, delta_q_u] = RodTimeStepping2_MIQP(qu, qa, qa_dot_d, h, z_start)
% assuming friction force between table and rod exists at both ends of the
% rod. 
% qa = yG(t) \in R^1
% qa_dot = d{yG(t)}/dt;

qu_l = qu;
qa_l = qa;
yc = qu_l(2);
zc = qu_l(3);
theta = qu_l(4);

% defining constants
r = 0; % m
l = 0.5; % m
mu = 0.6;
f = [0;0;-9.8;0];
f = [f;0]; 

nu = 5; % number of unactuated states
na = 1; % number of actuated states
n = nu + na;
nc = 4; % number of contacts
nf = [2, 2, 4, 4];
nd = sum(nf); % no. of columns in Wf

Wf1 = @(theta) [1 0 0 l/2*sin(theta);-1 0  0 -l/2*sin(theta)]';
Wf2 = @(theta) [1 0 0 -l/2*sin(theta);-1 0  0 l/2*sin(theta)]';
Wf3 = @(theta) [1 0 0 l/2*sin(theta);0 1 0 -l/2*cos(theta);
               -1 0 0 -l/2*sin(theta);0 -1 0 l/2*cos(theta)]';
Wf4 = @(theta) [1 0 0 -l/2*sin(theta);0 1 0 l/2*cos(theta);
               -1 0 0 l/2*sin(theta);0 -1 0 -l/2*cos(theta)]';
Wf = [Wf1(theta), Wf2(theta), Wf3(theta), Wf4(theta)];
Wf = [Wf;zeros(1,nd)];

n1 = @(theta) [0 1 0 -l/2*cos(theta)]';
n2 = @(theta) [0 1 0 l/2*cos(theta)]';
n3 = [0 0 1 0]';
n4 = [0 0 1 0]';
Wn = [n1(theta), n2(theta), n3, n4];
Wn = [Wn;[0 0 1 -1]];

Ja = zeros(nc+nd,na);
Ja(1:2) = -1;

E = zeros(nc, nd);
first = 1;
for i = 1:1:nc
  E(i, first:first+nf(i)-1) = 1;
  first = first + nf(i);
end
E = E';
U = diag([mu, mu, mu, mu]);
B = [zeros(nu,n) Wn Wf zeros(nu,nc);
     Wn' Ja(1:nc) zeros(nc, nc + nd + nc);
     Wf' Ja(nc+1:nc+nd) zeros(nd, nc + nd) E;
     zeros(nc, n) U  -E'  zeros(nc,nc)];
%%
% yG = @(t) t;
phi_n = zeros(4,1);
phi_n(1) = yc - l/2*sin(theta) - r - qa_l; 
phi_n(2) = yc + l/2*sin(theta) - r - qa_l;
phi_n(3) = zc; 
phi_n(4) = zc;
bn = phi_n; % 4*1
bf = zeros(nd,1); % 4*1
b = [f; bn; bf; zeros(nc,1)];

%% call solver
M = 50;
z = Mixed_LCP_As_MILP3(B,b,nu,na,h,qa_dot_d, M, z_start);
% delta_q_u = z(1:nu);

% QP to minimize delta_qu
epsilon = 1e-6;
delta_q_u = MinDeltaQu_QP(z, Wn, Wf, Ja, phi_n, nu, na, nc, nd, nf, epsilon);