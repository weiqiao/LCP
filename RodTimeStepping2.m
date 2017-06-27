function z = RodTimeStepping2(q, t, h)
% assuming friction force between table and rod exists at both ends of the
% rod. 
ql = q;
xc = ql(1);
yc = ql(2);
zc = ql(3);
theta = ql(4);

% defining constants
r = 0; % m
l = 0.5; % m
mu = 0.6;
f = [0;0;-9.8;0];
f = [f;0]; 

n = 4; % dimension of configuration q
nc = 4; % number of contacts

Wf1 = @(theta) [1 0 0 l/2*sin(theta);-1 0  0 -l/2*sin(theta)]';
Wf2 = @(theta) [1 0 0 -l/2*sin(theta);-1 0  0 l/2*sin(theta)]';
Wf3 = @(theta) [1 0 0 l/2*sin(theta);0 1 0 -l/2*cos(theta);
               -1 0 0 -l/2*sin(theta);0 -1 0 l/2*cos(theta)]';
Wf4 = @(theta) [1 0 0 -l/2*sin(theta);0 1 0 l/2*cos(theta);
               -1 0 0 l/2*sin(theta);0 -1 0 -l/2*cos(theta)]';
Wf = [Wf1(theta), Wf2(theta), Wf3(theta), Wf4(theta)];
nd = size(Wf,2); % no. of columns in Wf
Wf = [Wf;zeros(1,nd)];

n1 = @(theta) [0 1 0 -l/2*cos(theta)]';
n2 = @(theta) [0 1 0 l/2*cos(theta)]';
n3 = [0 0 1 0]';
n4 = [0 0 1 0]';
Wn = [n1(theta), n2(theta), n3, n4];
Wn = [Wn;[0 0 1 -1]];

E = zeros(nc, nd);
E(1,1:2) = 1;
E(2,3:4) = 1;
E(3,5:8) = 1;
E(4,9:12) = 1;
E = E';
U = diag([mu, mu, mu, mu]);
B = [zeros(n+1,n+1) Wn Wf zeros(n+1,nc);
     %zeros(1,n+1) [1 -1 0 0], zeros(1,nd+nc);
     Wn' zeros(nc, nc + nd + nc);
     Wf' zeros(nd, nc + nd) E;
     zeros(nc, n+1) U  -E'  zeros(nc,nc)];
%%
% yG = @(t) t;
phi_n = zeros(4,1);
phi_n(1) = yc - l/2*sin(theta) -r - yG(t); 
phi_n(2) = yc + l/2*sin(theta) -r - yG(t);
phi_n(3) = zc; 
phi_n(4) = zc;
bn = phi_n + h*[-1;-1;0;0]*0.1 - Wn'*ql; % 4*1
bf = h*zeros(nd,1) - Wf'*ql; % 4*1
b = [f; bn; bf; zeros(nc,1)];

N = length(b);
lb = zeros(N,1);
lb(1:n+1) = -Inf;
ub = Inf*ones(N,1);
%%
[z,mu] = pathlcp(B,b,lb,ub);