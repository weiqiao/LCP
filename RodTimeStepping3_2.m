function z = RodTimeStepping3_2(q, t, h)
% assuming no table and no gravity, i.e. table moves up to contact a rod in
% space. 
% "quasi-dynamic"

ql = q;
xc = ql(1);
yc = ql(2);
theta = ql(3);

n = 3; % dimension of configuration q
nc = 2; % number of contacts

% defining constants
r = 0.05; % m
l = 0.5; % m
mu = 0.6;
f = -[0, 9.8, 0]'*20;
M = eye(3);

Wf1 = @(theta) [-1 0 -l/2*sin(theta);1 0 l/2*sin(theta)]';
Wf2 = @(theta) [-1 0 l/2*sin(theta);1 0 -l/2*sin(theta)]';
Wf = [Wf1(theta), Wf2(theta)];

n1 = @(theta) [0 1 -l/2*cos(theta)]';
n2 = @(theta) [0 1 l/2*cos(theta)]';
Wn = [n1(theta), n2(theta)];
   
nd = size(Wf,2); % no. of columns in Wf, =4

E = zeros(nc, nd);
E(1,1:2) = 1;
E(2,3:4) = 1;
E = E';
U = diag([mu, mu]);
B = [-M/h Wn Wf zeros(n,nc); 
     Wn' zeros(nc, nc + nd + nc);
     Wf' zeros(nd, nc + nd) E;
     zeros(nc, n) U  -E'  zeros(nc,nc)];
%%
% yG = @(t) t;
yG = @(t) 0;
phi_n = zeros(2,1);
phi_n(1) = yc - l/2*sin(theta) - r - yG(t); 
phi_n(2) = yc + l/2*sin(theta) - r - yG(t);
% bn = phi_n + h*[-1;-1] - Wn'*ql; % 2*1
bn = phi_n + h*[0;0] - Wn'*ql; % 4*1
bf = h*zeros(nd,1) - Wf'*ql; % 4*1
b = [f*h + M/h*ql; bn; bf; zeros(nc,1)];

N = length(b);
lb = zeros(N,1);
lb(1:n) = -Inf;
ub = Inf*ones(N,1);
%%
[z,mu] = pathlcp(B,b,lb,ub);