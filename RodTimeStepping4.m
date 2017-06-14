function z = RodTimeStepping4(q, t, h)
% assuming no table and no gravity, i.e. table moves up to contact a rod in
% space. Point 2 is considered fixed in space. 
ql = q;
xc = ql(1);
yc = ql(2);
theta = ql(3);

% defining constants
r = 0.05; % m
l = 0.5; % m
mu = 0.6;
f = [0;0;0];

n = 3; % dimension of configuration q
nc = 5; % number of contacts

Wf1 = @(theta) [-1 0 -l/2*sin(theta);1 0 l/2*sin(theta)]';
% Wf2 = @(theta) [-1 0 l/2*sin(theta);1 0 -l/2*sin(theta)]';
Wf = [Wf1(theta)];

n1 = @(theta) [0 1 -l/2*cos(theta)]';
n2 = @(theta) [1 0 -l/2*sin(theta)]';
n4 = @(theta) [0 1 l/2*cos(theta)]';
Wn = [n1(theta), n2(theta), -n2(theta), n4(theta), -n4(theta)];
   
nd = size(Wf,2); % no. of columns in Wf, = 2

E = [1;1];
U = [mu 0 0 0 0];
B = [zeros(n,n) Wn Wf zeros(n,1); 
     Wn' zeros(nc, nc + nd + 1);
     Wf' zeros(nd, nc + nd) E;
     zeros(1, n) U  -E'  zeros(1,1)];
%%
yG = @(t) t;
phi_n = zeros(5,1);
phi_n(1) = yc - l/2*sin(theta) - r - yG(t); 
phi_n(2) = xc + l/2*cos(theta);
phi_n(3) = -phi_n(2);
phi_n(4) = yc + l/2*sin(theta);
phi_n(5) = -phi_n(4);
bn = phi_n + h*[-1;0;0;0;0] - Wn'*ql; % 2*1
bf = h*zeros(nd,1) - Wf'*ql; % 4*1
b = [f; bn; bf; zeros(1,1)];

N = length(b);
lb = zeros(N,1);
lb(1:n) = -Inf;
ub = Inf*ones(N,1);
%%
[z,mu] = pathlcp(B,b,lb,ub);