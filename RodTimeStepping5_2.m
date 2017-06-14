function z = RodTimeStepping5_2(q, t, h)
% assuming generalized friction between table and rod exists at the
% centroid of the rod.
% "quasi-dynamic"
%

ql = q;
xc = ql(1);
yc = ql(2);
zc = ql(3);
theta = ql(4);

% defining constants
r = 0.05; % m
l = 0.5; % m
mu = 0.6;
f = [0;0;-9.8;0];

m = 1;
M = diag([m, m, m, m*l^2/12]);

n = 4; % dimension of configuration q
nc = 3; % number of contacts

Wf1 = @(theta) [1 0 0 l/2*sin(theta);-1 0  0 -l/2*sin(theta)]';
Wf2 = @(theta) [1 0 0 -l/2*sin(theta);-1 0  0 l/2*sin(theta)]';
Wf3 = @(theta) [1 0 0 0;
                0 1 0 0;
               -1 0 0 0;
                0 -1 0 0;
                0 0 0 1;
                0 0 0 -1]';

Wf = [Wf1(theta), Wf2(theta), Wf3(theta)];

n1 = @(theta) [0 1 0 -l/2*cos(theta)]';
n2 = @(theta) [0 1 0 l/2*cos(theta)]';
n3 = [0 0 1 0]';
Wn = [n1(theta), n2(theta), n3];
   
nd = size(Wf,2); % no. of columns in Wf

E = zeros(nc, nd);
E(1,1:2) = 1;
E(2,3:4) = 1;
E(3,5:10) = 1;
E = E';
U = diag([mu, mu, 0.1]);
B = [-M/h Wn Wf zeros(n,nc); 
     Wn' zeros(nc, nc + nd + nc);
     Wf' zeros(nd, nc + nd) E;
     zeros(nc, n) U  -E'  zeros(nc,nc)];
%%
% yG = @(t) 0.25 + 0.1*t;
phi_n = zeros(3,1);
phi_n(1) = yc - l/2*sin(theta) - r -yG(t); 
phi_n(2) = yc + l/2*sin(theta) - r -yG(t);
phi_n(3) = zc; 
bn = phi_n + h*[-1;-1;0] - Wn'*ql; % 4*1
bf = h*zeros(nd,1) - Wf'*ql; % 4*1
b = [f*h + M/h*ql; bn; bf; zeros(nc,1)];

N = length(b);
lb = zeros(N,1);
lb(1:n) = -Inf;
ub = Inf*ones(N,1);
%%
[z,mu] = pathlcp(B,b,lb,ub);