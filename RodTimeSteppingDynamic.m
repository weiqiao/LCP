function [z, q_l1, v_l1] = RodTimeSteppingDynamic(ql, vl, h)
% assuming friction force between table and rod occurs at both ends of the
% rod. 

% defining constants
r = 0.05; % m
l = 0.5; % m
m = 1;
J = 0.002;
mu = 0.6;
f = [0;-9.8;0];
M = diag([m, m, J]);

n = 3; % dimension of configuration q
nc = 2; % number of contacts

% fake q(l+1), denoting q(l+1) by q
q = ql + h*vl;
v = vl;

yc_l = ql(2);
theta_l = ql(3);

n1 = @(theta) [0 1 -l/2*cos(theta)]';
n2 = @(theta) [0 1 l/2*cos(theta)]';

n1_l = n1(theta_l);
n2_l = n2(theta_l);

yc = q(2);
theta = q(3);

Wf1 = @(theta) [-1 0 -l/2*sin(theta);1 0 l/2*sin(theta)]';
Wf2 = @(theta) [-1 0 l/2*sin(theta);1 0 -l/2*sin(theta)]';
Wf = [Wf1(theta), Wf2(theta)];

Wn = [n1(theta), n2(theta)];
   
nd = size(Wf,2); % no. of columns in Wf

E = zeros(nc, nd);
E(1,1:2) = 1;
E(2,3:4) = 1;
E = E';
U = diag([mu, mu]);
B = [h^2*Wn'*(M\Wn), h^2*Wn'*(M\Wf), zeros(nc, nc); 
     h*Wf'*(M\Wn), h*Wf'*(M\Wf), E;
     U  -E'  zeros(nc,nc)];
%%
alpha0 = [-(yc_l - l/2*sin(theta_l) - r) + n1_l'*ql; -(yc_l + l/2*sin(theta_l) - r) + n2_l'*ql];
u = [Wn'*(ql + h*vl + h^2*(M\f))-alpha0; Wf'*(vl + h*(M\f)); [0;0]];

N = length(u);
lb = zeros(N,1);
ub = Inf*ones(N,1);
%%
[z,mu] = pathlcp(B,u,lb,ub);
cn = z(1:nc);
beta = z(nc+1:3*nc);
lambda = z(3*nc+1:4*nc);
v_l1 = vl + M\(h*Wn*cn + h*Wf*beta + h*f);
q_l1 = h*v_l1 + ql;
