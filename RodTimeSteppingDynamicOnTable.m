function [z, q_l1, v_l1] = RodTimeSteppingDynamicOnTable(ql, vl, h, t)
% assuming friction force between table and rod occurs at both ends of the
% rod. 

% defining constants
r = 0.05; % m
l = 0.5; % m
m = 1;
J = 1/12*m*l^2;
mu = 0.6;
f = [0;0;-9.8;0];
M = diag([m, m, m, J]);

n = 4; % dimension of configuration q
nc = 3; % number of contacts

% fake q(l+1), denoted by q
q = ql + h*vl;
v = vl;

yc_l = ql(2);
zc_l = ql(3);
theta_l = ql(4);
yc = q(2);
theta = q(3);

n1 = @(theta) [0 1 0 -l/2*cos(theta)]';
n2 = @(theta) [0 1 0 l/2*cos(theta)]';
n3 = [0 0 1 0]';
Wn = [n1(theta), n2(theta), n3];
Wn_l = [n1(theta_l), n2(theta_l), n3];

Wf1 = @(theta) [1 0 0 l/2*sin(theta);-1 0  0 -l/2*sin(theta)]';
Wf2 = @(theta) [1 0 0 -l/2*sin(theta);-1 0  0 l/2*sin(theta)]';
Wf3 = @(theta) [1 0 0 0;
                0 1 0 0;
               -1 0 0 0;
                0 -1 0 0;
                0 0 0 1;
                0 0 0 -1]';

Wf = [Wf1(theta), Wf2(theta), Wf3(theta)];
   
nd = size(Wf,2); % no. of columns in Wf

E = zeros(nc, nd);
E(1,1:2) = 1;
E(2,3:4) = 1;
E(3,5:10) = 1;
E = E';
U = diag([mu, mu, 0.1]);
B = [h^2*Wn'*(M\Wn), h^2*Wn'*(M\Wf), zeros(nc, nc); 
     h*Wf'*(M\Wn), h*Wf'*(M\Wf), E;
     U  -E'  zeros(nc,nc)];
%%
% yG = @(t) 0.25 + 0.1*t;
phi_n = zeros(3,1);
phi_n(1) = yc_l - l/2*sin(theta_l) - r -yG(t);
phi_n(2) = yc_l + l/2*sin(theta_l) - r -yG(t);
phi_n(3) = zc_l;
alpha0 = -(phi_n - Wn_l'*ql + h*[-1;-1;0]);
% alpha0 = [-(yc_l - l/2*sin(theta_l) - r) + n1_l'*ql; -(yc_l + l/2*sin(theta_l) - r) + n2_l'*ql];
u = [Wn'*(ql + h*vl + h^2*(M\f))-alpha0; Wf'*(vl + h*(M\f)); zeros(nc,1)];

N = length(u);
lb = zeros(N,1);
ub = Inf*ones(N,1);
%%
[z,mu] = pathlcp(B,u,lb,ub);
cn = z(1:nc);
beta = z(nc+1:nc + 2*2 + 6);
lambda = z(end-nc+1:end);
v_l1 = vl + M\(h*Wn*cn + h*Wf*beta + h*f);
q_l1 = h*v_l1 + ql;
