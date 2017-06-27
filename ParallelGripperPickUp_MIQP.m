function [z, delta_q_u] = ParallelGripperPickUp_MIQP(qu, qa, qa_dot_d, h)
% The ball (point mass) is placed on a table. 

nu = 2; % number of unactuated states
na = 3; % number of actuated states
n = nu + na;
nc = 3; % number of contacts
nf = [2, 2, 2];
nd = sum(nf); % no. of columns in Wf

qu_l = qu;
qa_l = qa;
x = qu_l(1);
y = qu_l(2);
xl = qa_l(1);
xr = qa_l(2);
yg = qa_l(3);

% defining constants
r = 0.05; % m
mu = 0.5;
f = [0;-9.8];

Wf = [0  0 0  0 1 -1; 
      1 -1 1 -1 0  0]; 

Wn = [1 -1 0;
      0  0 1];
Ja = zeros(nc+nd,na);
Ja(1,1) = -1;
Ja(2,2) = 1;
Ja(4:5,3) = [-1;1];
Ja(6:7,3) = [-1;1];

E = zeros(nc, nd);
first = 1;
for i = 1:1:nc
  E(i, first:first+nf(i)-1) = 1;
  first = first + nf(i);
end
E = E';
U = diag([mu, mu, mu]);
B = [zeros(nu,n) Wn Wf zeros(nu,nc);
     Wn' Ja(1:nc, :) zeros(nc, nc + nd + nc);
     Wf' Ja(nc+1:nc+nd, :) zeros(nd, nc + nd) E;
     zeros(nc, n) U  -E'  zeros(nc,nc)];
   
phi_n = zeros(nc,1);
phi_n(1) = x - xl - r; 
phi_n(2) = xr - x - r;
phi_n(3) = y - r;
bn = phi_n; % 4*1
bf = zeros(nd,1); % 4*1
b = [f; bn; bf; zeros(nc,1)];

% call solver
M= 50;
z = Mixed_LCP_As_MILP3(B,b,nu,na,h,qa_dot_d, M);

% QP to minimize delta_qu
epsilon = 1e-6;
delta_q_u = MinDeltaQu_QP(z, Wn, Wf, Ja, phi_n, nu, na, nc, nd, nf, epsilon);
