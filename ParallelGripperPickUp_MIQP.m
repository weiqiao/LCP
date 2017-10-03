function [z, delta_q_u, penetration1] = ParallelGripperPickUp_MIQP(qu, qa, qa_dot_d, h, penetration, is_MIQP)
% The ball (point mass) is placed on a table. 

nu = 2; % number of unactuated states
na = 3; % number of actuated states
n = nu + na;
nc = 3; % number of contacts
nf = [2, 2, 2];
nd = sum(nf); % no. of columns in Wf
n2 = nd + nc*2;

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
      1 -1 1 -1 0  0;]; 

Wn = [1 -1 0;
      0  0 1;];
Ja = zeros(nc+nd,na);
Ja(1,1) = -1;
Ja(2,2) = 1;
Ja(nc+1:nc+2, 3) = [-1;1];
Ja(nc+3:nc+4, 3) = [-1;1];

E = zeros(nc, nd);
first = 1;
for i = 1:1:nc
  E(i, first:first+nf(i)-1) = 1;
  first = first + nf(i);
end
E = E';
U = diag([mu, mu, mu]);
delta_q_a_d = qa_dot_d*h;

% call solver
if(is_MIQP)
  B = [zeros(nu,n) Wn Wf zeros(nu,nc);
     Wn' Ja(1:nc, :) zeros(nc, nc + nd + nc);
     Wf' Ja(nc+1:nc+nd, :) zeros(nd, nc + nd) E;
     zeros(nc, n) U  -E'  zeros(nc,nc)];
   
  phi_n = zeros(nc,1);
  phi_n(1) = x - xl - r; 
  phi_n(2) = xr - x - r;
  phi_n(3) = y - r;
  b = [f; phi_n; zeros(nd,1); zeros(nc,1)];
  
  M= 50;
  z = Mixed_LCP_As_MILP_w_normal_force_lower_bound(B,b,nu,na,nc,h,qa_dot_d, M, penetration);
  delta_q_u = z(1:nu);
  
else
  % solve by LCP, first nu entries of z is qu rather than delta_qu
  B = [zeros(nu,nu) Wn Wf zeros(nu,nc);
   Wn' zeros(nc, nc + nd + nc);
   Wf' zeros(nd, nc + nd) E;
   zeros(nc, nu) U  -E'  zeros(nc,nc)];
   
  phi_n = zeros(nc,1);
  phi_n(1) = x - xl - r; 
  phi_n(2) = xr - x - r;
  phi_n(3) = y - r;
  bn = phi_n + Ja(1:nc, :)*qa_dot_d*h;
  bf = Ja(nc+1:nc+nd, :)*qa_dot_d*h;
  b = [f; bn; bf; zeros(nc,1)];
  
  N = length(b);
  lb = zeros(N,1);
  lb(1:nu) = -Inf;
  ub = Inf*ones(N,1);
  [z,mu] = pathlcp(B,b,lb,ub);
  delta_q_u = z(1:nu);
  % adding zeros to match the length of z in the MIQP version
  z = [z(1:nu);delta_q_a_d;z(nu+1:end);zeros(n2,1)];
end
penetration1 = CalcPenetration(Wn, Ja(1:nc, :), delta_q_u, delta_q_a_d, phi_n, penetration);
% QP to minimize delta_qu
%epsilon = 1e-6;
%delta_q_u = MinDeltaQu_QP(z, Wn, Wf, Ja, phi_n, nu, na, nc, nd, nf, epsilon);
