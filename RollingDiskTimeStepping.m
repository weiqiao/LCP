function [z, delta_q_u] = RollingDiskTimeStepping(qu, qa, qa_dot_d, h)
% A disk contacts the table at one point and rolls on the table.

qu_l = qu; % qu: [y, theta]
qa_l = qa; % qa: x
y = qu_l(1);
theta = qu_l(2);

% defining constants
r = 0.5; % m
mu = 0.6;
f = [-9.8;0];

nu = 2; % number of unactuated states
na = 1; % number of actuated states
n = nu + na;
nc = 1; % number of contacts
nf = 2;
nd = sum(nf); % no. of columns in Wf

Wf = [0 r;0 -r]';
Wn = [1 0]';

% Ja = zeros(nc+nd,na);
Ja = [0; 1; -1];

E = zeros(nc, nd);
first = 1;
for i = 1:1:nc
  E(i, first:first+nf(i)-1) = 1;
  first = first + nf(i);
end
E = E';
U = diag([mu]);
B = [zeros(nu,n) Wn Wf zeros(nu,nc);
     Wn' Ja(1:nc) zeros(nc, nc + nd + nc);
     Wf' Ja(nc+1:nc+nd) zeros(nd, nc + nd) E;
     zeros(nc, n) U  -E'  zeros(nc,nc)];
%%
% yG = @(t) t;
phi_n = zeros(nc,1);
phi_n(1) = y-r;
bn = phi_n; % 4*1
bf = zeros(nd,1); % 4*1
b = [f; bn; bf; zeros(nc,1)];

%% call solver
M = 50;
z = Mixed_LCP_As_MILP3(B,b,nu,na,h,qa_dot_d, M);
% delta_q_u = z(1:nu);

% QP to minimize delta_qu
epsilon = 1e-9;
delta_q_u = MinDeltaQu_QP(z, Wn, Wf, Ja, phi_n, nu, na, nc, nd, nf, epsilon);