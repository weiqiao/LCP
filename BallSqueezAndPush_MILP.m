function z = BallSqueezAndPush_MILP(qu, qa, qa_dot_d, h)
% assuming friction force between table and rod exists at both ends of the
% rod. 
% qa = yG(t) \in R^1
% qa_dot = d{yG(t)}/dt;

qu_l = qu;
qa_l = qa;
x = qu_l(1);
y = qu_l(2);
xp = qa_l(1);
yp = qa_l(2);

% defining constants
r = 0.05; % m
mu = 0;
f = [0;-9.8];

nu = 2; % number of unactuated states
na = 2; % number of actuated states
n = nu + na;
nc = 3; % number of contacts

Wf = [0 0 0 0 1 -1 
      1 -1 1 -1 0 0];  
nd = size(Wf,2); % no. of columns in Wf

Wn = [1 -1 0
      0 0 1];

Ja = zeros(nc+nd,na);
Ja(2,1) = 1;
Ja(3,2) = -1;

E = zeros(nc, nd);
E(1,1:2) = 1;
E(2,3:4) = 1;
E(3,5:6) = 1;
E = E';
U = diag([mu, mu, mu]);
B = [zeros(nu,n) Wn Wf zeros(nu,nc);
     Wn' Ja(1:nc, :) zeros(nc, nc + nd + nc);
     Wf' Ja(nc+1:nc+nd, :) zeros(nd, nc + nd) E;
     zeros(nc, n) U  -E'  zeros(nc,nc)];
%%
% yG = @(t) t;
phi_n = zeros(3,1);
phi_n(1) = x - r; 
phi_n(2) = xp - (x + r);
phi_n(3) = y - yp - r;
bn = phi_n; % 4*1
bf = zeros(nd,1); % 4*1
b = [f; bn; bf; zeros(nc,1)];

%% call solver
M1 = 50;
M2 = 0.1;
z = Mixed_LCP_As_MILP2(B,b,nu,na,h,qa_dot_d, M1, M2);
z(1:n) = z(1:n) + [qu_l;qa_l];