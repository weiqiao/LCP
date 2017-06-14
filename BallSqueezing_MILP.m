function z = BallSqueezing_MILP(qu, qa, qa_dot_d, h)
% assuming friction force between table and rod exists at both ends of the
% rod. 
% qa = yG(t) \in R^1
% qa_dot = d{yG(t)}/dt;

qu_l = qu;
qa_l = qa;
x = qu_l(1);

% defining constants
r = 0.05; % m
mu = 0.6;
f = [0;-9.8];

nu = 2; % number of unactuated states
na = 1; % number of actuated states
n = nu + na;
nc = 2; % number of contacts

Wf = [0 0 0 0
      1 -1 1 -1];
nd = size(Wf,2); % no. of columns in Wf

Wn = [1 -1
      0 0];

Ja = zeros(nc+nd,1);
Ja(2) = 1;

E = zeros(nc, nd);
E(1,1:2) = 1;
E(2,3:4) = 1;
E = E';
U = diag([mu, mu]);
B = [zeros(nu,n) Wn Wf zeros(nu,nc);
     Wn' Ja(1:nc) zeros(nc, nc + nd + nc);
     Wf' Ja(nc+1:nc+nd) zeros(nd, nc + nd) E;
     zeros(nc, n) U  -E'  zeros(nc,nc)];
%%
% yG = @(t) t;
phi_n = zeros(2,1);
phi_n(1) = x - r; 
phi_n(2) = qa_l - (x + r);
bn = phi_n; % 4*1
bf = zeros(nd,1); % 4*1
b = [f; bn; bf; zeros(nc,1)];

%% call solver
z = Mixed_LCP_As_MILP2(B,b,nu,na,h,qa_dot_d);
z(1:n) = z(1:n) + [qu_l;qa_l];