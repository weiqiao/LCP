function q_next = TimeSteppingLCP(q, t, h)
% q: current state, q = [xp, yp, zp]'
% t: current time
% h: time step
% q_next: state at t + h
% example in Trinkle's 05 paper on quasistatic LCP

ql = q;
xp = ql(1);
yp = ql(2);
zp = ql(3);

mu = 0.6;
nc = 3; % number of contacts
n = 3; % dim(q)

Wn = [-1 1 0; 0 0 0; 0 0 1];
U = diag([mu mu mu]);
Wf1 = [0 -1 0; 0 1 0]';
Wf2 = [0 1 0; 0 -1 0]';
Wf3 = [1 0 0;0 1 0;-1 0 0;0 -1 0]';
Wf = [Wf1, Wf2, Wf3];
E = zeros(nc, 8);
E(1,1:2) = 1;
E(2,3:4) = 1;
E(3,5:8) = 1;
E = E';

B = [zeros(n,n), Wn, Wf, zeros(n, nc);
     Wn', zeros(nc, nc + 8 + nc);
     Wf', zeros(8, nc + 8), E;
     zeros(nc, n), U, -E', zeros(nc, nc)];
f = [0 0 -9.8]';
xf = @(t) 0.5 + 0.4*sin(t); % x_fence
phi_n = [1-xp; xp - xf(t); zp];
bn = phi_n + [0; -cos(t); 0]*h - Wn'*ql;
bf = [0; 0; -1; 1; 0; 0; 0; 0]*h - Wf'*ql;
b = [f;bn;bf;zeros(nc,1)];
lb = zeros(17,1);
ub = Inf*ones(17,1);
lb(1:3) = -Inf;

%%
[z,mu] = pathlcp(B,b,lb,ub);
q_next = z(1:n);



