function z = Mixed_LCP_As_MILP2(R, b, nu, na, h, qa_dot_d)
% force balance equality constraint replaced by objective.

M = 50; %big M
n1 = nu + na;
A = R(1:nu,1:n1);
B = R(1:nu,n1+1:end);
C = R(nu+1:end,1:n1);
D = R(nu+1:end, n1+1:end);
b1 = b(1:nu);
b2 = b(nu+1:end);
n2 = length(b) - nu;

% define constraints
model.A = sparse([C, D, zeros(n2, n2+2);
                  zeros(n2,n1), eye(n2), zeros(n2, n2+2);
                  C, D, -M*eye(n2), zeros(n2,2);
                  zeros(n2,n1), eye(n2), M*eye(n2), zeros(n2,2);
                  A, B, zeros(nu,n2), -ones(nu,1)*M, zeros(nu,1);
                  A, B, zeros(nu,n2), ones(nu,1)*M, zeros(nu,1);
                  zeros(na, nu), eye(na), zeros(na, 2*n2+1), -M*ones(na,1);
                  zeros(na, nu), eye(na), zeros(na, 2*n2+1), M*ones(na,1)]);
model.rhs = [-b2; zeros(n2,1); -b2; ones(n2,1)*M;-b1;-b1;qa_dot_d*h;qa_dot_d*h];

% define sense of constraints
Nc = n2*4 + nu*2 + na*2; % total number of constraints
charArray = char(zeros(1,Nc));
charArray(1:2*n2) = '>';
charArray(2*n2+1:4*n2) = '<';
charArray(4*n2+1:4*n2+nu) = '<';
charArray(4*n2+nu+1:4*n2+2*nu) = '>';
charArray(4*n2+2*nu+1:4*n2+2*nu+na) = '<';
charArray(4*n2+2*nu+na+1:4*n2+2*nu+2*na) = '>';
model.sense = charArray;

% set objective
N = n1 + n2*2 + 2; % number of decision variables
model.obj = [zeros(1,N-2), M, M];
Q = zeros(N);
Q(1:n1, 1:n1) = eye(n1)/h^2;
model.Q = sparse(Q);

% set lower bound
lb = zeros(N,1);
lb(1:n1) = -Inf;
model.lb = lb;

% set variable types
charArray = char(zeros(1,N));
charArray(1:n1+n2) = 'C';
charArray(n1+n2+1:end) = 'B';
model.vtype = charArray;

% set params
params.outputflag = 0;

result = gurobi(model, params);
z = result.x;

