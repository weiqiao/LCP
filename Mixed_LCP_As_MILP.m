function z = Mixed_LCP_As_MILP(R,b,n1)

M = 50; %big M
A = R(1:n1,1:n1);
B = R(1:n1,n1+1:end);
C = R(n1+1:end,1:n1);
D = R(n1+1:end, n1+1:end);
b1 = b(1:n1);
b2 = b(n1+1:end);
n2 = length(b) - n1;
model.A = sparse([A, B, zeros(n1,n2);
                  C, D, zeros(n2);
                  zeros(n2,n1), eye(n2), zeros(n2);
                  C, D, -M*eye(n2);
                  zeros(n2,n1), eye(n2), M*eye(n2)]);
model.rhs = [-b1; -b2; zeros(n2,1); -b2; ones(n2,1)*M];

% set objective
N = n1 + n2*2;
model.obj = zeros(N,1);
Q = zeros(N);
Q(1:n1, 1:n1) = eye(n1);
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

% define sense of constraints
Nc = n1 + n2 * 4; % total number of constraints
charArray = char(zeros(1,Nc));
charArray(1:n1) = '=';
charArray(n1+1:n1+2*n2) = '>';
charArray(n1+2*n2+1:end) = '<';
model.sense = charArray;

% set params
params.outputflag = 0;

result = gurobi(model, params);
z = result.x(1:n1+n2);

