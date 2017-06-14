function z = Mixed_LCP_As_MILP2(R, b, n1, h)
% force balance equality constraint replaced by objective.

M = 50; %big M
A = R(1:n1,1:n1);
B = R(1:n1,n1+1:end);
C = R(n1+1:end,1:n1);
D = R(n1+1:end, n1+1:end);
b1 = b(1:n1);
b2 = b(n1+1:end);
n2 = length(b) - n1;

% define constraints
model.A = sparse([C, D, zeros(n2, n2+1);
                  zeros(n2,n1), eye(n2), zeros(n2, n2+1);
                  C, D, -M*eye(n2), zeros(n2,1);
                  zeros(n2,n1), eye(n2), M*eye(n2), zeros(n2,1);
                  A, B, zeros(n1,n2), -ones(n1,1)*M;
                  A, B, zeros(n1,n2), ones(n1,1)*M]);
model.rhs = [-b2; zeros(n2,1); -b2; ones(n2,1)*M;-b1;-b1];

% define sense of constraints
Nc = n2*4 + n1*2; % total number of constraints
charArray = char(zeros(1,Nc));
charArray(1:2*n2) = '>';
charArray(2*n2+1:4*n2) = '<';
charArray(4*n2+1:4*n2+n1) = '<';
charArray(4*n2+n1+1:end) = '>';
model.sense = charArray;

% set objective
N = n1 + n2*2 + 1; % number of decision variables
model.obj = [zeros(1,N-1), M];
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
z = result.x(1:n1+n2);

