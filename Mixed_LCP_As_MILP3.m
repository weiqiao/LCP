function z = Mixed_LCP_As_MILP3(R, b, nu, na, h, qa_dot_d, M, z_start)
% objective:
% min ||q_a_desired(l+1)-q_a(l)|| 

n1 = nu + na;
A = R(1:nu,1:n1);
B = R(1:nu,n1+1:end);
C = R(nu+1:end,1:n1);
D = R(nu+1:end, n1+1:end);
b1 = b(1:nu);
b2 = b(nu+1:end);
n2 = length(b) - nu;

% define constraints
model.A = sparse([A, B, zeros(nu, n2);
                  C, D, zeros(n2, n2);
                  zeros(n2,n1), eye(n2), zeros(n2, n2);
                  C, D, -M*eye(n2);
                  zeros(n2,n1), eye(n2), M*eye(n2)]);
model.rhs = [-b1;-b2; zeros(n2,1); -b2; ones(n2,1)*M];

% define sense of constraints
Nc = n2*4 + nu; % total number of constraints
charArray = char(zeros(1,Nc));
charArray(1:nu) = '=';
charArray(nu+1:nu+2*n2) = '>';
charArray(nu+2*n2+1:nu+4*n2) = '<';
model.sense = charArray;

% set objective
N = n1 + n2*2; % number of decision variables
model.obj = zeros(N,1);
model.obj(nu+1:nu+na) = -2/h*qa_dot_d';
Q = zeros(N);
Q(nu+1:nu+na, nu+1:nu+na) = eye(na)/h^2;
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

% warm starting
model.start = z_start;

% set params
params.outputflag = 0;

result = gurobi(model, params);
% disp(result)
z = result.x;

