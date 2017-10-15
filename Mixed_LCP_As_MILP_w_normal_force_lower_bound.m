function z = Mixed_LCP_As_MILP_w_normal_force_lower_bound(R, b, nu, na, nc, h, qa_dot_d, M, penetration)
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
phi_n = b2(1:nc);

% integrator with anti-windup
K = 1000;
% lambda_n_lb = K*abs(penetration);
% for i = 1:length(lambda_n_lb)
%   if lambda_n_lb(i) > M
%     lambda_n_lb(i) = M;
%   end
% end
lambda_n_lb = zeros(length(penetration), 1);
for i = 1:length(lambda_n_lb)
  if penetration(i) < -1e-6
    lambda_n_lb(i) = 10;
  end
end
% big M for lambda_n and phi_n
M_phi = 1e5;
M_f = M_phi/K; 

% define constraints
model.A = sparse([A, B, zeros(nu, n2);
                  C, D, zeros(n2, n2);
                  zeros(n2,n1), eye(n2), zeros(n2, n2);
                  C, D, -M*eye(n2);
                  zeros(n2,n1), eye(n2), M*eye(n2);
                  zeros(nc,n1), eye(nc), zeros(nc, 2*n2-nc)]);  % lambda_lower_bound
model.rhs = [-b1;-b2; zeros(n2,1); -b2; ones(n2,1)*M; lambda_n_lb(:)];


% define sense of constraints
Nc = n2*4 + nu + nc; % total number of constraints
charArray = char(zeros(1,Nc));
charArray(1:nu) = '=';
charArray(nu+1:nu+2*n2) = '>';
charArray(nu+2*n2+1:nu+4*n2) = '<';
charArray(nu+4*n2+1:nu+4*n2+nc) = '>';
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
lb(1:n1) = -0.1;
ub = ones(N,1) * inf;
ub(1:n1) = 0.1;
model.lb = lb;
model.ub = ub;

% set variable types
charArray = char(zeros(1,N));
charArray(1:n1+n2) = 'C';
charArray(n1+n2+1:end) = 'B';
model.vtype = charArray;

% set params
params.outputflag = 0;

result = gurobi(model, params);
% disp(result)
z = result.x;

