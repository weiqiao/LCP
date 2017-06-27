function delta_q_u_opt = MinDeltaQu_QP(z, Wn, Wf, Ja, phi_n, nu, na, nc, nd, nf, epsilon)

n1 = nu + na;
delta_q_u = z(1:nu);
delta_q_a = z(nu+1:nu+na);
% integer variables corresponding to lambda_n. lambda_n(i) != 0 iff z_lambda_n(i) == 0
lambda_n = z(n1+1:n1+nc);
% integer variables corresponding to gamma. gamma(i)!=0 iff z_gamma(i)==0
gamma_n = z(n1+nc+nd+1:n1+nc+nd+nc);
delta_phi_f = Wf'*delta_q_u + Ja(nc+1:nc+nd, :)*delta_q_a;

Wf_half = zeros(nu,nd/2);
i_half_start = 1;
signs_n = char(zeros(1,nc));
signs_f = char(zeros(1,nd));
signs_e = char(zeros(1,nc)); % equality constraint when sliding occurs in more than one direction
n_f = []; % indices of rows of Wf' that enter constraints.
Eq = zeros(nu, nc);
count_equal = 0;

for i = 1:nc
  if abs(lambda_n(i)) > epsilon % if lambda_n(i) > 0
    signs_n(i) = '=';
    i_start = sum(nf(1:i-1)) + 1;
    i_end = i_start + nf(i)/2-1;
    n_f = [n_f, i_start:1:i_end];
    Wf_half(:,i_half_start:(i_half_start+nf(i)/2-1)) = Wf(:,i_start:i_end);
    
    % constraints 
    if abs(gamma_n(i)) < epsilon % gamma(i) == 0
      signs_f(i_half_start:(i_half_start+nf(i)/2-1)) = '=';
    else
      for j = 1:nf(i)/2
        if abs(delta_phi_f(i_start+j-1)) < epsilon % == 0 
          signs_f(i_half_start+j-1) = '=';
        elseif delta_phi_f(i_start+j-1) < 0;
          signs_f(i_half_start+j-1) = '<';
        else
          signs_f(i_half_start+j-1) = '>';
        end
      end
      
      % When sliding in two directions, the sliding speeds in both
      % directions are equal.
      % This implementation only works when a plane is positively spanned by FOUR vectors!
      if nf(i)/2>1 && prod(delta_phi_f(i_start:i_end))>epsilon
        count_equal = count_equal + 1;
        Eq(:,count_equal) = Wf(:,i_start) - Wf(:,i_end);
        signs_e(count_equal) = '=';
      end
    end    
    i_half_start = i_half_start + nf(i)/2;
  else
    signs_n(i) = '>';
  end
end

% trim constraints
signs_f = signs_f(1:length(n_f));
Wf_half = Wf_half(:, 1:length(n_f));
signs_e = signs_f(1:count_equal);
Eq = Eq(:, 1:count_equal);

% define constraints
model.A = sparse([Wn';Wf_half';Eq']);
model.rhs = -Ja([1:nc, n_f+nc],:)*delta_q_a - [phi_n;zeros(length(signs_f),1)];
model.rhs = [model.rhs; zeros(count_equal, 1)];
model.sense = [signs_n,signs_f,signs_e];

% set objective
model.Q = sparse(eye(nu));
model.obj = zeros(nu,1);

% set lower bound
lb = zeros(nu,1);
lb(:) = -Inf;
model.lb = lb;

% set variable types
charArray = char(zeros(1,nu));
charArray(:) = 'C';
model.vtype = charArray;

% set params
params.outputflag = 0;

result = gurobi(model, params);
% disp(result)
delta_q_u_opt = result.x;