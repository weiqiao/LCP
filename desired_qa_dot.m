function qa_dot = desired_qa_dot(i)

if i > 80
  qa_dot = [-0.1;0.1;0.1];
% elseif i == 80
%   qa_dot = [0;0;0.1];
else
  qa_dot = [0.1;-0.1;0.1];
end
