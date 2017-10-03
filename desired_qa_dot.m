function qa_dot = desired_qa_dot(i, is_MIQP)
if(is_MIQP)
  if i > 39
    qa_dot = [-0.1;0.1;0];
  % elseif i == 80
  %   qa_dot = [0;0;0.1];
  else
    qa_dot = [0.1;-0.1;0.1];
  end
else
  if i > 20
    qa_dot = [0;0;0.1];
  else
    qa_dot = [0.1;-0.1;0.1];
  end
end
