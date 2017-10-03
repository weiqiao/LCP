function penetration1 = CalcIntegratedPenetration(Wn, Jna, delta_q_u, delta_q_a, phi_n, penetration)
% phi_n: normal distance function at t = l*h


% normal distance functions at t = (l+1)h
phi_n1 = phi_n + Wn'*delta_q_u + Jna*delta_q_a;

% penetration at t = (l+1)h
penetration1 = penetration;
for i = 1:length(phi_n1)
  if phi_n1(i) >= 0
    penetration1(i) = 0;
  else
    penetration1(i) = penetration(i) + phi_n1(i);
  end
end

