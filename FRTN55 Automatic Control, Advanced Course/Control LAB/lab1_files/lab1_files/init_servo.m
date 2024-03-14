% Solution to problem with uninitialized LTI systems

if ~exist('F')
  F = tf(1);
end
if ~exist('C')
  C = tf(1);
end
if ~exist('P_sim')
  if exist('P')
    P_sim = P;
  else
    P_sim = tf(1);
  end
end
if ~exist('km')
  km = 2.95;
end
