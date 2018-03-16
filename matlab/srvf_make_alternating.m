
function [Qn, Tn] = srvf_make_alternating(Q, T)
  Qn = [Q(1)];
  Tn = [T(1)];

  for i=2:size(Q, 2)
    if sign(Q(i)) ~= sign(Qn(end))
      Qn = [Qn Q(i)];
      Tn = [Tn T(i)];
    end
  end

  Tn = [Tn T(end)];
end
