% libsrvf
% =======
%
% A shape analysis library using the square root velocity framework.
% 
% Copyright (C) 2012   FSU Statistical Shape Analysis and Modeling Group
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>
% --------------------------------------------------------------------


% Evaluates an SRVF at the given parameter values.
% 
% Inputs:
%  Q,T : the sample points and parameter values of the SRVF
%  t : the parameter values at which the SRVF will be evaluated.  
%      NOTE:  t must be non-decreasing.
%
% Outputs:
%  q : the SRVF function values for the parameter values in t
% --------------------------------------------------------------------------
function q = srvf_evaluate( Q, T, t )
  assert( size(T,1) == 1 );
  assert( size(Q,2) == size(T,2)-1 );
  assert( min(diff(T)) > 0);
  assert( min(diff(t)) >= 0 );
  epsval = 0.01*(T(end)-T(1));
  assert( t(1) > T(1)-epsval && t(end) < T(end)+epsval );

  q = zeros(size(Q,1),length(t));
  Qidx = 1;

  for qidx = 1:length(t)
    while ( Qidx < length(T)-1 && t(qidx) > T(Qidx+1) ) 
      Qidx = Qidx + 1; 
    end

    q(:,qidx) = Q(:,Qidx);
  end
end


%!test
%! Q=[0 1/2 1 -1];
%! T=[0 1/4 1/2 3/4 1];
%! t=[0 0.249 0.25 0.251 0.499 0.5 0.501 0.749 0.75 0.751 0.9 1];
%! qexp=[0 0 0 1/2 1/2 1/2 1 1 1 -1 -1 -1];
%! q=srvf_evaluate(Q,T,t);
%! assert(q,qexp,1e-3);
%!
%!#error
%! Q=[0 1/2 1 -1];
%! T=[0 1/4 1/2 3/4 1];
%! t=[-eps];
%! q=srvf_evaluate(Q,T,t);
%!
%!#error
%! Q=[0 1/2 1 -1];
%! T=[0 1/4 1/2 3/4 1];
%! t=[1+eps];
%! q=srvf_evaluate(Q,T,t);
%!
%!#error
%! Q=[0 1/2 1 -1];
%! T=[0 1/4 1/2 3/4 1];
%! t=[0 1/2 1/4 1];
%! q=srvf_evaluate(Q,T,t);
