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


% Evaluates a PLF at the given parameter values.  The PLF is taken to be 
% left-continuous at any jump discontinuities.
%
% Inputs:
%  F,T - the sample points and parameter values for the PLF
%  t - the parameter values at which the PLF will be evaluated
%
% Outputs:
%  f - the function values for the parameter values in t
% --------------------------------------------------------------------------
function f = plf_evaluate( F, T, t )
  if isscalar(t)
    % This 'feature' has unintended consequences and has been removed.
    %t = linspace(0,1,t);
  end

  f = zeros(size(F,1),length(t));
  for i=1:size(F,1)
    [Tunique, ind] = unique(T);
    f(i,:) = interp1(Tunique,F(i,ind),t,'linear','extrap');
  end

  % left-continuous at jump discontinuities
  jumps = find( diff(T) == 0 );
  for j=1:length(jumps)
    jump_tidx = find( t == T(jumps(j)) );
    for i=1:size(F,1)
      f(i,jump_tidx) = F(i,jumps(j));
    end
  end
end


%!test
%! F=[0 1 2 3 4 5 6;11 10 9 8 7 6 5];
%! T=[0 0 1/4 1/4 1/2 1 1];
%! tv=[0 1/8 1/4 3/8 1/2 3/4 1];
%! X=[0.0000   1.5000   2.0000   3.5000   4.0000   4.5000   5.0000;\
%!    11.0000    9.5000    9.0000    7.5000    7.0000    6.5000    6.0000];
%! A=plf_evaluate(F,T,tv);
%! assert(A,X);
