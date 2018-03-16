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


% Computes preimages of the values in xv under the non-decreasing 
% piecewise-linear function F,T.  The function must be a 1-D function, 
% and xv must be non-decreasing.
%
% If F is not a one-to-one function (i.e. F(j)==F(j+1) for some j), then
% there may be more than one preimage for a given value xv(i). In this case
% the xvi(i) is the midpoint of the first interval [T(j),T(j+1)], where
% F(j)==F(j+1). This should guarantee that if you do this
%
% xvi=plf_preimages(F,T,xv);
% Fxvi=plf_evaluate(F,T,xvi);
%
% then Fxvi will be equal (or at least close) to xv.
%
% Inputs
%  F,T:  the function.  F and T must be non-decreasing.
%  xv:   the sample points.  Must be non-decreasing.
%
% Outputs
%  xvi:  preimages of xv
% --------------------------------------------------------------------------
function xvi = plf_preimages( F, T, xv )
    assert( size(F,1) == 1 );
    assert( min(diff(F)) >= 0 );
    assert( min(diff(T)) >= 0 );
    
    epsilon = 1e-12;

    xvi = zeros(1,length(xv));
    j = 1;
    for i = 1:length(xv)
        % The value is at the beginning of T
        if j == 1 && xv(i) < F(j) - epsilon
            xvi(i) = -Inf;
            continue
        elseif j == 1 && xv(i) <= F(j) 
            % The value is borderline at the start
            xvi(i) = T(j);
            continue
        end

        % Finds smallest j such that F(j) <= xv(i) <= F(j+1)
        % Note: the condition xv(i) > F(j) is to catch the case when F is
        % constant on the interval [j, j+1]. In this case we want to stop
        % at the first such interval,
        while( j < length(F) && xv(i) > F(j) && xv(i) >= F(j+1) )
            j = j+1;
        end
        
        % The value lies above the image of F
        if j == length(F) && xv(i) > F(j) + epsilon
            xvi(i) = Inf;
            continue
        elseif j == length(F) % The value is borderline at the end
            xvi(i) = T(j);
            continue;
        end
        
        % F is constant on the interval [j, j+1]; return midpoint
        if F(j+1) - F(j) < epsilon
            xvi(i) = (T(j+1) + T(j)) / 2;
            continue
        end
        
        % The middle cases; F strictly increasing on [j, j+1]
        dF = F(j+1) - F(j);
        w1 = (F(j+1) - xv(i)) / dF;
        w2 = (xv(i) - F(j)) / dF;
        xvi(i) = w1 * T(j) + w2 * T(j+1);
    end
end


%!function v=_random_increasing_vector(v0,n)
%!  v=zeros(1,n);
%!  v(1)=v0;
%!  for i=2:n
%!    v(i)=v(i-1)+rand();
%!  end
%! end
%!
%!test
%! F=linspace(0,1,5);
%! T=linspace(0,1,5);
%! xv=linspace(0,1,100);
%! xvi=plf_preimages(F,T,xv);
%! assert(xvi,xv,1e-5);
%!
%!test
%! F=[0 0 1/2 1/2 1 1];
%! T=linspace(0,1,6);
%! tv=linspace(0,1,100);
%! xv=[0 0.2 0.499 0.5 0.501 0.999 1];
%! xvie=[0 0.28 0.3996 0.4 0.6004 0.7996 0.8];
%! xvi=plf_preimages(F,T,xv);
%! assert(xvi,xvie,1e-4);
%!
%!test
%! F=_random_increasing_vector(-100*rand(),50);
%! T=_random_increasing_vector(-5*rand(),50);
%! tv=linspace(T(1),T(end),100);
%! xv=interp1(T,F,tv);
%! xvi=plf_preimages(F,T,xv);
%! assert(xvi,tv,1e-4)
%!
%!test
%! F=_random_increasing_vector(0,1000);
%! F=F/F(end);
%! T=linspace(0,1,length(F));
%! tv=linspace(T(1),T(end),2000);
%! xv=interp1(T,F,tv);
%! xvi=plf_preimages(F,T,xv);
%! Fxvi=plf_evaluate(F,T,xvi);
%! assert(Fxvi,xv,1e-4)
