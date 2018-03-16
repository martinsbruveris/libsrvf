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


% Plots the matching between two PLFs.
%
% Input
%  F1,T1 : the first PLF
%  F2,T2 : the second PLF
%  fmt1 : a format string accepted by plot()
%  fmt2 : a format string accepted by plot()
% --------------------------------------------------------------------------
function plot_registration( F1, T1, F2, T2, fmt1, fmt2 )
  if ( nargin < 6 )
    fmt2 = 'r';
  end
  if ( nargin < 5 )
    fmt1 = 'b';
  end
  Tr = unique( [T1 T2] );

  dim = size(F1,1);
  nsamps = length(Tr);

  F1r = plf_evaluate(F1,T1,Tr);
  F2r = plf_evaluate(F2,T2,Tr);

  if ( dim == 1 )
    plot( Tr, F1r, fmt1, Tr, F2r, fmt2 );
  elseif ( dim == 2 )
    figure();
    hold on;

    plot( F1(1,:), F1(2,:), fmt1 );
    plot( F2(1,:), F2(2,:), fmt2 );

    for i=1:3:nsamps
      plot( [F1r(1,i) F2r(1,i)], [F1r(2,i) F2r(2,i)], 'k' );
    end
  elseif ( dim == 3 )
    figure();
    hold on;

    plot3( F1(1,:), F1(2,:), F1(3,:), fmt1 );
    plot3( F2(1,:), F2(2,:), F2(3,:), fmt2 );

    for i=1:3:nsamps
      plot3( [F1r(1,i) F2r(1,i)], ...
             [F1r(2,i) F2r(2,i)], ...
             [F1r(3,i) F2r(3,i)], 'k' );
    end
  else
    error 'Unsupported dimension';
  end
end


%!demo
%! load ../tests/data/horse-1.csv
%! load ../tests/data/horse-2.csv
%! X1 = horse_1;
%! X2 = horse_2;
%! T1 = linspace(0,1,length(X1));
%! T2 = linspace(0,1,length(X2));
%! plot_registration(X1,T1,X2,T2);

