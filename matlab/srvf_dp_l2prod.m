% libsrvf 
% =======
%
% A shape analysis library using the square root velocity framework.
% 
% Copyright (C) 2018   Martins Bruveris
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

% Helper function for srvf_dp_optimal_reparam. Calculates the L^2-inner
% product between Q1 and Q2 * G where G is the straight line between
% (i0, j0) and (i1, j1).
%
% Inputs
%  Q1, T1:     First SRVF
%  Q2, T2:     Second SRVF
%  i0, j0, i1, j1: Start and vertex where to compute the product
%
% Outputs
%  E:  Value of the L^2-inner product
function E = srvf_dp_l2prod(Q1, T1, Q2, T2, i0, j0, i1, j1)

% Rename input variables to something more familiar
x = T1;
u = Q1;
y = T2;
v = Q2;

% Crossing vertical grid lines
xCross = x(i0:i1);
xCross2 = [];

if j1 > j0 + 1
    % Crossing horizontal grid lines
    xCross2 = x(i0) + ...
        (y(j0+1:j1-1) - y(j0)) / (y(j1) - y(j0)) * (x(i1) - x(i0));
end
xCross = sort([xCross xCross2]);
yCross = y(j0) + (xCross - x(i0)) / (x(i1) - x(i0)) * (y(j1) - y(j0));

% Midpoints of intervals for finding correct indices
xMid = (xCross(1:end-1) + xCross(2:end)) / 2;
yMid = (yCross(1:end-1) + yCross(2:end)) / 2;

% Finding corresponding indices
xInd = interp1(x(i0:i1), i0:i1, xMid, 'previous');
yInd = interp1(y(j0:j1), j0:j1, yMid, 'previous');

% Construct integrand
g1p = diff(xCross);
g2p = diff(yCross);
ug = u(:, xInd);
vg = v(:, yInd);
wg = sum(ug .* vg, 1);

% Compute integral
E = sum( wg .* sqrt(g1p) .* sqrt(g2p) );

end