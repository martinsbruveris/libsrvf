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

% Helper function for srvf_klr_optimal_reparam. Constructs the N-segment
% starting from (i0, j0) to (i1, j1).
%
% We first proceed horizontally and then vertically.
%
% Inputs
%  i0, j0:     Initial position for N-segment
%  x, y:       Corresponds to T1, T2 in optimal_reparam
%  i1, j1:     Endpoint of N-segment.
%
% Outputs
%  gamma:    Points where the N-segment intersects the xy-grid.
%  gammaInd: If gammaInd(k) = [i; j], then 
%              x(i) <= gamma(1, k) <= x(i+1)
%              y(j) <= gamma(2, k) <= y(j+1)
% -------------------------------------------------------------------------
function [gamma, gammaInd] = srvf_klr_follow_n_seg(i0, j0, x, y, i1, j1)

m = size(x, 1);
n = size(y, 1);

gamma = zeros(2, m+n-1);
gamma(:, 1) = [x(i0); y(j0)];
gammaInd = zeros(2, m+n-1);
gammaInd(:, 1) = [i0; j0];
numGamma = 1;

for i=i0+1:i1
    gamma(:, numGamma+1) = [x(i); y(j0)];
    gammaInd(:, numGamma+1) = [i; j0];
    numGamma = numGamma + 1;
end

for j=j0+1:j1
    gamma(:, numGamma+1) = [x(i1); y(j)];
    gammaInd(:, numGamma+1) = [i1; j];
    numGamma = numGamma + 1;
end

gamma = gamma(:, 1:numGamma);
gammaInd = gammaInd(:, 1:numGamma);

end