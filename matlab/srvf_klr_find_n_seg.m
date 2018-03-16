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

% Helper function for srvf_klr_optimal_reparam. At a grid position (i0, j0)
% it finds all possible N-segments starting from this point.
%
% See the paper for the definition of an N-segment.
%
% Inputs
%  i0, j0:  Initial position for N-segment
%  w:       Weight matrix containing the inner products
%             w(i, j) = Q1(1)' * Q2(j)
%  w_eps:   Threshold below which we consider and entry of W to be
%           nonpositive.
%
% Outputs
%  endPoints: (2, numPts) matrix containing all possible endpoints.
% -------------------------------------------------------------------------
function endPoints = srvf_klr_find_n_seg(i0, j0, w, w_eps)

m = size(w, 1) + 1;
n = size(w, 2) + 1;

endPoints = zeros(m*n, 2);
numPoints = 0;

% Explore horizontally
i = i0;
j = j0;
continueHor = 1;
while continueHor
    if i == m % We have reached the end
        continueHor = 0;
    elseif (j == n || w(i, j) <= w_eps) && ...
           (j == 1 || w(i, j-1) <= w_eps)
        i = i + 1;
    else
        continueHor = 0;
    end
end
imax = i;

i = i0;
j = j0;
continueVer = 1;
while continueVer
    if j == n
        continueVer = 0;
    elseif i == 1 || w(i-1, j) <= w_eps
        j = j+1;
    else
        continueVer = 0;
    end
end
jmax = j;

for i = i0:imax
    j = j0;
    continueVer = 1;
    while continueVer && j <= jmax && j <= n
        if j == n || i == m
            endPoints(numPoints+1,1) = i;
            endPoints(numPoints+1,2) = j;
            numPoints = numPoints + 1;
            
            j = j + 1;
        elseif (i < m && w(i, j) <= w_eps) || (i==m && j < jmax)
            j = j + 1;
        else
            endPoints(numPoints+1,1) = i;
            endPoints(numPoints+1,2) = j;
            numPoints = numPoints + 1;
            
            jmax = j-1;
    
            continueVer = 0;
        end
    end
end

endPoints = endPoints(1:numPoints,:);