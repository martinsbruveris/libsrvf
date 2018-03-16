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
% it follows the P-segment with slope al assuming that the segment ends at
% the point (iend, jend).
%
% The algorithm uses the change of slopes relation from the paper.
%
% Inputs
%  i0, j0:     Initial position for P-segment
%  x, y:       Corresponds to T1, T2 in optimal_reparam
%  w:          Weight matrix containing the inner products
%                w(i, j) = Q1(1)' * Q2(j)
%  al:         Initial angle for the P-segment
%  iend, jend: Endpoint of the P-segment. If not known in advance, set
%                iend=0, jend=0
%              But then one needs to adjust gridpoint detection.
%  w_eps:      Threshold below which we consider and entry of W to be
%              nonpositive.
%
% Outputs
%  gamma:    Points where the P-segment intersects the xy-grid.
%  gammaInd: If gammaInd(k) = [i; j], then 
%              x(i) <= gamma(1, k) <= x(i+1)
%              y(j) <= gamma(2, k) <= y(j+1)
% -------------------------------------------------------------------------

function [gamma, gammaInd] = ...
    srvf_klr_follow_p_seg(i0, j0, x, y, w, al, iend, jend, w_eps)

% u(i) = q1 restricted to [x(i), x(i+1)]
% v(j) = q2 restricted to [y(j), y(j+1)]
%x = q1(i0:end,1);
%u = q1(i0:end-1, 2:end);
%y = q2(j0:end,1);
%v = q2(j0:end-1, 2:end);

if i0 == 8 && j0 == 5 && iend==19 && jend==19
    disp('Break');
end

m = size(x, 2);
n = size(y, 2);

% Follow the guess for a P-segment with slope h
alNow = al;
wNow = w(i0, j0);

i = i0;
j = j0;
xNow = x(i0);
yNow = y(j0);

gamma = -1 * ones(2, m+n-1);
gamma(:,1) = [x(i0); y(j0)];
gammaInd = zeros(2, m+n-1);
numGamma = 1;

crossType = zeros(1, m+n-1);
alHist = zeros(1, m+n-1);

alPrev = al; % Last non-zero slope
wPrev = wNow;

while i < m && j < n
    % I have crossed the grid at the point (xNow, yNow) and I know that the
    % last non-zero slope is hPrev. The point (xNow, yNow) lies in the
    % rectangle G(i,j) = (x(i),x(i+1)) x (y(i),y(i+1)).
    
    % We add the current rectangle to gamma
    gammaInd(1, numGamma) = i;
    gammaInd(2, numGamma) = j;
    
    alHist(numGamma) = alNow;
    
    % We know a priori that the next step will get us to where we want
    if iend ~= 0 && jend ~= 0 && ...
            i+1 == iend && j+1 == jend
        gamma(:, numGamma+1) = [x(iend); y(jend)];
        gammaInd(:, numGamma+1) = [iend; jend];
        numGamma = numGamma + 1;
        % disp('Segment found.');
        break
    end
    
    wNow = w(i, j);
    if wNow > w_eps
        % Find slope that would get us to the corner [i+1,j+1]
        alCorner = atan2(y(j+1) - yNow, x(i+1) - xNow);
        
        % Determine how we cross into next cell
        if abs(alNow - alCorner) < eps
            gamma(1, numGamma+1) = x(i+1);
            gamma(2, numGamma+1) = y(j+1);
            gammaInd(1, numGamma+1) = i+1;
            gammaInd(2, numGamma+1) = j+1;
            numGamma = numGamma + 1;
            % disp('Segment found!');
            break
        elseif alNow < alCorner
            crossToRight = 1;
            xNext = x(i+1);
            yNext = yNow +(xNext - xNow) * tan(alNow);
        else
            crossToRight = 0;
            yNext = y(j+1);
            xNext = xNow + (yNext - yNow) * tan(pi/2 - alNow);
        end
    elseif crossFromRight
        % Cross horizontally
        xNext = x(i+1);
        yNext = yNow;            
        crossToRight = 1;
        
        crossType(numGamma) = 1;
    else % crossFromRight = 0
        % Cross vertically
        xNext = xNow;
        yNext = y(j+1);
        crossToRight = 0;
        
        crossType(numGamma) = -1;
    end
    
    if crossToRight
        i = i + 1;    
    else
        j = j + 1;
    end
    
    if i == m || j == n
        disp('Segment not found.');
        break
    end
    
    % Add point to gamma
    gamma(:, numGamma+1) = [xNext; yNext];
    gammaInd(:, numGamma+1) = [i; j];
    numGamma = numGamma + 1;
    
    % Update last nontrivial angle
    if wNow > w_eps
        alPrev = alNow;
        wPrev = wNow;
    end
    
    % Calculate next angle; note i and j have been updated
    wNext = w(i, j);
    if wNext <= w_eps
        alNow = Inf;
    elseif crossToRight
        alNow = atan2(wNext^2 * sin(alPrev), wPrev^2 * cos(alPrev));
    else % crossToRight=0
        alNow = atan2(wPrev^2 * sin(alPrev), wNext^2 * cos(alPrev));
    end
        
    % Update current position
    xNow = xNext;
    yNow = yNext;
    crossFromRight = crossToRight;
    
end

gamma = gamma(:, 1:numGamma);
gammaInd = gammaInd(:, 1:numGamma);

assert(gammaInd(1,end) == iend);
assert(gammaInd(2,end) == jend);

end