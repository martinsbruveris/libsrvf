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
% it finds the next P-segment from the point with slope larger than almin.
%
% The algorithm adopts the rescaling of the grid approach detailed in the
% paper.
%
% Inputs
%  i0, j0:  Initial position for P-segment
%  x, y:    Corresponds to T1, T2 in optimal_reparam
%  w:       Weight matrix containing the inner products
%             w(i, j) = Q1(1)' * Q2(j)
%  almin:   Angle determining the lower bound for the slope of the
%           P-segment
%  al_eps:  Smallest step we want to make with alpha
%  w_eps:   Threshold below which we consider and entry of W to be
%           nonpositive.
%
% Outputs
%  al:       Slope of the P-segment.
%  gamma:    Points where the P-segment intersects the xy-grid.
%  gammaInd: If gammaInd(k) = [i; j], then 
%              x(i) <= gamma(1, k) <= x(i+1)
%              y(j) <= gamma(2, k) <= y(j+1)
% -------------------------------------------------------------------------
function [al, gamma, gammaInd] = ...
    srvf_klr_find_next_p_seg(i0, j0, x, y, w, almin, al_eps, w_eps)

m = length(x);
n = length(y);

% We start with angle al0 and see where following the rules for a P-segment
% would lead us. At the same time we reparametrize the domain so that the
% P-segment becomes a straight line.
al0 = almin;
w0 = w(i0, j0); % Assumption: w0 > 0; Otherwise we have an N-segment
if w0 <= w_eps
    al = -Inf;
    gamma = [];
    return
end

% Our starting point
i = i0;
j = j0;
xNow = x(i0);
yNow = y(j0);

% The gridpoints of the rescaled grid
xTilde = -1 * ones(1, m);
yTilde = -1 * ones(1, n);
xTilde(i0) = x(i0); % The first rectangle remains the same
xTilde(i0+1) = x(i0+1);
yTilde(j0) = y(j0);
yTilde(j0+1) = y(j0+1);

% Here we save the scaling factors going between the variables
% Note: W_new(i,j) = W(i,j) * sqrt(xScale(i)) * sqrt(yScale(j))
xScale = ones(1, m-1);
yScale = ones(1, n-1);

% Save gridpoints above the line for search afterwards
iAbove = zeros(1, m);
jAbove = zeros(1, m);
indAbove = zeros(1, m);
numAbove = 0;

crossPtsTilde = zeros(2, n+m-1);
crossType = zeros(1, n+m-1); % 1 .. to right, 0 .. to top
numPts = 0;

% How did we cross into this rectangle. It does not matter for the first
% one and so we choose that we came from the side.
crossFromRight = 1;

% In case we are lucky and we find a P-segment immediatly
segmentFound = 0;
iend = 0;
jend = 0;

while i < m && j < n
    % I have crossed the grid at the point (xNow, yNow). The point 
    % (xNow, yNow) lies in the rectangle 
    % G(i,j) = (x(i),x(i+1)) x (y(j),y(j+1)).
    
    wNow = w(i, j) * sqrt(xScale(i)) * sqrt(yScale(j));
    
    if wNow > w_eps
        % Adjust width of cell to keep slope constant
        gaPrimeInv = (wNow / w0)^2;
        
        if crossFromRight
            xTilde(i+1) = xTilde(i) + (x(i+1) - x(i)) * gaPrimeInv;
            xScale(i) = 1 ./ gaPrimeInv;
        else
            yTilde(j+1) = yTilde(j) + (y(j+1) - y(j)) * gaPrimeInv;
            yScale(j) = 1 ./ gaPrimeInv;
        end

        % Find slope that would get us to the corner [i+1,j+1]
        alCorner = atan2(yTilde(j+1) - yNow, xTilde(i+1) - xNow);

        % Determine how we cross into next cell
        if abs(alCorner - al0) < eps
        %if abs(yNext - yTilde(j+1)) < tol
            % disp('Luck.');
            segmentFound = 1;
            iend = i+1;
            jend = j+1;
            
            numPts = numPts + 1;
            crossPtsTilde(1, numPts) = xTilde(i+1);
            crossPtsTilde(2, numPts) = yTilde(j+1);
            crossType(numPts) = 1;
            
            break
        elseif al0 < alCorner
            crossToRight = 1;
            xNext = xTilde(i+1);
            yNext = yNow +(xNext - xNow) * tan(al0);
        else
            crossToRight = 0;
            yNext = yTilde(j+1);
            xNext = xNow + (yNext - yNow) * tan(pi/2 - al0);
        end
    elseif crossFromRight
        % If slope is zero we collapse the interval in the new
        % coordinates.
        xTilde(i+1) = xTilde(i);
        xScale(i) = Inf;

        xNext = xTilde(i+1);
        yNext = yNow;            
        crossToRight = 1;
    else % crossFromRight = 0
        % If slope is infinity we collapse the interval.
        yTilde(j+1) = yTilde(j);
        yScale(j) = Inf;
    
        xNext = xNow;
        yNext = yTilde(j+1);
        crossToRight = 0;
    end
    
    % Perform the actual crossing
    if crossToRight
        if wNow > 0
            % Add vertex above crossing point to list of candidates for
            % endpoints of closest P-segment
            iAbove(numAbove + 1) = i + 1;
            jAbove(numAbove + 1) = j + 1;
            indAbove(numAbove + 1) = numPts + 1;
            numAbove = numAbove + 1;
        end
        
        i = i + 1;
    else
        % When we cross above there is no candidate for a P-segment.
        j = j + 1;
    end
    
    % Update new point
    xNow = xNext;
    yNow = yNext;
    crossFromRight = crossToRight;
    
    numPts = numPts + 1;
    crossPtsTilde(1, numPts) = xNow;
    crossPtsTilde(2, numPts) = yNow;
    crossType(numPts) = crossToRight;
end

% For debugging: in the new coordinates following the P-segment with the
% original slope h0 should give a straight line.
%
% % The new coordinates won't span the whole grid and we have to add one 
% % one more boundary point.
% if i < m
%     i = i + 1;
% end
% if j < n
%     j = j + 1;
% end
% 
% % Perform rescaling to q1, q2 as effected by the change of parametrization
% un = [u .* xScale; Inf];
% q1n = [xTilde, un];
% q1n = q1n(i0:i,:); % Take only those we have looked at
% 
% vn = [v .* yScale; Inf];
% q2n = [yTilde, vn];
% q2n = q2n(j0:j,:);
% 
% % Now follow the segment
% gamma2 = srvf_klr_follow_p_seg(1, 1, xTilde(i0:i,:), yTilde(j0:j,:), ...
%                                w (????), h, tol);
% 
% % This should give a straight line
% figure;
% plotGamma(xTilde(i0:i), yTilde(j0:j), gamma2);

numPtsSeg = 0;

if segmentFound
    al = al0;
    numPtsSeg = numPts;
elseif numAbove == 0
    % In this case there is no P-segment with slope larger than h0.
    al = Inf;
else
    % Find grid point closest to the constructed P-segment
    iAbove = iAbove(1:numAbove);
    jAbove = jAbove(1:numAbove);
    cornerSlopes = atan2(yTilde(jAbove) - yTilde(j0), ...
                             xTilde(iAbove) - xTilde(i0));

    % We don't want the mathematical minimum. If two slopes are very close
    % together we prefer the corner that is closer to our point, hence we
    % require a new candidate for the minimum to be at least al_eps
    % smaller.
    indMin = 1;
    for k = 2:numAbove
        if cornerSlopes(k) < cornerSlopes(indMin) - al_eps
            indMin = k;
        end
    end

    % Slope for P-segment
    iend = iAbove(indMin);
    jend = jAbove(indMin);
    numPtsSeg = indAbove(indMin);
    al = atan2(yTilde(jend) - yTilde(j0), xTilde(iend) - xTilde(i0));
end

% Adjust crossing points to new angle
crossPts = zeros(2, numPtsSeg);
for k = 1:numPtsSeg
    if crossType(k) % Cross to right --> adjust y coordinate
        crossPts(1, k) = crossPtsTilde(1, k);
        crossPts(2, k) = yTilde(j0) + ...
                            (crossPtsTilde(1, k) - xTilde(i0)) * tan(al);
    else % Cross to top --> adjust x coordinate
        crossPts(1, k) = xTilde(i0) + ...
                            (crossPtsTilde(2, k) - yTilde(j0)) * tan(pi/2-al);
        crossPts(2, k) = crossPtsTilde(2, k);
    end
end

if isinf(al)
    gamma = [];
    gammaInd = [];
    return
end

% Adjust scaling back to original values
i = i0;
j = j0;
gamma = zeros(2, numPtsSeg+1);
gammaInd = zeros(2, numPtsSeg+1);
gamma(:,1) = [x(i0); y(j0)];
gammaInd(:, 1) = [i0; j0];

for k = 1:numPtsSeg
    if ~isinf(xScale(i)) && ~isinf(yScale(j))
        gamma(1, k+1) = x(i) + xScale(i) * (crossPts(1, k) - xTilde(i));
        gamma(2, k+1) = y(j) + yScale(j) * (crossPts(2, k) - yTilde(j));
    elseif isinf(xScale(i))
        gamma(1, k+1) = x(i+1);
        gamma(2, k+1) = gamma(2, k);
    else
        gamma(1, k+1) = gamma(1, k);
        gamma(2, k+1) = y(j+1);
    end
    
    gamma(1,k+1) = max(gamma(1,k), gamma(1,k+1));
    gamma(2,k+1) = max(gamma(2,k), gamma(2,k+1));
    
    if crossType(k)
        i = i + 1;
    else
        j = j + 1;
    end
    
    gammaInd(:, k+1) = [i; j];
end

gammaInd(:, end) = [iend; jend];

assert(min(diff(gamma(1,:))) >= 0);
assert(min(diff(gamma(2,:))) >= 0);

end