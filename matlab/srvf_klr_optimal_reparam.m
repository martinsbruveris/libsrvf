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

% Computes the optimal reparametrization between the SRVFs Q1 and Q2 using
% the algorithm described in
%   Sayani Lahiri, Daniel Robinson and Eric Klassen. Precise matching of PL
%   curves in R^N in the square root velocity framework. Geometry, Imaging
%   and Computing, 2(3), 133-186, 2015.
%
% Inputs
%  Q1, T1: The first SRVF
%  Q2, T2: The second SRVF
%
% Outputs
%  G1, GT1: Reparametrization of Q1
%  G2, GT2: Reparametrization of Q2
%  segPts:  (2,numPts) matrix containing the endpoints of the P-segments
%           and N-segments as described in the algorithm.
%  distOpt: Matrix of inner products of the optimal paths from point (1,1) 
%           to points (i,j) as built up during dynamic programming.
% -------------------------------------------------------------------------

function [G1, GT1, G2, GT2, segPts, distOpt] = ...
            srvf_klr_optimal_reparam(Q1, T1, Q2, T2)

% Some flags to change behaviour of algorithm
debug_mode = 0;
ignore_slope_constraints = 0;

% al_eps is the increment of alpha after a P-segment has been found
al_eps = 1e-12;
% w_eps is the threshold below which entries in the weight matrix w=q1'*q2
% will be counted as negative
w_eps = 1e-6;

% Rename input variables to something more familiar
x = T1;
u = Q1;
y = T2;
v = Q2;

m = length(x);
n = length(y);

% Weight matrix
w = u' * v;
% Very small weights lead to problems when looking for P-segments because of
% the rescaling. This is covered by the parameter w_eps.

% The order in which we iterate over vertices
xList = repmat(1:m, [n, 1])';
yList = repmat(1:n, [m, 1]);
vertexList = [xList(:)'; yList(:)']; % size=(2, n*m)
numVertices = size(vertexList, 2);
pair2list = @(a, b) m*(b-1) + a; % Mapping between (i,j) and vertexList

% Matrix saving the inner product of the optimal path from point (1,1) to
% point (i,j)
distOpt = -Inf * ones(m, n);
distOpt(1, 1) = 0;

% Keeping track where the last P-segment of the optimal path ending at
% vertex (i,j) ended and what the final slope of this vertex is. This is
% used for the slope constraints between two consecutive P-segments.
prevPSegEnd = zeros(1, numVertices);
prevPSegEndSlope = zeros(2, numVertices);

% Keeping track of where the last segment ending at (i,j) originated and
% what its initial slope is. Slope=Inf for N-segments. Used for
% constructing optimal matching at the end.
prevVertex = zeros(1, numVertices);
prevSlope = zeros(1, numVertices);

for k = 1:numVertices    
    % if mod(k, m) == 0
    %     disp(k);
    % end

    i = vertexList(1, k);
    j = vertexList(2, k);
    
    % No P- or N-segment has terminated at this vertex. We can skip it.
    if isinf(distOpt(i,j))
        continue
    end
    
    if i < m && j < n && w(i,j) > w_eps
        % We are looking for P-segments.
        
        % We are at the boundary already. No P-segments any more.
        if i == m || j == n
            continue
        end
        
        % Find slope constraints on new P-segment from this vertex
        if prevPSegEnd(k) == 0 || prevPSegEnd(k) == 1
            % We don't know where or how the previous P-segment ended.
            % Hence, no constraints on slope.
            almin = al_eps;
            almax = pi/2;
        else
            % This is where the previous P-segment ended and its slope
            kp = prevPSegEnd(k);
            ip = vertexList(1, kp);
            jp = vertexList(2, kp);
            dxp = prevPSegEndSlope(1,k);
            dyp = prevPSegEndSlope(2,k);
            
            % Calculate constants from the paper
            A = w(ip-1, jp-1);
            B = w(i,j);
            C = w(ip-1, j);
            D = w(i, jp-1);
            
            if debug_mode
                assert(A > w_eps);
                assert(B > w_eps);
            end
            
            % In this case there are no optimal P-segments from here
            if C > 0 && D > 0 && C*D >= A*B + 1e-12
               continue
            end
            
            % Use constraints on the change of slope from the paper
            al_eps2 = sqrt(al_eps);
            if C > w_eps && D > w_eps
                if C*D - A*B < 1e-12
                    % If C*D = A*B, then the interval collapses to a point.
                    % We need to numerically stabilize this operation.
                    % We use al_eps2 here to make the interval large
                    % enough.
                    almin = atan2(D^2 * dyp, C^2 * dxp) - al_eps2;
                    almax = almin + al_eps2;
                else
                    almin = atan2(D^4 * dyp, A^2*B^2 * dxp) - al_eps;
                    almax = atan2(A^2*B^2 * dyp, C^4 * dxp) + al_eps;
                end
            elseif D <= w_eps &&  C > w_eps
                almin = al_eps;
                almax = atan2(A^2*B^2 * dyp, C^4 * dxp) + al_eps;
            elseif D > w_eps && C <= w_eps
                almin = atan2(D^4 * dyp, A^2*B^2 * dxp) - al_eps;
                almax = pi/2;
            else
                almin = al_eps;
                almax = pi/2;
            end
        end % End of looking at the history for [almin, almax]
            
        if debug_mode
            assert(almin < almax);
        end
        
        % Set initial and final slope for searching.
        al = max(almin, al_eps);
        almax = min(pi/2, almax);
        
        % In case using the constraints from the paper leads to problems
        if ignore_slope_constraints
            al = almin;
            almax = pi/2;
        end
            
        % Iterate over all permissible P-segments
        segmentsLeft = 1;
        while segmentsLeft           
            [alNext, gammaNext, gammaNextInd] = ...
                srvf_klr_find_next_p_seg(i, j, x, y, w, al, al_eps, w_eps);
            
            if debug_mode && al > al_eps && al < almax && ...
                    alNext <= al
                disp('When searching for P-seg, alNext=al.');
            end
            
            % Have we found a valid segment?
            if alNext >= almax
                segmentsLeft = 0;
                continue;
            end
            
            in = gammaNextInd(1, end);
            jn = gammaNextInd(2, end);
            kn = pair2list(in, jn);
            
            % Inner product of newly found segment
            ESeg = srvf_klr_l2prod_gamma(u, v, ...
                gammaNext(1,:), gammaNext(2,:), gammaNextInd(1,:), gammaNextInd(2,:));
            
            % if i == 10 && j == 8 && in == 23 && jn == 20
            % disp(['(', num2str(i), ', ', num2str(j), ') ', ...
            %       '(', num2str(in), ', ', num2str(jn), ') ', ...
            %       num2str(ESeg), ' ', num2str(alNext), ' ', num2str(al)]);
            % end
            
            % We have found a new maximum
            if distOpt(i, j) + ESeg > distOpt(in, jn)
                distOpt(in, jn) = distOpt(i, j) + ESeg;
                prevPSegEnd(kn) = kn;
                
                dy = gammaNext(2,end) - gammaNext(2,end-1);
                dx = gammaNext(1,end) - gammaNext(1,end-1);
                prevPSegEndSlope(:, kn) = [dy; dx];
                
                prevVertex(kn) = k;
                prevSlope(kn) = alNext;
            end
            
            al = alNext + al_eps;
            if al > almax || al >= pi/2
                segmentsLeft = 0;
            end
        end
    else
        % We are looking for N-segments.
        
        % No two consecutive N-segments are allowed
        % Except that we always allow an N-segment to start from (1,1)
        if k ~= 1 && prevPSegEnd(k) ~= k
            continue;
        end
        
        % Construct all possible endpoints
        endPoints = srvf_klr_find_n_seg(i, j, w, w_eps);
        numPts = size(endPoints, 1);
        
        % Visit and update endpoints
        for l = 1:numPts
            i1 = endPoints(l, 1);
            j1 = endPoints(l, 2);
            k1 = pair2list(i1, j1);
            
            if distOpt(i, j) > distOpt(i1, j1)
                % N-segments have zero cost
                distOpt(i1, j1) = distOpt(i, j);
                prevPSegEnd(k1) = k;
                prevPSegEndSlope(:, k1) = prevPSegEndSlope(:,k);
                
                prevVertex(k1) = k;
                prevSlope(k1) = Inf;
            end
        end
    end
end

% Now we perform backtracking to build gamma
segSeq = zeros(1, m+n-1);
segSlope = zeros(1, m+n-1);
segSeq(1) = pair2list(m, n);
numSeg = 1;
while segSeq(numSeg) > 1
    segSeq(numSeg+1) = prevVertex(segSeq(numSeg));
    segSlope(numSeg+1) = prevSlope(segSeq(numSeg));
    numSeg = numSeg + 1;
end

% Reverse order
segSeq = flip(segSeq(1:numSeg));
segSlope = flip(segSlope(1:numSeg));

% Initialize gamma
gamma = zeros(2, m+n-1);
gamma(:, 1) = [x(1); y(1)];
numGamma = 1;

segPts = zeros(2, numSeg);
segPts(:,1) = [x(1); y(1)];

for k = 1:numSeg-1
    % Start and end points of segment
    i0 = vertexList(1, segSeq(k));
    j0 = vertexList(2, segSeq(k));
    i1 = vertexList(1, segSeq(k+1));
    j1 = vertexList(2, segSeq(k+1));
    segPts(:, k+1) = [x(i1); y(j1)];
    
    al = segSlope(k);
    
    if ~isinf(al)
        % This is a P-segment
        [gammaSeg, ~] = ...
            srvf_klr_follow_p_seg(i0, j0, x, y, w, al, i1, j1, w_eps);
        numGammaSeg = size(gammaSeg, 2) - 1;
    else
        % This is an N-segment
        [gammaSeg, ~] = srvf_klr_follow_n_seg(i0, j0, x, y, i1, j1);
        numGammaSeg = size(gammaSeg, 2) - 1;
    end
    
    % We skip the first point because it coincides with last point
    gamma(:,numGamma+1:numGamma+numGammaSeg) = gammaSeg(:, 2:end);
    numGamma = numGamma + numGammaSeg;
end

gamma = gamma(:,1:numGamma);
GT = linspace(x(1), x(end), numGamma);

G1 = gamma(1,:);
GT1 = GT;
G2 = gamma(2,:);
GT2 = GT;

end