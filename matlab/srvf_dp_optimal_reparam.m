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
% the dynamical programming algorithm. Pure matlab implementation in case
% compiling the C++ code creates difficulties.
%
% Inputs
%  F1, T1: The first function
%  F2, T2: The second function
%
% Outputs
%  G, TG: Reparametrization of F2
% -------------------------------------------------------------------------
function [G, TG] = srvf_dp_optimal_reparam(F1, T1, F2, T2)

% Hard-coded neighborhood for faster access
nbrs = [1, 1;  1, 2;  2, 1;  2, 3;  3, 2;  1, 3;  3, 1; ...
        1, 4;  3, 4;  4, 3;  4, 1;  1, 5;  2, 5;  3, 5; ...
        4, 5;  5, 4;  5, 3;  5, 2;  5, 1;  1, 6;  5, 6; ...
        6, 5;  6, 1];
numNbrs = size(nbrs, 1);

% Compute SRVFs
Q1 = plf_to_srvf(F1, T1);
Q2 = plf_to_srvf(F2, T2);

% Rename input variables to something more familiar
x = T1;
u = Q1;
y = T2;
v = Q2;

m = length(x);
n = length(y);

% Weight matrix
w = u' * v;

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

% Keeping track of where the last segment ending at (i,j) originated.
prevVertex = zeros(1, numVertices);

for k = 1:numVertices    
    
    i = vertexList(1, k);
    j = vertexList(2, k);
    
    % We cannot reach this vertex for some reason
    if isinf(distOpt(i,j))
        continue
    end
    
    for l = 1:numNbrs
        in = i + nbrs(l, 1);
        jn = j + nbrs(l, 2);
        
        % We don't want to run out of bounds
        if in > m || jn > n
            continue
        end
        
        E = srvf_dp_l2prod(Q1, T1, Q2, T2, i, j, in, jn);
        
        % We have found a new maximum
        if distOpt(i, j) + E > distOpt(in, jn)
            distOpt(in, jn) = distOpt(i, j) + E;
            
            kn = pair2list(in, jn);
            prevVertex(kn) = k;
        end
    end
end

% Backtracking to construct reparametrization
segSeq = zeros(1, m);
segSeq(1) = pair2list(m, n);
numSeg = 1;
while segSeq(numSeg) > 1
    segSeq(numSeg+1) = prevVertex(segSeq(numSeg));
    numSeg = numSeg + 1;
end

% Reverse order
segSeq = flip(segSeq(1:numSeg));

% Indices corresponding to gamma
ig = vertexList(1, segSeq);
jg = vertexList(2, segSeq);

% Construct piecewise linear gamma
TG = x(ig);
G = y(jg);
