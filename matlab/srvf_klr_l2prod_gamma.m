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

% Helper function for srvf_klr_optimal_reparam. Calculates the L^2-inner
% product between Q1 * G1 and Q2 * G2 with the SRVF action of G1, G2 on Q1,
% Q2.
%
% Inputs
%  Q1, Q2:     Values of the SRVFs
%  G1, G2:     Values of the reparamterization. Note that the inner product
%              is independent of the parametrizations TG1, TG2. Only the
%              path taken by G1, G2 is important.
%  G1ind, G2ind: To speed up evaluation. Q1(G1ind(k)) is the value of Q1
%              evaluated at G1(k).
%
% Outputs
%  E:  Value of the L^2-inner product
% -------------------------------------------------------------------------
function E = srvf_klr_l2prod_gamma(Q1, Q2, G1, G2, G1ind, G2ind)

G1p = diff(G1);
G2p = diff(G2);
Q1g = Q1(:, G1ind(1:end-1));
Q2g = Q2(:, G2ind(1:end-1));
wg = sum(Q1g .* Q2g, 1);

E = sum( wg .* sqrt(G1p) .* sqrt(G2p) );

end