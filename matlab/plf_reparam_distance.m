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
%
% Computes the distance between F1 and F2 after optimizing over
% reparametrizations.
%
% Inputs:
%  F1, T1 : a plf
%  F2, T2 : a plf
%  method : one of {'dp', 'robinson', 'klr'}
%           Note: 'robinson only works on 1D data
%
% Outputs
%  d : the distance between the equivalence classes of F1 and F2
% -------------------------------------------------------------------------
function d = plf_reparam_distance(F1, T1, F2, T2, method)
    if nargin < 5
        method = 'dp';
    end
    
    % Checking input arguments
    assert( size(F1,1) == size(F2,1) );
    assert( size(T1,1) == 1 && size(T1,2) == size(F1,2) );
    assert( size(T2,1) == 1 && size(T2,2) == size(F2,2) );
    assert( min(diff(T1)) > 0 );
    assert( min(diff(T2)) > 0 );
    assert( abs(T1(1)-T2(1)) < 1e-12 );
    assert( abs(T1(end)-T2(end)) < 1e-12 );
    % Robinson's method only works on functional data
    assert( ~(strcmp(method, 'robinson') && size(F1, 1) > 1) );
    
    if strcmp(method, 'klr') % Exact matching following the paper by
                             % Klassen, Lahiri and Robinson
        % Compute SRVFs
        Q1 = plf_to_srvf(F1, T1);
        Q2 = plf_to_srvf(F2, T2);
        
        % Find optimal matching
        [G1, TG1, G2, TG2] = ...
            srvf_klr_optimal_reparam(Q1, T1, Q2, T2);
        
        % Apply optimal matching
        [Q1r, TQ1r] = srvf_gamma_action(Q1, T1, G1, TG1);
        [Q2r, TQ2r] = srvf_gamma_action(Q2, T2, G2, TG2);
        
        % Compute L^2-distance
        d = srvf_l2distance(Q1r, TQ1r, Q2r, TQ2r);
        
    elseif strcmp(method, 'robinson') % Exact matching for 1D data as 
                                      % described in Robinsons' thesis
        % srvf_fa_optimal_reparam assumes that the functions are 
        % scaled to unit arclength, and have constant-speed parametrizations
        L1 = plf_arclength(F1, T1);
        L2 = plf_arclength(F2, T2);
        F1s = F1 ./ L1;
        F2s = F2 ./ L2;

        [Gn1, TGn1] = plf_constant_speed_reparam(F1s, T1);
        [Gn2, TGn2] = plf_constant_speed_reparam(F2s, T2);

        [F1n, TF1n] = plf_compose(F1s, T1, Gn1, TGn1);
        [F2n, TF2n] = plf_compose(F2s, T2, Gn2, TGn2);

        Q1 = plf_to_srvf(F1n, TF1n);
        Q2 = plf_to_srvf(F2n, TF2n);

        % srvf_fa_optimal_reparam also assumes that the SRVF values
        % alternate between 1 and -1 on adjacent intervals
        [Q1n, TQ1n] = srvf_make_alternating(Q1, TF1n);
        [Q2n, TQ2n] = srvf_make_alternating(Q2, TF2n);

        [G1, TG1, G2, TG2] = srvf_fa_optimal_reparam(Q1n, TQ1n, Q2n, TQ2n);
        
        % Apply optimal matching
        [Q1r, TQ1r] = srvf_gamma_action(Q1n, TQ1n, G1, TG1);
        [Q2r, TQ2r] = srvf_gamma_action(Q2n, TQ2n, G2, TG2);
        
        % Compute L^2-distance, remembering that we rescaled F1 and F2
        % Note that taking the SRVFs scale with sqrt(L).
        dsq = L1 * srvf_l2norm(Q1r, TQ1r)^2 + ...
                -2 * sqrt(L1) * sqrt(L2) * srvf_l2product(Q1r, TQ1r, Q2r, TQ2r) + ...
                L2 * srvf_l2norm(Q2r, TQ2r)^2;
        d = sqrt(dsq);
        
    elseif strcmp(method, 'dp') % Dynamic programming
        % Compute SRVFs
        Q1 = plf_to_srvf(F1, T1);
        Q2 = plf_to_srvf(F2, T2);
        
        % Find optimal matching using dynamic programming
        [G, TG] = srvf_optimal_reparam(F1, T1, F2, T2);
        
        % Apply optimal matching
        [Q2r, TQ2r] = srvf_gamma_action(Q2, T2, G, TG);
        
        % Compute L^2-distance
        d = srvf_l2distance(Q1, T1, Q2r, TQ2r);
    else
        error('Unknown method.');
    end

end