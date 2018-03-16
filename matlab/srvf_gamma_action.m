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


% Computes the action of a 1-D diffeomorphism on the given SRVF.
%
% Inputs:
%  Q  - SRVF function values
%  TQ - SRVF parameter values
%  G  - diffeomorphism function values.  Must be nondecreasing, and must 
%       satisfy TQ(1)<=G(1) and TQ(end)>=G(end).
%  TG - diffeomorphism parameter values.
%
% Outputs:
%  Qr - the new SRVF function values
%  Tr - the new SRVF parameter values
% --------------------------------------------------------------------------
function [Qr, Tr] = srvf_gamma_action(Q, TQ, G, TG)
    % Construct finer subdivision where Q*G will be piecewise constant
    TGx = plf_preimages(G,TG,TQ);
    Tr = unique([TG TGx]);
    
    % Evaluate Q on G(Tr)
    Trmid = (Tr(1:end-1) + Tr(2:end)) / 2;
    Grmid = plf_evaluate(G, TG, Trmid);
    Qr = srvf_evaluate(Q, TQ, Grmid);
    
    % Multiply by sqrt(G')
    Gr = plf_evaluate(G, TG, Tr);
    DGr = (Gr(2:end)-Gr(1:(end-1))) ./ (Tr(2:end)-Tr(1:(end-1)));
    Qr = Qr .* repmat(sqrt(DGr), size(Qr,1), 1);
end

%!test
%! Q=[1];
%! TQ=[0 1];
%! G=[0 1/3 1];
%! TG=[0 1/2 1];
%! [Qr Tr]=srvf_gamma_action(Q,TQ,G,TG);
%! Qrexp=[sqrt(2/3) sqrt(4/3)];
%! Trexp=[0 1/2 1];
%! assert(Qr,Qrexp,1e-6);
%! assert(Tr,Trexp,1e-6);
%!
%!test
%! Q=[1 -1];
%! TQ=[0 1/2 1];
%! G=[0 1/3 1];
%! TG=[0 2/3 1];
%! [Qr Tr]=srvf_gamma_action(Q,TQ,G,TG);
%! Qrexp=[sqrt(1/2) sqrt(2) -sqrt(2)];
%! Trexp=[0 2/3 3/4 1];
%! assert(Qr,Qrexp,1e-6);
%! assert(Tr,Trexp,1e-6);
