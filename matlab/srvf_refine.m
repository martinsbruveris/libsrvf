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


% Given a piecewise-constant SRVF Q with changepoint vector T, 
% and a new changepoint vector Tr which is a refinement of T, 
% compute the values qr corresponding to Tr.  Basically, this 
% just amounts to duplicating an element of Q whenever the corresponding 
% interval is split.
%
% Inputs
%  Q,T : the SRVF
%  Tr : a refinement of T
%
% Outputs
%  Qr : the SRVF values corresponding to the new changepoint vector
% --------------------------------------------------------------------------
function Qr = srvf_refine( Q, T, Tr )
    assert( all(ismember(T, Tr)) ); % Note, this check might be slow.
    assert( min(diff(T)) > 0 );
    assert( min(diff(Tr)) > 0 );

    tv = (Tr(1:end-1)+Tr(2:end))/2;
    Qr = srvf_evaluate(Q,T,tv);
end


%!test
%! Q=[1];
%! T=[0 1];
%! Tr=[0 1/4 1/2 3/4 1];
%! Qrexp=[1 1 1 1];
%! Qr=srvf_refine( Q, T, Tr );
%! assert( Qr, Qrexp );
%!
%!test
%! Q=[1 -1 1 -1 1;\
%!    1 -1 1 -1 1];
%! T=[0 1/5 2/5 3/5 4/5 1];
%! Tr=unique([T 0.19 0.2001 0.65 0.999]);
%! Qrexp=[1 1 -1 -1 1 -1 -1 1 1;\
%!        1 1 -1 -1 1 -1 -1 1 1];
%! Qr=srvf_refine( Q, T, Tr );
%! assert( Qr, Qrexp );
