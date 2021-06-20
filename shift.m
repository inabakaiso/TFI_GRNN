function y = shift(x, nsr, nsc)
%******************************************************************
%SHIFT Shift indices of a vector or matrix.    ====indice‚Æ‚ÍõˆøH====
%   Y = SHIFT(X, NSR) returns the vector with indices shifted by NSR.
%
%   Y = SHIFT(X, NSR, NSC) returns the matrix with row indices shifted
%   by NSR and column indices shifted by NSC.
%
%   NSR and NSC can be both positive and negative.
%
%   ex1: x = [1:9]
%       shift(x, 2)  = [3 4 5 6 7 8 9 0 0]
%       shift(x, -2) = [0 0 1 2 3 4 5 6 7]
%
%   ex2: x = [1:3;1:3];
%       shift(x,1,-1) = [0 1 2; 0 0 0]
%
% Copyright (c) Chang-Jun Ahn 2000 (jun@sasase.ics.keio.ac.jp)
%*******************************************************************

if nargin == 1
     nsr = 0;
     nsc = 0;
end;

if nargin == 2
     nsc = 0;
end;

[r,c] = size(x);
if nsr > 0
     if r == 1
         y = [x((nsr+1):c) zeros(r,nsr)];
     elseif c == 1
         y = [x((nsr+1):r); zeros(nsr,c)];
     elseif nsc > 0
         y = [[x((nsr+1):r,(nsc+1):c); zeros(nsr,(c-nsc))] zeros(r,nsc)];
     else
         y = [zeros(r,abs(nsc)) [x((nsr+1):r,1:(c+nsc)); zeros(nsr,(c+nsc))]];
     end;
else
     if r == 1
         y = [zeros(r,abs(nsr)) x(1:(c+nsr))];
     elseif c == 1
         y = [zeros(abs(nsr),c); x(1:(r+nsr))];
     elseif nsc > 0
         y = [[zeros(abs(nsr),(c-nsc)); x(1:(r+nsr),(nsc+1):c)] zeros(r,nsc)];
     else
         y = [zeros(r,abs(nsc)) [zeros(abs(nsr),(c+nsc)); 
x(1:(r+nsr),1:(c+nsc))]];
     end;
end;