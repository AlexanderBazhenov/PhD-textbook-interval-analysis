## Created: 2024-05-24
function retval = Dist1AB22 (X1, X0)
% 2x2 matrix distance norm1
retval = max( abs(inf(X1(1,1))-inf(X0(1,1))),  abs(sup(X1(1,1))-sup(X0(1,1)))) + ...
max( abs(inf(X1(1,2))-inf(X0(1,2))),  abs(sup(X1(1,2))-sup(X0(1,2)))) + ...
max( abs(inf(X1(2,1))-inf(X0(2,1))),  abs(sup(X1(2,1))-sup(X0(2,1)))) + ...
max( abs(inf(X1(2,2))-inf(X0(2,2))),  abs(sup(X1(2,2))-sup(X0(2,2))));
endfunction
