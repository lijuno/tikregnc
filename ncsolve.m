function varargout = ncsolve(A, b, lambda, L)
% Examples: 
%    X = ncsolve(A, b, lambda, L);
%    [X, res, smn] = ncsolve(A, b, lambda, L);
% Solve the linear least square program with x>=0 constraint
% min ||C*x -d|| subject to A*x<=b, A_eq*x = b_eq, lb<=x<=ub
% Input args: 
%    A: Discretization matrix
%    b: data (the transient current)
%    lambda: the regularization parameter
%    L: derivative matrix
% (Optional) output args: 
%    X: the solution
%    rmn: fitting residue
%    smn: smoothness

x0 = zeros(size(L, 2),1);  % default 
N = size(A, 2);
C = [A; lambda*L];
d = [b; lambda* L*x0];
AA = -eye(N, N);
bb = zeros(N, 1);
opts = optimset('MaxIter', 2000);
[X, resnorm, residual, exitflag, output] = lsqlin(C, d, AA, bb, [], [], [],[], [], opts); 
res = (A*X - b)'*(A*X -b);   % fitting residue
smn =(L*X)'*(L*X);     % smoothness

if nargout >=1
    varargout=cell(nargout,1);
    varargout{1} = X; 
    if nargout >=2
       varargout{2} = res; 
       if nargout >=3
           varargout{3} = smn; 
       end
    end
end
