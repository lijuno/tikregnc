function varargout = tikregnc(t, b, dconfig, rconfig) 
% A rewritten solver of Tikhonov regularization with non-negativity constraint
% Usage: [X] = tikregnc(t, b, dstruc, rconfig) 
%        [X, f] = tikregnc(t, b, dstruc, rconfig)
%        [X, f, A, res, smn] = tikregnc(t, b, dstruc, rconfig)
% Input arguments:
%   t: time axis, a column vector
%   b: observed data, a column vector
%   dconfig: discretization configuration, a struct
%       dconfig.type: 'linear', 'log' or 'glq'
%       dconfig.N: discretization points in f space
%       dconfig.param1: [lower_bound, upper_bound] for 'linear' or 'log' type; 
%                       or scaling factor (a scalar) for 'glq' type
%   rconfig: regularization configuration, a struct
%         rconfig.type: 'External', 'Morozov' (Morozov principle), 'L-curve1' (L-curve criterion)
%         rconfig.order: L_order (the order of L matrix
%         rconfig.param1: depending on rconfig.type as follows
%             'External' type: regularzation parameter (a scalar)
%             or 'Morozov' type: delta_b (noise of b, a scalar) 
%             or 'L-curve*' type: l_curve_plot, whether to plot L-curve (1 for yes, 0 for no); default: 0
%   
% Output arguments:
%   X: regularized solutions 
%   f_out: discretized f space
%   A: discretization matrix
%   lambda_out: regularization parameter used in the end
% 
% Written by Lijun Li, Oct, 2013
% 
% Update log: 


glq_alpha_ratio_tol = 1e-4;  % The tolerence of searching for the alpha value in 'glq'

% Get discretization matrix A and f-space axis f
if strcmpi(dconfig.type, 'linear') || strcmpi(dconfig.type, 'log')
    [A, f] = get_A(t, dconfig.type, dconfig.N, dconfig.param1);
elseif strcmpi(dconfig.type, 'glq')
    sf = dconfig.param1;
    [A, f] = get_A(t*sf, dconfig.type, dconfig.N);  % scale the t axis
end

% Solve the problem
L = full(get_l(length(f), rconfig.order));
if strcmpi(rconfig.type, 'External')
    % Use the externally provided regularization parameter
    [sol] = ncsolve(A, b, rconfig.param1, L);
    
elseif strcmpi(rconfig.type, 'Morozov')
    % Use Morozov discrepancy principle to determine the regularization
    % parameter
    delta_b = rconfig.param1;  % noise amplitude of observed data b
    lambda_opt = fsolve(@(lambda) lambda_func(lambda, A, b, L, delta_b), 1);
    [sol] = ncsolve(A, b, lambda_opt, L);
    
elseif strcmpi(rconfig.type, 'L-curve1')
    % Use L-curve criterion to determine the regularization parameter
    
    l_curve_plot = rconfig.param1; 
    [U, sm, XX, V] = cgsvd(A, L); % Generalized singular value decomposition in compact form
    figure(1)
    [lc0, rho, eta, reg_params] = l_curve(U,sm,b, 'Tikh');
    % lc0 is the regularization parameter from the  unconstrained problem, 
    % which can be used as the initial value for the search of 
    % the regularization parameter for the constrained problem
    close all
    
    X0 = ncsolve(A, b, lc0, L);   % initial value 
    dlambda = lc0 * 1e-2;
    res0 = sum((A*X0 - b).^2);
    smn0 = sum((L*X0).^2);
    alpha = lc0^2;
    dalpha_ratio = (res0/smn0 - lc0^2)/lc0^2;
    disp(['Init alpha = ', num2str(lc0^2)])
    iter = 1;
    while (abs(dalpha_ratio) > glq_alpha_ratio_tol)
        alpha1_arr(iter) = alpha;  % will become a row vector
        [X1, res1, smn1] = ncsolve(A, b, sqrt(alpha), L);
        dalpha_ratio = (alpha - res1/smn1)/alpha;
        alpha = res1/smn1;
        res1_arr(iter) = res1;
        smn1_arr(iter) = smn1;
        disp(['iteration=', num2str(iter), '; alpha=', num2str(alpha)])
        iter = iter + 1;
    end
    lambda_opt = sqrt(alpha);
    sol = ncsolve(A, b, lambda_opt, L);
    
    if l_curve_plot
        l_curve_lambda2 =lambda_opt*logspace(-2, log10(20), 10);
        res2_arr = zeros(1, length(l_curve_lambda2));
        smn2_arr = zeros(1, length(l_curve_lambda2));
        for ii=1:length(l_curve_lambda2)
            [Xtmp, res2_arr(ii), smn2_arr(ii)] = ncsolve(A, b, l_curve_lambda2(ii), L);
        end
        
        res_arr = [res1_arr, res2_arr];
        smn_arr = [smn1_arr, smn2_arr];
        [l_curve_lambda, ind] = sort([sqrt(alpha1_arr), l_curve_lambda2]);
        figure(10)
        loglog(res_arr(ind), smn_arr(ind), '-', res1_arr(end), smn1_arr(end), 'o', 'linewidth', 2)
        axis_min = min(min(res_arr), min(smn_arr));
        axis_max = max(max(res_arr), max(smn_arr));
        axis([axis_min axis_max axis_min axis_max])
        xlabel('Log(||C*x-d||^2)', 'fontsize', 16)
        ylabel('Log(||L*x||^2)', 'fontsize', 16)
        title('L-curve', 'fontsize', 12)
    end
else 
    error('Regularization type unrecognized!')
end

if strcmp(dconfig.type, 'glq')
    % If GLQ discretization is used, scale the solution and f output
    sol = sol/sf^2;
    f = f*sf;   
end

if nargout >=1
    varargout=cell(nargout,1);
    varargout{1} = sol; 
    if nargout >=2
        varargout{2} = f;
        if nargout >=3
            varargout{3} = A;
            if nargout >=4
                varargout{4} = lambda_opt;
            end
        end
    end
end

end

function y = lambda_func(lambda, A, b, L, delta_b)
% The affiliated function when solving ||A*x -b|| = delta_b when using
% Morozov discrepancy principle to select regularization parameters
%
% Input args: 
% lambda: regularization parameter
% A, b, L: same as those defined in tikregnc.m
% delta_b: noise amplitude in b, a scalar
% 
x = ncsolve(A, b, lambda, L);
y = (A*x-b)'*(A*x-b) - delta_b^2;
end
