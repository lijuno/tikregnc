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
%       dconfig.param1: scaling factor (a scalar). A few words about
%           scaling -- scaling is done for A matrix during the calculation
%           for the purpose of adjusting the scale of regularization parameter 
%           because different value of regularization parameter tends to
%           give very different computing time. If A becauses A*sf, then
%           the corresponding soultion x becomes x/sf, and regularization
%           paramter will be sf times larger (since it is actually lambda^2
%           before the ||x||_2^2 term. 
%           update: seems like the scaling isn't really helping speed up
%           the program
%       dconfig.param2: [lower_bound, upper_bound] for 'linear' or 'log' type
%   rconfig: regularization configuration, a struct
%         rconfig.type: model selection type
%             'External': single regularization parameter provided externally
%             'Morozov': Morozov principle
%             'L-curve-ext': L-curve with regularization parameter provided externally
%             'OCV-ext': ordineary cross validation with lambda externally provided
%         rconfig.order: L_order (the order of L matrix
%         rconfig.param1: depending on rconfig.type as follows
%             for 'External' type: regularzation parameter (a scalar)
%             or for 'Morozov' type: delta_b (noise of b, a scalar) 
%             or for 'L-curve-ext' type: array of regularization parameter, and plot L-curve
%             or for 'L-curve1': initial guess of regularization parameter
%             or for 'OCV-ext' type: an array of externally provided regularization parameters, 
%                  and find the optimal lambda by fitting CV curve
%         rconfig.param2: 
%             for 'L-curve1': whether to plot L-curve 
%             for 'OCV-*': the n for leave-n-out CV, or the number of test data points 
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
%     Mar, 2014: Add ordinary cross-validation for model selection


glq_alpha_ratio_tol = 1e-4;  % The tolerence of searching for the alpha value in 'glq'

t = force_column_vector(t);
b = force_column_vector(b);

[A, f_orig] = get_A(t, dconfig.type, dconfig.N, dconfig.param2); 
% note that f_orig is the real f grid for 'log' scaling, not log(f) grid
sf = dconfig.param1;  % scaling factor
A = A*sf; 

% Solve the problem
L = full(get_l(length(f_orig), rconfig.order));
if strcmpi(rconfig.type, 'External')
    % Use the externally provided regularization parameter
    sol_scaled = ncsolve(A, b, rconfig.param1, L);
    lambda_opt = rconfig.param1;   % for the completeness of the output args
    
elseif strcmpi(rconfig.type, 'Morozov')
    % Use Morozov discrepancy principle to determine the regularization
    % parameter
    delta_b = rconfig.param1;  % noise amplitude of observed data b
    lambda_opt = fsolve(@(lambda) lambda_func(lambda, A, b, L, delta_b), 1);
    [sol] = ncsolve(A, b, lambda_opt, L);

elseif strcmpi(rconfig.type, 'OCV-ext')
    % Ordinary cross validation with regularization parameters externally
    % provided
    % Note that the calling times of l1_ls sove is length(lambda_arr) * k_fold
    
    lambda_arr = rconfig.param1;  % regularization parameter array provided externally
    n_testdata = rconfig.param2;   % the number of test data points, or the data point of each subsample
    k_fold = length(b) /n_testdata;    % k for k-fold ordinary cross-validation
    
    ocv_error_mean = zeros(length(lambda_arr), 1);
    ocv_error_std = zeros(length(lambda_arr), 1);
    ocv_error_iter = zeros(k_fold, 1);   % OCV error in each iteration below
    
    disp(['Going to do ', num2str(k_fold), '-fold cross-validation for ', num2str(length(lambda_arr)), ' regularization parameters'])
    disp(['Data points in each sub-sample = ', num2str(n_testdata)])
    disp(['Discretization scaling factor = ', num2str(dconfig.param1)])
    disp('--------------------------------')
    
    for ii = 1:length(lambda_arr)
        lambda_iter = lambda_arr(ii);
        disp(['Iteration #', num2str(ii), ', lambda = ', num2str(lambda_iter)])
        for jj = 1:k_fold
            
            % First compute the indices of training data and testing data
            ind_testdata = zeros(length(n_testdata), 1);
            for kk = 1:n_testdata
                ind_testdata(kk) = jj + k_fold*(kk-1);
            end
            ind_whole = 1:length(b);
            ind_traindata = ind_whole;  
            ind_traindata(ind_testdata) = [];  % the way to get complementary data set from an array
            
            % Now extract the training and testing data sets
            b_test = b(ind_testdata);
            t_test = t(ind_testdata);
            b_train = b(ind_traindata);
            t_train = t(ind_traindata);
            
            % Note the bumbers of columns of A_train and A_test are still the same
            [A_train, f_train_orig] = get_A(t_train, dconfig.type, dconfig.N, dconfig.param2);
            [A_test, f_test_orig] = get_A(t_test, dconfig.type, dconfig.N, dconfig.param2);
            A_train = A_train*sf;
            A_test = A_test*sf;
            
            % Now solve the model based on training data
            x_train = ncsolve(A_train, b_train, lambda_iter, L);  % note that this is the scaled intermittent solution
            
            % Compute the OCV error
            % doesn't need to scale back to the original axis since both A
            % and x are scaled, and the final product just make the scaling
            % factor cancel out each other. 
            ocv_error_iter(jj) = 1/n_testdata * norm(A_test*x_train - b_test, 2)^2;
        end
        
        % Now compute the mean and variance of OCV error for this lambda_iter
        ocv_error_mean(ii) = mean(ocv_error_iter);
        ocv_error_std(ii) = std(ocv_error_iter);
        
        disp(['OCV error mean = ', num2str(ocv_error_mean(ii)), ', OCV error standard deviation =', num2str(ocv_error_std(ii))])
        disp('OCV error for this iteration: ')
        disp(ocv_error_iter)
    end
    
    
    % Save the lambda path data
    data_out = [force_column_vector(lambda_arr), ocv_error_mean, ocv_error_std];  
    save('-ascii', 'ocv_lambda_path.dat', 'data_out')
    disp('Lambda path data saved to "ocv_lambda_path.dat".')
    
    % Now determine the optimal lambda based on those results
    % Results are also plotted in this function
    lambda_opt = find_lambda_opt_ocv(lambda_arr, ocv_error_mean, ocv_error_std);

    % Compute the final solution
    sol_scaled = ncsolve(A, b, lambda_opt, L);
    
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
    sol_scaled = ncsolve(A, b, lambda_opt, L);
    
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
        title(strcat('L-curve, \lambda=', num2str(lambda_opt)), 'fontsize', 12)
        print('-dpng', 'L_curve.png')   % save L-curve plot
        L_curve_data_out = [res_arr', smn_arr'];
        save('-ascii', 'L_curve.dat', 'L_curve_data_out')  % save L-curve data
    end
else 
    error('Regularization type unrecognized!')
end


% Now scale the solution back to its original axis
% !! Need to change this line for different models!!

sol = sol_scaled*sf;   % for capacitance DLTS model; this only works for linear discretization, or log discretization with sf ==1

if nargout >=1
    varargout=cell(nargout,1);
    varargout{1} = sol; 
    if nargout >=2
        varargout{2} = f_orig;
        if nargout >=3
            varargout{3} = A;
            if nargout >=4
                varargout{4} = lambda_opt;
            end
        end
    end
end

end


% function [A, f_scaled, f_orig] = get_AA(t, dconfig)
% % A wrapper based on get_A.m
% % Get discretization matrix A and f-space axis f
% 
% sf = dconfig.param1;   % scaling factor 
% if strcmpi(dconfig.type, 'linear')
%     [A, f_scaled] = get_A(t, dconfig.type, dconfig.N, dconfig.param2);
%     
%     f_orig = f_scaled * sf;  % original f-space before scaling
% elseif strcmpi(dconfig.type, 'log')
%     % Need to think about how the scaling affect the logarithmic discretization 
%     [A, f_orig] = get_A(t, dconfig.type, dconfig.N, dconfig.param2);
%     A = A*sf;
%     f_scaled = f_orig;   % note that scaling for log discretization is done w.r.t. the A matrix
% elseif strcmpi(dconfig.type, 'glq')
%     [A, f_scaled] = get_A(t, dconfig.type, dconfig.N);  % scale the t axis 
%     f_orig = f_scaled * sf;  % original f-space before scaling
% end
% % A = sparse(A);  % l1_ls program can take the advantage of sparse matrix
% 
% end


function v_out = force_column_vector(v_in)
% Convert the input to column vector if it is not
if ~iscolumn(v_in)
    v_out = v_in';
else
    v_out = v_in;
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