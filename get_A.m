function [A, varargout] = get_A(t, dtype, n, varargin)
% [A, varargout] = get_A(t, I, dtype, n, varargin)
% Generate the discretization matrix for PICTS
% Input args: 
%     t:  x-axis of experimental data (e.g., time), a column vector
%     dtype: discretization type, should be 'glq' (generalized Gauss-Laguerre quadrature), 'linear', or 'log'
%     N: discretization points in f space
%     varargin{1}: [lower_bound, upper_bound] for 'linear' or 'log'; 'glq' doesn't need this input 
%
% Output args: 
%     A: dicretized matrix
%     varargout{1}: f_out, sampled f-space
% 
% Written by Lijun Li
% Update log: 
%     Created: Oct, 2013

if length(varargin) >= 1
    dlimit = varargin{1};  % discretization limit for 'linear' or 'log'
end

A = zeros(length(t), n);

if strcmpi(dtype, 'glq')
    if length(varargin) >= 1
        warning('Ignore input argument starting from the 4th')
    end
    
    % Run generalized Gauss-Laguerre quadrature routine, with alpha=1
    % f: abscissa of Laguerra polynomial; a column vector
    % wt: weight of generalized Gauss-Laguerre quadrature; a column vector
    [f, wt] = gen_laguerre_rule2(n, 1, 0, 1);
    for ii = 1:length(t)
        for jj = 1:n
            A(ii, jj) = wt(jj)*exp(-f(jj)*(t(ii)-1));
        end
    end
    f_out = f;
    
elseif strcmpi(dtype, 'linear') || strcmpi(dtype, 'log')
    if length(dlimit) == 2
        lb = dlimit(1);
        ub = dlimit(2);
    else
        error('Linear or logarithmic discretization limit should be a 2-by-1 vector like [lower_bound, upper_bound]!')
    end
    if strcmpi(dtype, 'linear')
        f = linspace(lb, ub, n+1);
    else
        f = logspace(log10(lb), log10(ub), n+1);
    end
    df = zeros(n, 1);
    f2 = zeros(n,1);
    for ii = 1:n
        df(ii) = f(ii+1) - f(ii);
    end
    for ii = 1:n
        f2(ii) = (f(ii) + f(ii+1))/2;
    end
    for ii = 1:length(t)
        for jj = 1:n
            A(ii, jj) = f2(jj)*exp(-f2(jj)*t(ii))*df(jj);
        end
    end
    f_out = f2;
else
    error('Unrecognized discretization type!')
end

if nargout >=2
    varargout{1} = f_out;
end
end