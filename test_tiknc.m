% Test Tikhonov regularization 

clear
close all

data_source = 1;   % 1 for Gaussian test data; 0 for real experimental data
delta_b = 0.5e-4;   % noise in b if using Gaussian test data, reference to the max value of input signal
% data_file = 'data_filename';

if data_source
    % Generate test data from Gaussian spectrum
    % The corresponding function is test_func.m
    N_t = 100;  
    t_end = 0.6;   % The length of the time-domain signal
    t = logspace(log10(1/N_t), log10(t_end), N_t)';
    mu = [8, 25, 80, 200];
    sigma =[1, 2, 4, 6];
    C = [1, 5, 7, 5];
    gaussian_ind = 1:4;
    
    b_bar =0;
    for ii = gaussian_ind
        b_bar = b_bar + C(ii)*t_func(t,mu(ii),sigma(ii));  % generate test data
    end
    b_bar = b_bar/b_bar(1);   % normalize the transient to its first point
    b_noise = randn(length(t), 1)*delta_b;
    b = b_bar + b_noise;
else
    data = importdata(data_file);
    t = data(:,1);
    b = data(:,2);
end

dstruct = struct();
rstruct = struct();

% Dicretization configuration
% dstruct.type = 'glq';
% dstruct.N = 120;
% dstruct.param1 = 1;
dstruct.type = 'log';
dstruct.N = 100;
dstruct.param1 = [0.1, 500];

% Regularization configuration
rstruct.type = 'L-curve1';
rstruct.order = 0;
rstruct.param1 = 0;
% rstruct.type = 'Morozov';
% rstruct.order = 2;
% rstruct.param1 = delta_b ;

[sol, f, A, lambda_opt] = tikregnc(t, b, dstruct, rstruct);


% Original spectrum (made up of Gaussians)
if data_source
    NT_original = 0;
    for ii = gaussian_ind
        NT_original = NT_original + 1/sqrt(2*pi)*C(ii)/sigma(ii)*exp(-(f-mu(ii)).^2/(2*sigma(ii).^2));
    end
end

fs1 = 16;
fs2 = 20;
figure(1)
if data_source
    semilogx(f, sol/max(sol), f, NT_original/max(NT_original))
else
    semilogx(f, data(:,2)/max(data(:,2)))
end
xlabel('f (ab. unit)', 'fontsize', fs1)
ylabel('N_T (ab. unit)', 'fontsize', fs1)
set(gca, 'fontsize', fs1, 'linewidth', 2)

figure(2)
semilogy(t,b/b(1), 'linewidth', 3)
xlabel('Time (ab. unit)', 'fontsize', fs2)
ylabel('I (ab. unit)', 'fontsize', fs2)
set(gca, 'fontsize', fs2, 'linewidth', 2)
axis([0,1,1e-6, 1])