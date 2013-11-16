function yt = t_func(t, mu, sigma)
% Calculate the time-domain data from a Gaussian spectrum in f-space
% t: input time axis, a vector
% mu, sigma: Gaussian mean and standard deviation, must be scalars

a = (mu-sigma^2*t)/(sqrt(2)*sigma);

A = exp(-mu*t + (sigma*t).^2/2);
B = (exp(-a.^2) + sqrt(pi)*a.*(1+erf(a)));
yt = sigma/sqrt(2*pi) * A.*B;
end