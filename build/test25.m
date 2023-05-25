clc
clear all
clc

format long;

a = [sqrt(2), sqrt(3), sqrt(4)];
nu = [0, 0];
mu = linspace(-0.9+e/50,0,10);
m = 0;

addpath("./858")
% ALGORITHM 858: BESSELINT (Van Deun and Cools, 2006)
res_858 = zeros(length(mu), 1);
time_858 = zeros(length(mu), 1);
funeval_858 = zeros(length(mu), 1);

addpath("./935")
% ALGORITHM 935: IIBPF (Lucas et al.)
res_935 = zeros(length(mu), 1);
time_935 = zeros(length(mu), 1);
funeval_935 = zeros(length(mu), 1);

for k=1:length(mu)
  tic
  [f, res_858(k), error_858(k), funeval_858(k)] = fri(a, [mu(k), nu(1), nu(2)], m, 1, [eps, eps]);
  time_858(k) = toc;
  disp("\n\n")
  
  
  tic
  fx = @(x) besselj(mu(k), a(1)*x);
  func = @(x) fx(x).*besselj(nu(1), a(2)*x).*besselj(nu(2), a(3)*x);
  [res_935(k), error_935(k), funeval_935(k)]=dqagea(func, 0, 1, eps, eps);
  time_935(k) = toc;
endfor