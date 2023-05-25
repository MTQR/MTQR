clc
clear all
clc

format long;

a = [0.1, 10];
nu = [1/2, -1/2];
m = nu(2)-nu(1)+1;
u = linspace(0.001,0.1,20);

addpath("./AEAH")
% ALGORITHM AEAH-CPC: BESSELINTR (Van Deun and Cools, 2008)
res_AEAH = zeros(length(u), 1);
error_AEAH = zeros(length(u), 1);
time_AEAH = zeros(length(u), 1);
funeval_AEAH = zeros(length(u), 1);

addpath("./935")
% ALGORITHM 935: IIBPF (Lucas et al.)
res_935 = zeros(length(u), 1);
error_935 = zeros(length(u), 1);
time_935 = zeros(length(u), 1);
funeval_935 = zeros(length(u), 1);

for j=1:length(u)
  
  tic
  [res_AEAH(j), error_AEAH(j), funeval_AEAH(j)] = frir(a, nu, m, u(j), 1, [eps, eps]);
  time_AEAH(j) = toc;
  
  tic
  fx = @(x) x.^m ./ (x.^2 + u(j).^2);
  func = @(x) fx(x).*besselj(nu(1), a(1)*x).*besselj(nu(2), a(2)*x);
  [res_935(j), error_935(j), funeval_935(j)]=dqagea(func, 0, 1, eps, eps);
  time_935(j) = toc;
endfor


