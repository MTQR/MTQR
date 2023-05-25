clc
clear all
clc

format long;

a = [1, 3/2;
     1, 3;
     1, 1;
     1, 1];

nu = [0, 1;
      -1/3, 0;
      -1/2, -1/3; 
      0, -pi/4];
  
m = [-1/2, 1/6, 0, 0];

addpath("./858")
% ALGORITHM 858-TOMS: BESSELINT (Van Deun and Cools, 2006)
res_858 = zeros(length(m), 1);
time_858 = zeros(length(m), 1);
funeval_858 = zeros(length(m), 1);

addpath("./935")
% ALGORITHM 935: IIBPF (Lucas et al.)
res_935 = zeros(length(m), 1);
time_935 = zeros(length(m), 1);
funeval_935 = zeros(length(m), 1);

for k=1:length(m)
  tic
  [f, res_858(k), error_858(k), funeval_858(k)] = fri(a(k,:), nu(k,:), m(k), 1, [eps, eps]);
  time_858(k) = toc;
  disp("\n\n")
  
  
  tic
  fx = @(x) x.^m(k);
  func = @(x) fx(x).*besselj(nu(k,1), a(k,1)*x).*besselj(nu(k,2), a(k,2)*x);
  [res_935(k), error_935(k), funeval_935(k)]=dqagea(func, 0, 1, eps, eps);
  time_935(k) = toc;
endfor