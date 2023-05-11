clc
clear all
clc

format long;

a=-1/2;
b=-1/3;
m=0;

addpath("./858")
%     F = BESSELINT(A, NU, M)
%
%     X^M * J_NU(1)(A(1)*X) * ... * J_NU(k)(A(k)*X)

%I=besselint(1, 0, 2)
%I=besselint(1, 1/3, 0)
%I=besselint([5 1],[0 1],0)
%I=besselint([3 1], [-1/2 1/3], 1/6)

%     [f,err,eflag] = fri(a,nu,m,x0,tol)

tic
[res_858, error_858, eflag] = fri([1, 1], [a, b], m, 1, [eps, eps])
toc
disp("\n\n")


addpath("./935")
%    [sol,reterr,evals]=IIPBF(f,rhoI,tauI,a,b,abserr,relerr,typeI)
%
% Input
%          abserr = required minimum absolute error
%          relerr = required minimum relative error
%          a = non negative integer
%          b = non negative integer
%          rho = positive real number
%          tau = positive real number
%          type = 'JJ', 'JY', or 'YY' (refer to accompanying paper)
%
% Output
%          sol = estimated solution to integral
%          reterr = predicted error
%          evals = number of function evaluations

%fx = @(x) 1;
%[sol,reterr,evals]=IIPBF(fx,5,1,0,1,1e-10,1e-11,'JJ')

%     [result,err,neval]=dqagea(func,a,b,abserr,relerr)

%tic
%fx = @(x) x.^m;
%func = @(x) fx(x).*besselj(a, x).*besselj(b, x);
%[res_935, error_935, neval]=dqagea(func, 0, 1, eps, eps)
%toc