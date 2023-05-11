function perprofr

% To add new examples:
%  exX_X = {a, nu, m, c, tol, flag, exact}
% where 'flag' is 0 for relative tolerance or 1 for absolute
%  and 'exact' is a vector containing the exact result for each
%  entry in 'c'
% Then add 'exX_X' to the list in variable 'exples' below

% ----- Example 1.1 (relative tolerance)-------------------------------

rr = linspace(1,4,500);
nu = 4/9;
tol = 500*eps;
ex1_1 = {1, nu, nu+1, rr, tol, 0, rr.^nu.*besselk(nu,rr)};

% ----- Example 1.2 (absolute tolerance)-------------------------------

ra = linspace(1,25,500);
ex1_2 = {1, nu, nu+1, ra, tol, 1, ra.^nu.*besselk(nu,ra)};

% ----- Example 2.1 (relative tolerance)-------------------------------

nu = -1/2;
ex2_1 = {1, nu, nu+1, rr, tol, 0, rr.^nu.*besselk(nu,rr)};

% ----- Example 2.2 (absolute tolerance)-------------------------------

ex2_2 = {1, nu, nu+1, ra, tol, 1, ra.^nu.*besselk(nu,ra)};

% ----- Example 3.1 (relative tolerance)-------------------------------

b = 2; nu = 0;
ex3_1 = {[1,b], nu, 1, rr, tol, 0, ...
    real(pi*i/2*besselj(nu,i*rr).*besselh(nu,1,i*rr*b))};

% ----- Example 3.2 (absolute tolerance)-------------------------------

ex3_2 = {[1,b], nu, 1, ra, tol, 1, ...
    real(pi*i/2*besselj(nu,i*ra).*besselh(nu,1,i*ra*b))};

% ----- Example 4.1 (relative tolerance)-------------------------------

nu = -8/9;
ex4_1 = {[1,b], nu, 1, rr, tol, 0, ...
    real(pi*i/2*besselj(nu,i*rr).*besselh(nu,1,i*rr*b))};

% ----- Example 4.2 (absolute tolerance)-------------------------------

ex4_2 = {[1,b], nu, 1, ra, tol, 1, ...
    real(pi*i/2*besselj(nu,i*ra).*besselh(nu,1,i*ra*b))};

% ----- Example 5.1 (relative tolerance)-------------------------------

b = [1,1.1,1.2,1.3]; mu = [0,0,0,0];
a = 5; nu = 1; rho = 5;
fx = -besselk(nu,rr*a).*rr.^(rho-2);
for j = b,
    fx = fx.*besseli(0,j*rr);
end
ex5_1 = {[b,a], [mu,nu], rho-1, rr, tol, 0, fx};

% ----- Example 5.2 (absolute tolerance)-------------------------------

fx = -besselk(nu,ra*a).*ra.^(rho-2);
for j = b,
    fx = fx.*besseli(0,j*ra);
end
ex5_2 = {[b,a], [mu,nu], rho-1, ra, tol, 1, fx};

% ---------------------------------------------------------------------

exples = {ex1_1, ex1_2, ex2_1, ex2_2, ex3_1, ex3_2, ...
    ex4_1, ex4_2, ex5_1, ex5_2};
disp(sprintf('Nr\tReq\tUEst\tMnDev\tMaxDev\tFalse\tLoss\tMax'));
for j = 1:length(exples),
    par = exples{j};
    r = par{4}; tol = par{5};
    [t1,m1a,m1b,t2,t3,m3,e,ex] = testperf(par{1},par{2},par{3},r,tol,...
        par{6},par{7});
    if par{6}, s = 'A'; else s = 'R'; end
    disp(sprintf('%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d',ceil(j/2),s,t1,m1a,...
        m1b,t2,t3,m3));
    figure; semilogy(r,e,'b',r,ex,'r',r,tol*ones(size(r)),'g');
end

% ---------------------------------------------------------------------

function [t1,m1a,m1b,t2,t3,m3,e,ex] = testperf(a,nu,m,r,tol,flag,fx)

n = length(r);
t1 = 0; t2 = 0; t3 = 0;
m1a = []; m1b = 0;
m3 = 0;

for i = 1:n,
   if flag,
       [f,err] = besselintr(a,nu,m,r(i),[],tol);
       ex(i) = abs(f-fx(i));
       e(i) = err(2);
   else,
       [f,err] = besselintr(a,nu,m,r(i),tol);
       ex(i) = max(eps,abs(1-f/fx(i)));
       e(i) = err(1);
   end
   if (ex(i) > e(i)),
      t1 = t1 + 1;
      m1a = [m1a,ceil(log10(ex(i))-log10(e(i)))];
      m1b = max(ceil(log10(ex(i))-log10(e(i))),m1b);
   end 
   if (e(i) > tol) & (ex(i) <= tol),
      t2 = t2 + 1;
   end
   if ex(i) > tol,
      t3 = t3 + 1;
      m3 = max(ceil(log10(ex(i))-log10(tol)),m3);
   end
end

t1 = round(t1/n*100); t2 = round(t2/n*100); t3 = round(t3/n*100);
if isempty(m1a), m1a = 0; else m1a = round(mean(m1a)); end
