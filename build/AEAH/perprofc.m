function perprofc

% To add new examples:
%  exX_X = {a, nu, m, c, tol, flag, exact}
% where 'flag' is 0 for relative tolerance or 1 for absolute
%  and 'exact' is a vector containing the exact result for each
%  entry in 'c'
% Then add 'exX_X' to the list in variable 'exples' below

% ----- Example 1.1 (relative tolerance)-------------------------------

c = linspace(0.001,25,500);
tol = 500*eps;
ex1_1 = {1, 0, 0, c, tol, 0, 1./sqrt(c.^2+1)};

% ----- Example 1.2 (absolute tolerance)-------------------------------

ex1_2 = {1, 0, 0, c, tol, 1, 1./sqrt(c.^2+1)};

% ----- Example 2.1 (relative tolerance) ------------------------------

nu = -1/3;
ex2_1 = {1, nu, nu, c, tol, 0, 2^nu*gamma(nu+1/2)/sqrt(pi)...
    ./(1+c.^2).^(nu+1/2)};

% ----- Example 2.2 (absolute tolerance) ------------------------------

ex2_2 = {1, nu, nu, c, tol, 1, 2^nu*gamma(nu+1/2)/sqrt(pi)...
    ./(1+c.^2).^(nu+1/2)};

% ----- Example 3.1 (relative tolerance) ------------------------------

nu = 1/3;
ex3_1 = {1, nu, -1, c, tol, 0, (sqrt(1+c.^2)-c).^nu/nu};

% ----- Example 3.2 (absolute tolerance) ------------------------------

ex3_2 = {1, nu, -1, c, tol, 1, (sqrt(1+c.^2)-c).^nu/nu};

% ----- Example 4.1 (relative tolerance) ------------------------------

c = linspace(0.001/2,25/2,500);
[k,e] = ellipke(1./(1+c'.^2)); k = k'; e = e';
ex4_1 = {[1,1], [0,1], 1, 2*c, tol, 0, (k-e)./(2*pi*sqrt(1+c.^2))};

% ----- Example 4.2 (absolute tolerance) ------------------------------

ex4_2 = {[1,1], [0,1], 1, 2*c, tol, 1, (k-e)./(2*pi*sqrt(1+c.^2))};

% ----- Example 5.1 (relative tolerance) ------------------------------

ex5_1 = {[1,1], 0, 0, 2*c, tol, 0, k/pi./sqrt(1+c.^2)};

% ----- Example 5.2 (absolute tolerance) ------------------------------

ex5_2 = {[1,1], 0, 0, 2*c, tol, 1, k/pi./sqrt(1+c.^2)};

% ---------------------------------------------------------------------

exples = {ex1_1, ex1_2, ex2_1, ex2_2, ex3_1, ex3_2, ex4_1, ex4_2, ...
    ex5_1, ex5_2};
disp(sprintf('Nr\tReq\tUEst\tMnDev\tMaxDev\tFalse\tLoss\tMax'));
for j = 1:length(exples),
    par = exples{j};
    c = par{4}; tol = par{5};
    [t1,m1a,m1b,t2,t3,m3,e,ex] = testperf(par{1},par{2},par{3},c,tol,...
        par{6},par{7});
    if par{6}, s = 'A'; else s = 'R'; end
    disp(sprintf('%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d',ceil(j/2),s,t1,m1a,...
        m1b,t2,t3,m3));
    figure; semilogy(c,e,'b',c,ex,'r',c,tol*ones(size(c)),'g');
end

% ---------------------------------------------------------------------

function [t1,m1a,m1b,t2,t3,m3,e,ex] = testperf(a,nu,m,c,tol,flag,fx)

n = length(c);
t1 = 0; t2 = 0; t3 = 0;
m1a = []; m1b = 0;
m3 = 0;

for i = 1:n,
   if flag,
       [f,err] = besselintc(a,nu,m,c(i),[],tol);
       ex(i) = abs(f-fx(i));
       e(i) = err(2);
   else,
       [f,err] = besselintc(a,nu,m,c(i),tol);
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
