function [f,err,eflag,neval] = fri(a,nu,m,x0,tol)
%FRI Finite range integration
neval = 0;

% get rough approximation to number of zeros in integrand
n = max(ceil(sum(x0*a/pi+1/4-nu/2)),0)+1;
h = x0/n; % length of subintervals
z = [0:h:x0];

f = 0; err1 = 0;
% treat possible algebraic singularity in 0 (first interval)
p = sum(nu) + m; % degree of singularity
if abs(fix(p)-p) > abs(p)*50*eps, % fractional degree -> extrapolate
    xs = z(2); z = z(2:end);
    rho = 1/4; f1 = 0; i = 0;
    [f2,err,eflag,funceval1] = nqf(@fun,[xs/2,xs],tol,a,nu,m);
    neval = neval + funceval1;
    I1(1) = f2; xs = xs/2;
    while ~eflag & abs(f2-f1) > tol(1)*abs(f2) & abs(f2-f1) > tol(2),
        i = i + 1;
        [I2(1),err,eflag,funceval2] = nqf(@fun,[xs*rho,xs],tol,a,nu,m);
        neval = neval + funceval2;
        I2(1) = I1(1) + I2(1);
        for j = 1:i, % extrapolate in sense of Richardson
            I2(j+1) = (I2(j) - rho^(p+j)*I1(j)) / (1 - rho^(p+j));
        end;
        f1 = I1(i); f2 = I2(i+1); I1 = I2; xs = xs*rho;
    end
    f = f2;
    err1 = abs(f-f1);
    if length(z) == 1,
        err = err1;
        return
    end
end

z = [z(1:2:end-1) x0]; % take 2 intervals together to save work

% split integrand at equidistant points to integrate
V = [z(1:end-1)' z(2:end)'];
[f2,err,eflag,funceval3] = nqf(@fun,V,tol,a,nu,m);
neval = neval + funceval3;
f = f + f2;
err = err + err1;

% ----------------------------------------------------------------------

function f = fun(x,a,nu,m)
%FUN Integrand

f = x.^m;
for i = 1:length(a),
    f = f.*besselj(nu(i),a(i)*x);
end

% ----------------------------------------------------------------------

function [I,err,eflag,funceval] = nqf(fun,X,tol,varargin)
%NQF Numerical quadrature formula
funceval = 0;
% some hard coded Gauss-Legendre rules
xw5=[0,                         0.28444444444444444444;
     0.53846931010568309104,    0.23931433524968323402;
     0.90617984593866399280,    0.11846344252809454375;];

xw7=[0,                         0.20897959183673469388;
     0.40584515137739716691,    0.19091502525255947247;
     0.74153118559939443986,    0.13985269574463833395;
     0.94910791234275852453,    0.064742483084434846630;];

xw11=[0,                        0.13646254338895031536;
      0.26954315595234497233,   0.13140227225512333109;
      0.51909612920681181593,   0.11659688229599523996;
      0.73015200557404932409,   0.093145105463867125710;
      0.88706259976809529908,   0.062790184732452312314;
      0.97822865814605699280,   0.027834283558086833246;];

xw15=[0,                        0.10128912096278063644;
      0.20119409399743452230,   0.099215742663555788230;
      0.39415134707756336990,   0.093080500007781105516;
      0.57097217260853884754,   0.083134602908496966776;
      0.72441773136017004742,   0.069785338963077157225;
      0.84820658341042721620,   0.053579610233585967505;
      0.93727339240070590431,   0.035183023744054062354;
      0.98799251802048542849,   0.015376620998058634177;];

xw19=[0,                        0.080527224924391847988;
      0.16035864564022537587,   0.079484421696977173823;
      0.31656409996362983199,   0.076383021032929833389;
      0.46457074137596094572,   0.071303351086803305889;
      0.60054530466168102347,   0.064376981269668113838;
      0.72096617733522937862,   0.055783322773666997358;
      0.82271465653714282498,   0.045745010811224999731;
      0.90315590361481790164,   0.034522271368820613291;
      0.96020815213483003085,   0.022407113382849800168;
      0.99240684384358440319,   0.0097308941148632385170;];

xw23=[0,                        0.066827286093053087676;
      0.13325682429846611093,   0.066231019702348308687;
      0.26413568097034493053,   0.064452861094041074987;
      0.39030103803029083142,   0.061524542153364765234;
      0.50950147784600754969,   0.057498320111205682472;
      0.61960987576364615639,   0.052446045732270705037;
      0.71866136313195019446,   0.046457883030017573738;
      0.80488840161883989215,   0.039640705888359477462;
      0.87675235827044166738,   0.032116210704262926063;
      0.93297108682601610235,   0.024018835865542334286;
      0.97254247121811523196,   0.015494002928489722153;
      0.99476933499755212352,   0.0067059297435708860456;];

xw={xw5,xw7,xw11,xw15,xw19,xw23};

if (tol(1) <= 50*eps) & (tol(2) <= 50*eps),
    k = 4; % avoid taking not enough points
else
    k = 1;
end
N = size(X,1);
eflag = 0;
I1 = 0; % dummy value
for j = k:length(xw),
    If = 0;
    w = [xw{j}(end:-1:2,2); xw{j}(:,2)];
    for i = 1:N, % integrate over different subintervals
        L = (X(i,2)-X(i,1))/2; M = (X(i,2)+X(i,1))/2;
        x = L*[-xw{j}(end:-1:2,1);xw{j}(:,1)]' + M; % maps the nodes in the actual interval
        If = If + 2*L*feval(fun,x,varargin{:})*w;
        funceval = funceval + length(w);
    end
    I2 = I1;
    I1 = If;
    err = abs(I2-I1);
    if (j > k) & ((err <= tol(1)*abs(I1)) | (err <= tol(2))),
        I = I1; % precision reached
        return
    end
end
% requested precision not reached
I = I1; eflag = 1;