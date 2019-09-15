function [lyap,matcontlyap] = firstlyap(t,x,odefile,Jac,Hess,varargin)
% i_firstlyap
%
% Computation of first Lyapunov coefficient at equilibrium. The sign of the first Lyapunov
% coefficient determines the sub- or supercritical nature of a Hopf bifurcation.
% (adapted from COCO example, added third derivative)
% res.lyapunov is the lyapunov coefficient.
% formula see http://www.scholarpedia.org/article/Andronov-Hopf_bifurcation

% compute eigenvector
if nargin==0
    %if no arguments are given, the current state is taken
    global g_grind;  %#ok<TLEV>
    t=0;
    x=i_initvar;
    parslist = sprintf(',%s',g_grind.pars{:});
    odefile = i_getodehandle('singlestep', parslist);
    Jac = i_getodehandle('Jacobian', parslist);
    Hess = i_getodehandle('Hessian', parslist);
    varargin=num2cell(evalin('base',sprintf('[%s];',parslist)));
end

A = Jac(t,x,varargin{:});
[X, D]   = eig(A);
[~, idx] = min(abs(real(diag(D))));
v  = X(:,idx);
om = imag(D(idx,idx));
vb = conj(v);
%if neutral saddle
if abs(om)<1e-6
    lyap=nan;
    matcontlyap=nan;
    return
end

% compute adjoint eigenvector
[X, D]   = eig(A');
[~, idx] = min(abs(real(diag(D)))); 
w = X(:,idx);

if om*imag(D(idx,idx))>0
    w = conj(w);
end
w = w/conj(w'*v);

% compute tensor products
B  = Hess(t,x,varargin{:});
B1 = zeros(numel(x),1);
B3 = zeros(numel(x),1);
for i=1:numel(x)
    Bmat  = reshape(B(i,:,:),[numel(x),numel(x)]);
    B1(i) = transpose(v)*Bmat*v;
    B3(i) = transpose(v)*Bmat*vb;
end
t1 = (2*sqrt(-1)*om*eye(numel(x))-A)\B1;
t2 = A\B3;
B2 = zeros(numel(x),1);
B4 = zeros(numel(x),1);
for i=1:numel(x)
    Bmat  = reshape(B(i,:,:),[numel(x),numel(x)]);
    B2(i) = transpose(vb)*Bmat*t1;
    B4(i) = transpose(v)*Bmat*t2;
end
%C(q,q,q) is missing: taken from MATCONT:
ten3Increment = (1e-5)^(3.0/5.0);
C = multilinear3(odefile,v,v,conj(v),x,ten3Increment,varargin{:}); %from MATCONT
matcontlyap = real(w'*(C+B2-2*B4))/2;  %this is the value MATCONT gives
lyap = matcontlyap/om; %real lyaphttp://www.scholarpedia.org/article/Andronov-Hopf_bifurcation
end

function vec3 = multilinear3(odefile,q1,q2,q3,x0,increment,varargin)
%--------------------------------------------------------------
% This file computes the multilinear function C(q1,q2,q3) where
% C = D^3(F(x0)), the 3rd derivative of the vectorfield wrt to phase
% variables only.
%--------------------------------------------------------------
if (q1==q2)
    if (q1==q3)
        vec3 = Cvvv(odefile,q1,x0,increment,varargin{:});
    else
        part1 = Cvvv(odefile,q1+q3,x0,increment,varargin{:});
        part2 = Cvvv(odefile,q1-q3,x0,increment,varargin{:});
        part3 = Cvvv(odefile,q3,x0,increment,varargin{:});
        vec3 = (part1 - part2 - 2.0*part3)/6.0;
    end
else
    part1 = Cvvv(odefile,q1+q2+q3,x0,increment,varargin{:});
    part2 = Cvvv(odefile,q1+q2-q3,x0,increment,varargin{:});
    part3 = Cvvv(odefile,q1-q2+q3,x0,increment,varargin{:});
    part4 = Cvvv(odefile,q1-q2-q3,x0,increment,varargin{:});
    vec3 = (part1 - part2 - part3 + part4)/24.0;
end

end

%----------------------------------------------------

function tempvec = Cvvv(odefile,vq,x0,increment,varargin)
vs=vq/norm(vq);
f1 = x0 + 3.0*increment*vs;
f2 = x0 +     increment*vs;
f3 = x0 -     increment*vs;
f4 = x0 - 3.0*increment*vs;
f1 = feval(odefile, 0, f1,varargin{:});
f2 = feval(odefile, 0, f2,varargin{:});
f3 = feval(odefile, 0, f3,varargin{:});
f4 = feval(odefile, 0, f4,varargin{:});
tempvec = (f1 - 3.0*f2 + 3.0*f3 - f4)/8.0*(norm(vq)/increment)^3;
end




