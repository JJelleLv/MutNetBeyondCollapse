%DWIENER   Add stochasticity to a differential equation (Wiener process)
%  Consider the stochastic differential equation:
%  dX = f(X) dt + g(X) dW
%  in which
%     X is the state variable
%     f(X) the "drift term" as function of X
%     g(X) the "diffusive term" as function of X (or a constant in case of additive noise)
%     dW the increment of the Wiener process or Brownian motion
% 
%  This model can be defined in GRIND as:
%  X' = f(X) + dwiener(g(X),g'(X))
%  in which g'(X) is the derivative of g(X) to X. (thus zero for additive noise)
%  If this function is added to a differential equation, the solver is set to Euler integration
%  and the model is solved with the explicit Milstein method (if g' is defined) or with the Euler-Maruyama
%  scheme for additive noise. (both are based on the Ito integral). <a href="matlab:help solver">solver -n</a> can be used 
%  to get positive solutions only.
%  The function can return a vector if the sigmas is a vector. The elements of the vector are independent. 
%  Optionally they can be correlated approximately by a certain correlation if the 'corr' argument is used.
%  If the parameter 'alpha'(=stability) is between 0 and 2, a heavy tailed "alpha-stable" distribution is used instead of 
%  the Normal distribution. (by default the other parameters are: beta(=skewness)=0, gamma(=scaling)=1, delta(=location)). 
%  The result is then a Lévy flight (for a Poisson process see <a href="matlab:help djump">djump</a>). For the method see:
%  Siegert, S. and R. Friedrich. 2001. Modelling of nonlinear Lévy processes by data analysis. Physical Review E - 
%  Statistical, Nonlinear, and Soft Matter Physics 64:411071-4110712.
%
%
%
%  Usage:
%  DWIENER(SIGMA) - adding additive noise with Normal distribution with standard deviation SIGMA.
%  DWIENER(GX,GXACCENT) - adding multiplicative noise solved with Milstein method.
%  DWIENER(GX) - if GXACCENT is ignored the less efficient Euler-Maruyama method is used.
%  DWIENER('argname',argvalue,...) - Valid argument name-value pairs [with type]:
%     'alpha' [number>0 and number<=2] - stability parameter [0 2] of the alpha-stable distribution (default=2=Normal)
%     'beta' [number>=-1 and number<=1] - skewness parameter [-1 1] of the alpha-stable distribution, negative is 
%    skewed to negative values (default=0=symmetric) Note: if alpha=2, beta has no effect
%     'corr' [string] - enter a required correlation between the random variables (or a correlation matrix) (can be a function or parameter)
%     'dfundx' [string] - the derivative of the GX function
%     'fun' [string] - the function or value of the standard deviation SIGMA.
%     'number' [integer] - number of the dwiener function
%  DWIENER('-opt1','-opt2',...) - Valid command line options:
%     '-i' - option that is used internally before each run, to update the g_grind.solver.dwiener field
%
%  See also modelpanel, model, rednoise, djump, euler, stochast_heun csolvers
%
%   Reference page in Help browser:
%      <a href="matlab:commands('dwiener')">commands dwiener</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function res=dwiener(at,gY,dgY_dY,No,varargin)
global g_grind t;
if nargin>4||ischar(at)
    if nargin==1
        if isfield(g_grind,'solver')&&isfield(g_grind.solver,'dwiener')
            for i=1:length(g_grind.solver.dwiener.args)
                g_grind.solver.dwiener.args(i).L=[];
            end
        end
        if strcmp(at,'-i')
            return;
        end
        args={at};
    elseif nargin==2
        args={at,gY};
    elseif nargin==3
        args={at,gY,dgY_dY};
    elseif nargin>=4
        args=[{at,gY,dgY_dY,No} varargin];
    end
    fieldnams={'fun', 's', 'the function or value of the standard deviation SIGMA.','sigma';...
        'dfundx', 's', 'the derivative of the GX function','';...
        'corr', 's', 'enter a required correlation between the random variables (or a correlation matrix) (can be a function or parameter)',0;...
        'number', 'i', 'number of the dwiener function',1;...
        'alpha', 'n>0&n<=2', 'stability parameter [0 2] of the alpha-stable distribution (default=2=Normal)',2;...
        'beta', 'n>=-1&n<=1', 'skewness parameter [-1 1] of the alpha-stable distribution, negative is',1}';
    res1=i_parseargs(fieldnams,'fun,dfundx','-i',args);
    %     if isfield(res1,'corr')
    %         g_grind.solver.dwiener.args.corr=res1.corr;
    %     else
    %         g_grind.solver.dwiener.args.corr=[];
    if ~isfield(res1,'alpha')
        res1.alpha=2;
    end
    if ~isfield(res1,'beta')
        res1.beta=0;
    end
    %     end
    if nargout>0
        res=res1;
    end
    return;
end
if nargin==3
    No=dgY_dY;
    dgY_dY=0;
end
if ~isempty(g_grind.solver.dwiener.args(No).corr)&&(isempty(g_grind.solver.dwiener.args(No).L)||size(g_grind.solver.dwiener.args(No).L,1)~=size(gY,1))
    %update the dwiener.args.L field that is used to make the noise sources
    %correlated at a certain level. Note that we get a problem with the
    %chol function if the covariance matrix is not positive definite
    %currently we cannot use several different dwiener sets.
    if numel(gY)==1&&~g_grind.statevars.vector
        warning('grind:dwiener:corr','Correlation has been defined in scalar DWIENER, but it only works for vectors');
    end
    if numel(gY)==1&&g_grind.statevars.vector
        siz=g_grind.statevars.dims{g_grind.solver.dwiener.args(No).var}.dim1*g_grind.statevars.dims{g_grind.solver.dwiener.args(No).var}.dim2;
    else
        siz=numel(gY);
    end
    corr=evalin('base',g_grind.solver.dwiener.args(No).corr);
    if numel(corr)==1
        if corr==1
            g_grind.solver.dwiener.args(No).L=zeros(siz);
            g_grind.solver.dwiener.args(No).L(1,:)=1;
        else
            g_grind.solver.dwiener.args(No).L=chol(setdiagon(zeros(siz)+corr,1));
        end
    else
        g_grind.solver.dwiener.args(No).L=chol(corr);
    end
end
alpha=g_grind.solver.dwiener.args(No).alpha;
% if at==1&&~any(strcmpi(solver('name'),{'euler','c.euler'}))
%    if g_grind.solver.isdiffer
%       error('GRIND:dwiener:diffEquation','dwiener cannot be used for a difference equation, use randn() instead');
%    end
%    error('GRIND:dwiener:NoEuler','dwiener process needs Euler integration');
% end
if size(at,1)>1
    %remake a similar (BUT NOT THE SAME!) data set if used in a function.
    gY=transpose(gY);
    dgY_dY=transpose(dgY_dY);
    h=1;
    if alpha==2
        dW=randn(length(at),length(gY));
    else
        dW=stable_rnd(alpha,g_grind.solver.dwiener.args(No).beta,1,0,[length(at),length(gY)]);
    end
    if ~isempty(g_grind.solver.dwiener.args(No).L)
        dW=dW*g_grind.solver.dwiener.args(No).L;
    end
    warning('GRIND:dwiener:newseries','Note that dwiener generated a new (but similar) noise series');
elseif size(at,2)>1
    %remake a similar (BUT NOT THE SAME!) data set if used in a function.
    gY=transpose(gY);
    dgY_dY=transpose(dgY_dY);
    h=1;
    if alpha==2
        dW=randn(length(gY),length(at));
    else
        dW=stable_rnd(alpha,g_grind.solver.dwiener.args(No).beta,1,0,[length(at),length(at)]);
    end
    if ~isempty(g_grind.solver.dwiener.args(No).L)
        dW=transpose(g_grind.solver.dwiener.args(No).L)*dW;
    end
else
    h=g_grind.solver.opt.StepSize;
    if at== t + 0.00001 %null uses this
        dW=zeros(size(gY));
    else
        if alpha==2
            dW=randn(size(gY)).*sqrt(h);
        else
            dW=stable_rnd(alpha,g_grind.solver.dwiener.args(No).beta,1,0,size(gY)).*h./h.^(1/alpha);
        end
        if ~isempty(g_grind.solver.dwiener.args(No).L)
            if numel(dW)==1&&g_grind.statevars.vector
                siz=[g_grind.statevars.dims{g_grind.solver.dwiener.args(No).var}.dim1 g_grind.statevars.dims{g_grind.solver.dwiener.args(No).var}.dim2];
                dW=dW+zeros(siz);
            else
                siz=size(gY);
            end
            dW=reshape(transpose(g_grind.solver.dwiener.args(No).L)*dW(:),siz);
        end
    end
end

%Milstein scheme (reduces to Euler-Maruyama if d(gY)/dY=0)
res=(gY.*dW+gY.*dgY_dY.*(dW.^2-h))./h;% divide by h as the Euler routine multiplies with h
