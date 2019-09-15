function [tout,yout] = bv78(odefun,tspan,y0,options)
% Solve initial value problems y' = f(t,y).  The user interface is much
% like the solvers of the Matlab ODE Suite.  If VECTORIZED is FALSE, the
% program allows f(t,y) coded in a conventional way, but the default is
% TRUE.  The function must then be coded so that f(T,Y) returns YP
% with YP(:,k) = f(T(k),Y(:,k)) for each entry in T and corresponding
% column in Y.
% 
% reference
% Vectorized Solution of ODEs in Matlab with Control of Residual and Error
% Recent Advances in Computational and Applied Mathematics pp 259-271
% L.F. Shampine, 2011
% Mathematics Department
% Southern Methodist University
% Dallas, TX 75275
% shampine@smu.edu
% Abstract
% Vectorization is very important to the e±ciency of computation in the popular
% problem-solving environment Matlab. Here we develop an explicit Runge{
% Kutta (7,8) pair of formulas that exploits vectorization. Conventional Runge{
% Kutta pairs control local error at the end of a step. The new method controls
% the extended local error at 8 points equally spaced in the span of a step. This is
% a byproduct of a control of the residual at these points. A new solver based on
% this pair, odevr7, not only has a very much stronger control of error than the
% recommended Matlab solver ode45, but competes well at modest tolerances
% and is notably more e±cient at stringent tolerances.
% Keywords
% Matlab, vectorization, ordinary di®erential equations, initial value problems
if nargin<4||~isfield(options,'RelTol')||isempty(options.RelTol)
    re = 1e-3;
else
    re=options.RelTol;
end

if nargin<4||~isfield(options,'AbsTol')||isempty(options.AbsTol)
    ae = 1e-6;
else
    ae=options.AbsTol;
end

if nargin<4||~isfield(options,'Vectorized')||isempty(options.Vectorized)
    vectorized=false;
else
    vectorized=strcmp(options.Vectorized,'on');
end

haveOutputFcn=~isempty(options.OutputFcn);
outputFcn=options.OutputFcn;
if haveOutputFcn
  feval(outputFcn,tspan,y0,'init');
end
t = tspan(1);
tend = tspan(end);
dir = sign(tend - t);
tout_given = length(tspan) > 2;
y = y0(:);
yp = odefun(t,y);
if re < 100 * eps 
  re = 100 * eps;
  warning('MATLAB:REIncrease','RE has been increased to %g.',re)
end
if ae <= 0
  error('MATLAB:AENotPos', 'AE must be positive.');
end
   
% Compute an initial step size h using y'(t).
hmax = 0.1*abs(tend - t);
hmin = 16*max(eps(t),eps(tend));
rh = norm((yp ./ max(abs(y),(ae/re))),inf) / (0.8*re^(1/8));
absh = hmax;
if absh*rh > 1
    absh = 1/rh;
end
h = dir*max(absh,hmin);

neq = length(y);
if tout_given
    Lout = length(tspan);
    tout = tspan(:);
    yout = zeros(Lout,neq);
else
    Lout = 100;
    tout = zeros(Lout,1);
    yout = zeros(Lout,neq);
end

% Initialize output.
nout = 1;
tout(nout) = t;
yout(nout,:) = transpose(y);

done = false;
while ~done
    
    % Adjust step size near the end of the interval.
    if 1.1*abs(h) >= abs(tend - t)     % stretch
        h = tend - t;
        done = true;
    elseif 2*abs(h) >= abs(tend - t)   % look-ahead
        h = (tend - t)/2;
    end

    failed = false;
    while true  % Reduce h until the step to t + h succeeds.
        T = t + h*[1/6, 2/6, 3/6, 4/6, 5/6, 1]; 
           
       M = h*[ 2713/15120, 47/189, 27/112, 232/945, 725/3024, 9/35 ;    ...
          -15487/120960, 11/7560, 387/4480, 64/945, 2125/24192, 9/280 ; ...
           293/2835, 166/2835, 17/105, 752/2835, 125/567, 34/105;       ...
          -6737/120960, -269/7560, -243/4480, 29/945, 3875/24192, 9/280;...
           263/15120, 11/945, 9/560, 8/945, 235/3024, 9/35; ...
          -863/362880, -37/22680, -29/13440, -4/2835, -275/72576, 41/840 ];
             
        temp = y(:,ones(1,6)) + yp*h*[19087/362880, 1139/22680, ...
                                137/2688, 143/2835, 3715/72576, 41/840];
        D = yp(:,ones(1,6))*M;
        Y = temp + D;  % Locally second order approximation.

        for k = 1:7
            % For k = 1:6, the local order is k+2 at each point in 
            % the block and the global order is k+1. When k is increased 
            % from 6 to 7, the local order remains 8 at all the interior 
            % points, but the local order at the end of the step is 9, 
            % so the global order is 8.
            if vectorized
                YP = odefun(T,Y);
            else
                for col = 1:6
                    YP(:,col) = odefun(T(col),Y(:,col));
                end
            end
            D = YP*M;
            Y = temp + D;
            if k == 6
                % Remove "temp" from error estimate analytically
                % to reduce loss of significance.
                d7 = D(:,6);    
                y7 = Y(:,6);
            end
        end

        wt = ae/re + 0.5*(abs(Y(:,6)) + abs(y7));
        errnrm = norm((D(:,6) - d7)./wt,inf);
        absh = abs(h);
        if errnrm <= re  % successful step
            break
        else             % reduce h, try again
            if absh < 16*max(eps(t),eps(tend))
                error('grind:bv78:steptoosmall','Step size needed is too small for the precision.')
            end
            h = h*max(0.5,0.8*(re/errnrm)^(1/8));
            failed = true;
            done = false;
        end  
    end  % Loop until step succeeds.
    
    % If TOUT is given, see whether any output points lie in (t,t+h].
    % If so, evaluate the continuous extension to get solution there.
    if tout_given
        indices = find( (dir*(tout - t) > 0) & (dir*(t+h - tout) >= 0) );
        if ~isempty(indices)
            yint = continuous_extension(tout(indices));
            numout = length(indices);            
            yout(nout+1:nout+numout,:) = transpose(yint);
            nout = nout + numout;
        end
    else
        if nout + 6 > Lout       % Allocate more storage as needed.
            Lout = Lout + 100;
            tout(Lout) = 0;
            yout(Lout,1:neq) = 0;
        end
        tout(nout+1:nout+6) = transpose(T);
        yout(nout+1:nout+6,:) = transpose(Y);
        nout = nout + 6;
        
    end
    if haveOutputFcn
        if feval(outputFcn,T(end),transpose(Y(:,end)),'')
            done = true;
        end
    end
    if ~done
        t = T(end);
        y = Y(:,end);
        yp = odefun(t,y);
        if ~failed           % No increase after a failed attempt.
            absh = absh/max(0.2,1.25*(errnrm/re)^(1/8));
            h = sign(h)*min(absh,hmax);
        end
    end
      
end    

% Trim output arrays.
tout = tout(1:nout);
yout = yout(1:nout,:);


%===Nested function========================================================
function yce = continuous_extension(tint)
    s = ( transpose(tint(:)) - t )/h; 
    s2 = s .* s; s3 = s .* s2; s4 = s .* s3; s5 = s .* s4;
    s6 = s .* s5; s7 = s .*s6;
    yce = y(:,ones(size(tint))) ...
        + (324*h/5)*yp*((1/7)*s7-(7/12)*s6+(35/36)*s5-(245/288)*s4+...
          (203/486)*s3-(49/432)*s2+(5/324)*s) ...
        - (1944*h/5)*YP(:,1)*((1/7)*s7-(5/9)*s6+(31/36)*s5-...
          (145/216)*s4+(29/108)*s3-(5/108)*s2)...
        + (972*h)*YP(:,2)*((1/7)*s7-(19/36)*s6+(137/180)*s5-...
          (461/864)*s4+(13/72)*s3-(5/216)*s2)...
        - (1296*h)*YP(:,3)*((1/7)*s7-(1/2)*s6+(121/180)*s5-...
          (31/72)*s4+(127/972)*s3-(5/324)*s2)...
        + (972*h)*YP(:,4)*((1/7)*s7-(17/36)*s6+(107/180)*s5-...
          (307/864)*s4+(11/108)*s3-(5/432)*s2)...
        - (1944*h/5)*YP(:,5)*((1/7)*s7-(4/9)*s6+(19/36)*s5-...
          (65/216)*s4+(1/12)*s3-(1/108)*s2)...
        + (324*h/5)*YP(:,6)*((1/7)*s7-(5/12)*s6+(17/36)*s5-...
          (25/96)*s4+(137/1944)*s3-(5/648)*s2);
end % continuous_extension
%==========================================================================

end % BV78
