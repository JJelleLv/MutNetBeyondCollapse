function ux=i_exitprob(xpoints,toleft)
%exit probability if there are two absorbing boundaries (it is not usefull
%to do this analysis if one or two of the boundaries is reflecting as there is
%then only one exit so the prob will be 0 for the reflecting and 1 for the 
%absorbing boundary 

%generate initial values
global g_grind;
solinit = bvpinit(xpoints,@mat4init);
g_grind.hfun.curr=i_getodehandle(0,'');
% The BVP solver returns the structure 'sol'. 
   %for probability u(x)
try
%   options= bvpset('RelTol', 1e-7,'AbsTol',1e-6);
   sol = bvp4c(@prob4ode,@prob4bc,solinit);
   ux = deval(sol,xpoints); %only if we want to interpolate between the mesh
   ux=ux(1,:);
   %points
catch err
   disp(err.message);
   ux=NaN;
end



  function res = prob4bc(ya,yb)
  % Boundary conditions. lambda is a required argument.
       
  %Probability of exit through the left boundary if we are in the middle basin  
  if toleft
  res = [ ya(1)-1
          yb(1)];
  else
  %Probability of exit through the right boundary if we are in the middle basin     
    res = [ ya(1)
          yb(1)-1];    
  end
  % -----------------------------------------------------------------------
  end
  % -----------------------------------------------------------------------
% -------------------------------------------------------------------------
% Auxiliary function -- initial guess for the solution

function yinit = mat4init(x)
%   a=1/(xpoints(end)-xpoints(1));
%   if ~toleft
%       yinit=[-a*x+a
%           -a];
%   else
%       yinit=[a*x-a
%              a];
%   end
% 
   yinit = [   cos(4*x)
             -4*sin(4*x) ];
end

end  % mat4bvp  
    

% -------------------------------------------------------------------------


function g_dydx = prob4ode(g_x,g_y)
global g_grind
if ischar(g_grind.fokkerplanck.sigma)
       for g_i=1:length(g_grind.fokkerplanck.sigmavars)
           g_ix=i_getno(g_grind.fokkerplanck.sigmavars{g_i});
           if g_ix.isvar
               eval(sprintf('%s=g_x(%d);',g_grind.fokkerplanck.sigmavars{g_i},g_ix.no));
           else
               eval(sprintf('global %s;',g_grind.fokkerplanck.sigmavars{g_i}));
           end

       end

       g_sigma = eval(g_grind.fokkerplanck.sigma);
else
   g_sigma = g_grind.fokkerplanck.sigma;
end

g_Fx = feval(g_grind.hfun.curr, 1, g_x);
 g_dydx = [              g_y(2)
          -g_Fx*g_y(2)/(g_sigma^2/2) ];
end
  % -----------------------------------------------------------------------

  % -----------------------------------------------------------------------
      
 
