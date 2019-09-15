function Tx=i_meanexittime(xpoints,leftabsorb,rightabsorb,ux)
global g_grind;
if numel(ux)==1
    g_grind.fokkerplanck.tmp=ux;
else
    g_grind.fokkerplanck.tmp=[xpoints;ux];
end

g_grind.hfun.curr=i_getodehandle(1,'');
%generate initial values
solinit = bvpinit(xpoints,@mat4init);
%options= bvpset('RelTol', 1e-7,'AbsTol',1e-6);
% The BVP solver returns the structure 'sol'. 
%for exit time T(x)
try
   sol = bvp4c(@exit4ode,@exit4bc,solinit);
   Tx=deval(sol,xpoints);
   Tx=Tx(1,:);
catch err
   disp(err.message);
   Tx=NaN;
end

%Y = deval(sol,xpoints); %only if we want to interpolate between the mesh
%points

g_grind.fokkerplanck=rmfield(g_grind.fokkerplanck,'tmp');

  % -----------------------------------------------------------------------
  function res = exit4bc(ya,yb)
  % Boundary conditions. lambda is a required argument.
  %Left basin
  res = [ ya(2)   %reflecting
          yb(2)];  
  if leftabsorb
      res(1)=ya(1);%absorbing
  end

  if rightabsorb
      res(2)=yb(1);
  end

  end
  % -----------------------------------------------------------------------
      
   
% -------------------------------------------------------------------------
% Auxiliary function -- initial guess for the solution

function yinit = mat4init(x)
%   a=1/(xpoints(end)-xpoints(1));
%   if leftabsorb&&rightabsorb
%       %parabola
%       yinit=[(x-xpoints(end))*(x-xpoints(1))
%             2*x-xpoints(end)-xpoints(1)];
%   elseif leftabsorb
%       yinit=[-a*x+a
%           -a];
%   elseif rightabsorb
%       yinit=[a*x-a
%              a];
%   else
%       yinit=[0
%           0];
%   end
   yinit = [   cos(4*x)
             -4*sin(4*x) ];
end
end  % mat4bvp
% -------------------------------------------------------------------------
function g_dydx = exit4ode(g_x,g_y)
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

if numel(g_grind.fokkerplanck.tmp)==1
    g_ux=g_grind.fokkerplanck.tmp;
else
    g_ux=interp1(g_grind.fokkerplanck.tmp(1,:),g_grind.fokkerplanck.tmp(2,:),g_x);
end

g_Fx = feval(g_grind.hfun.curr, 1, g_x);
g_dydx = [              g_y(2)
       -(g_Fx*g_y(2)+g_ux)/(g_sigma^2/2)];
end
