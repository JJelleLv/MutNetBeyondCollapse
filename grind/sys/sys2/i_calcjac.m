function g_Jacobian = i_calcjac(donumerical,~, g_l_N0, g_l_nlag)
global g_grind t; %t is needed
% g_l = g_grind.statevars.dim;
% g_Jacobian = zeros(g_l);
if ~donumerical&&~isempty(g_grind.syms.Jacobian)
%    g_h=evalin('base', sprintf('@(t,g_X1)reshape([%s],[%d,%d]);',sprintf('%s',g_grind.syms.Jacobian{:}),g_grind.statevars.dim,g_grind.statevars.dim));
%    g_Jacobian=g_h(t,g_l_N0);

   if ~isempty(g_grind.pars)
      eval(i_globalstr(g_grind.pars));
   end

   g_X1=g_l_N0; 
   if ~isempty(g_grind.externvars)
      for g_l_m = 1:length(g_grind.externvars)
         eval(sprintf('%s=externvar(%d, %s,t);',g_grind.externvars{g_l_m}.name,g_l_m, g_grind.externvars{g_l_m}.default));
      end

   end
   if isa(g_grind.syms.Jacobian,'function_handle')
       jac=i_getodehandle('Jacobian');
       g_Jacobian=jac(1,g_l_N0);
   else
       g_Jacobian=eval( sprintf('reshape([%s],[%d,%d]);',sprintf('%s;',g_grind.syms.Jacobian{:}),g_grind.statevars.dim,g_grind.statevars.dim));
   end
%    for g_l_i = 1:g_l
%       for g_l_j = 1:g_l
%          g_Jacobian(g_l_i, g_l_j) = eval(g_grind.syms.Jacobian{g_l_i, g_l_j});
%       end
%    end
   if g_grind.solver.backwards
      g_Jacobian = -g_Jacobian;
   end

else
   if ~isfield(g_grind.solver.opt,'numjac_fac')
       g_grind.solver.opt.numjac_fac=[];
   end

   if ~isfield(g_grind.solver.opt,'numjac_thresh')||isempty(g_grind.solver.opt.numjac_thresh)
       g_grind.solver.opt.numjac_thresh=1e-6+zeros(size(g_l_N0));
   end

   if g_grind.solver.haslags
       if nargin<4
           g_l_nlag=0;
       end

       lags_1 = repmat(g_l_N0, [1, length(g_grind.dde.lags)]);
       odehandle1=i_getodehandle(0,'');
       odehandle=@(t,g_l_lags0)lagodehandle(t,g_l_lags0,g_l_nlag,g_l_N0,lags_1,odehandle1);
   else
      odehandle=i_getodehandle(1,'');
   end

      %numjac(F,T,Y,FTY,THRESH,FAC,VECTORIZED) %numjac uses a rather precise algorithm 
      [g_Jacobian,g_grind.solver.opt.numjac_fac] = numjac(odehandle,t,g_l_N0,odehandle(t,g_l_N0),g_grind.solver.opt.numjac_thresh,g_grind.solver.opt.numjac_fac,false);
%       N0=repmat(g_l_N0',g_l+1,1);
%       N0(2:end,:)=N0(2:end,:)+eye(g_l)*delta;
%       Nres=i_runsinglestep(1,N0,true);
%       g_Jacobian=((Nres(2:end,:)-repmat(Nres(1,:),g_l,1))./delta)';
end
function res=lagodehandle(t,X1,nlag,N0,lags,odehandle)
if nlag>0
   lags(:,nlag)=X1;
   res=odehandle(t,N0,lags);
else
   res=odehandle(t,X1,lags);
end    

