%sets a parameter and sets NO if iX is a state variable
function N1 = i_setparvar(iX, N0, aval)
global g_grind;
N1=N0;
if ~isempty(aval)
   if ~isempty(iX.transform)
       if nargin(iX.transform.invfun)==2
          aval=iX.transform.invfun(aval,N0(iX.no));
       else
          aval=iX.transform.invfun(aval);
       end
   end
   if iX.ispar
      if isempty(iX.ndx)
         assignin('base', g_grind.pars{iX.no}, aval);
      else
         multassignin('base',sprintf('%s(%d)',g_grind.pars{iX.no},iX.ndx),aval);
      end

   elseif iX.isext
      if isempty(iX.ndx)
         assignin('base', g_grind.externvars{iX.no}.name, aval);
      else
         multassignin('base',sprintf('%s(%d)',g_grind.externvars{iX.no}.name,iX.ndx),aval);
      end

   else
      N1(iX.no) = aval;
   end

end

