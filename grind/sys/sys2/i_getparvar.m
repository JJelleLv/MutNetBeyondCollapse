%gets a parameter and gets [] if iX is a state variable
function aval=i_getparvar(iX)
global g_grind;
if iX.ispar
   if isempty(iX.ndx)
      aval = evalin('base', g_grind.pars{iX.no});
   else
      aval =evalin('base',sprintf('%s(%d)',g_grind.pars{iX.no},iX.ndx));
   end

elseif iX.isext
   if isempty(iX.ndx)
      aval = evalin('base', g_grind.externvars{iX.no}.name);
   else
      aval =evalin('base',sprintf('%s(%d)', g_grind.externvars{iX.no}.name,iX.ndx));
   end

else
   aval=[];
end

