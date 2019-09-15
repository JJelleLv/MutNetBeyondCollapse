function resetpars(silent)
global g_grind;
if ~isempty(g_grind.commands)
 %  g_l_linefeed = sprintf('\n');
   g_l_s = '';
   g_l_ok = 1;
   %first try to initialize all parameters (to avoid errors if other commands are given)
   for g_l_i = 1:size(g_grind.commands, 2)
   %    if ~isempty(strfind(g_grind.commands{g_l_i},'defextern '))
   %        error('GRIND:model:SyntaxError','It is not allowed to use "defextern" in the parameter panel, move it to model definition');
   %    end
  %    for g_l_i = 1:size(g_grind.commands, 2)
         if ~isempty(strtrim(g_grind.commands{g_l_i}))
            g_l_s = strtrim(char(g_grind.commands{g_l_i}));
            f1 = strfind(g_l_s, '%');
            f2=strfind(g_l_s, '=');
            if ~isempty(f2) && (isempty(f1) || f1(1) > f2(1))
               g_l_val = strtrim(g_l_s(f2(1):end));
               f1 = strfind(g_l_val, '%');
               if ~isempty(f1)
                  g_l_val=g_l_val(1:f1(1));
               end
               if ~isempty(str2num(g_l_val)) %#ok<ST2NM>
                  evalin('base', g_l_s);
               end
            end
         end
   %   end
   end
   msg='';
   for g_l_i = 1:size(g_grind.commands, 2)
      if ~isempty(strtrim(g_grind.commands{g_l_i}))
         %   g_l_s = [g_l_s char(g_grind.commands{g_l_i}) g_l_linefeed];
         g_l_s = char(g_grind.commands{g_l_i});
         %some changed global variables
%          g_l_s=strrep(g_l_s,'g_Jacobian','g_grind.syms.Jacobian');
%          g_l_s=strrep(g_l_s,'g_ndays','g_grind.ndays');
%          g_l_s=strrep(g_l_s,'g_tstep','g_grind.tstep');
%          g_l_s=strrep(g_l_s,'g_solver.drawnow','g_grind.drawnow');
%          g_l_s=strrep(g_l_s,'g_pen','g_grind.pen');
%          g_l_s=strrep(g_l_s,'g_solver.diffto','g_grind.diffto');
%          g_l_s=strrep(g_l_s,'g_statevars','g_grind.statevars');
         try
            evalin('base', g_l_s);
         catch err
            msg=err.message;
          %  disp(err.message);
            if strcmp('GRIND:loaddata:NoFile',err.identifier)
               g_l_ok = 1;
               fprintf(2,'Warning: Error reading data in parameter initialisation, line %d:\n%s\n', g_l_i, g_l_s);
               fprintf(2,'%s\n',err.message);
            else
               g_l_ok = 0;
               fprintf(2,'line %d: %s\n', g_l_i, g_l_s);
%               fprintf(2,'Error in parameter initialisation, line %d:\n%s\n', g_l_i, g_l_s);
%               fprintf(2,'%s\n',err.message);
            end
         end
      end
   end
   if ~g_l_ok
      error('GRIND:resetpars:ParInitError','Error in parameter initialisation (line %d): %s',g_l_i,msg);
   else
      n = 0;
      for i = 1:length(g_grind.pars)
         x = evalin('base', char(g_grind.pars{i}));
         if isempty(x)
            n = n + 1;
         end
      end
      if n > 0
         %parameter initiations have wrong crosslinks
         for i = 1:3
            evalin('base', g_l_s);
         end
      end
   end
   %  try
   %     evalin('base', g_l_s);
   %  catch
   %     for g_l_i = 1:size(g_grind.commands, 2)
   %
   %  end
end
N0=i_initvar;
g_grind.solver.iscomplex=any(~isreal(N0));
if (nargin==0) || ~silent
   disp('Parameters reset to the default values:');
   par;
end

% function parno = i_parno(apar)
% global g_grind;
% parno = [];
% f =  strfind(apar, '(');
% if ~isempty(f)
%    apar = apar(1:f(1) - 1);
% end
% for k = 1:size(g_grind.pars, 2)
%    if strcmp(apar, char(g_grind.pars{k}))
%       parno = k;
%       return;
%    end
% end
