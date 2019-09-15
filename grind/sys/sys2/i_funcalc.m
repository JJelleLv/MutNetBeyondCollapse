function [g_l_ran,g_l_res,g_s_s11] = i_funcalc(g_l_afun,g_l_avar,g_l_axrange,g_npoints,g_loc)
%function funplot(afun,g_l_avar,[exchange],[axrange, ayrange],logxy, npoints);
%
global g_grind g_Y; 
if nargin>4
   evalin('base',g_loc.comm);
end

eval(i_globalstr({g_l_avar}));
if isfield(g_grind, 'pars')
   if (~isempty(g_grind.pars))
      eval(i_globalstr(g_grind.pars));
   end

end

g_l_filter={'pi','inf','Inf','nan','NaN','eps'};
if ~isempty(g_grind)
   if isfield(g_grind, 'funcs')
      g_l_funcs = str2cell(g_grind.funcs);
      for g_l_i = 1:length(g_l_funcs)
         g_l_f=strfind(g_l_funcs{g_l_i}, '=');
         if ~isempty(g_l_f)
            g_l_filter = [g_l_filter, {strtrim(g_l_funcs{g_l_i}(1:g_l_f - 1))}];
         end

      end

   end

end

g_s_s11 = i_symvar(g_l_afun, g_l_filter);
if ~isempty(g_s_s11)
   for g_l_i = length(g_s_s11):-1:1
      if strcmp(g_s_s11{g_l_i}, g_l_avar)
         g_s_s11 = {g_s_s11{1:g_l_i - 1} g_s_s11{g_l_i + 1:length(g_s_s11)}};
      end

   end

   g_l_prompt = {};
   g_l_unknownpars = {};
   if ~isempty(g_s_s11)
      g_s_s12 = i_globalstr(g_s_s11);
      evalin('base', g_s_s12);
      eval(g_s_s12);
      for g_l_i = 1:length(g_s_s11)
         if eval(sprintf('isempty(%s)', g_s_s11{g_l_i}))
            g_l_prompt = [g_l_prompt {sprintf('Enter a value for "%s"', g_s_s11{g_l_i})}];
            g_l_unknownpars = [g_l_unknownpars g_s_s11(g_l_i)];
         end

      end

   end

   if ~isempty(g_l_prompt)
      Answer = inputdlg(g_l_prompt, 'Unknown parameter(s)', 1);
      for g_l_i = 1:length(g_l_prompt)
         assignin('base', g_l_unknownpars{g_l_i}, str2num(Answer{g_l_i}));  %#ok<ST2NM>
      end

   end

end

if isfield(g_grind, 'statevars')
   if ~isempty(g_grind.statevars) && g_grind.statevars.dim>0
      try
         NNN0 = i_initvar;
      catch  %#ok<CTCH>
         NNN0 = [];
      end

      if ~isempty(NNN0) && ~g_grind.statevars.vector || (g_grind.statevars.dims{1}.dim2 == 1)
         if g_grind.statevars.vector
            eval([g_grind.statevars.vectnames{1} ' = g_Y;']);
         else
            for g_l_i = 1:g_grind.statevars.dim
               eval([char(g_grind.statevars.names{g_l_i}) '= NNN0(' num2str(g_l_i) ');']);
            end

         end
      end

   end

end

g_l_ran = g_l_axrange(1):(g_l_axrange(2) - g_l_axrange(1))  / g_npoints:g_l_axrange(2);
if exist(g_l_avar,'var') && ~isempty(g_l_avar)
   eval(['g_l_old = ' g_l_avar ';']);
else
   g_l_old = []; 
end

try
   eval([g_l_avar '=g_l_ran;']);
   if isfield(g_grind, 'funcs')
      if ~isempty(g_grind.funcs)
         allfuncs=str2cell(g_grind.funcs);
         if ~strcmp(g_l_avar, 't')
            t = 1; 
            %temporary variable t
         end
         funcs=allfuncs(~strcmp(g_l_avar,g_grind.funcnames.names));
      %   g_l_s1 = g_grind.funcs; 
         try
            i_evalarray(funcs);
         catch %#ok<CTCH>
         end

      end

   end

   %set all operands to array operands
  % g_l_s1=['g_l_res=' g_l_afun ';']; 
   g_l_res=i_evalarray(g_l_afun);
   eval([g_l_avar '=g_l_old;']);
catch err
%   err=lasterror;
   eval([g_l_avar '=g_l_old;']);
   rethrow(err);
end

%++++++++++FUNCTIONS++++++++++
%function [x, y] = exchange(y, x)


