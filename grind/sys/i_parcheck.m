function [Res, S] = i_parcheck(nodlg)
global g_grind;
if nargin == 0
   nodlg = 0;
end
res = 1;
s = '';
startinggrind = isfield(g_grind, 'starting');
LF =  sprintf('\n');  %#ok<SPRINTFN>

if isfield(g_grind, 'pars') && ~startinggrind
   if any(strcmp('eps', g_grind.pars))
      i_warningdlg('GRIND:model:eps','"eps" is a reserved word in MATLAB, do not use for parameter name');
   end
   for i = 1:length(g_grind.pars)
      x = evalin('base', char(g_grind.pars{i}));
      if isempty(x)
         res = 0;
         s = sprintf('%sParameter %s is not initialized\n', s, char(g_grind.pars{i}));
      end
   end
end
if res
   if isempty(g_grind)
      s = 'No model defined, use <a href="matlab: model">model</a> or <a href="matlab:vismod">vismod</a> to define a model';
      res = 0;
   elseif ~isfield(g_grind, 'statevars') || isempty(g_grind.statevars) || ...
         (isempty(g_grind.statevars.names)&&isempty(g_grind.statevars.vectnames))
      s = 'No state variables defined';
      res = 0;
   else
      if ~g_grind.solver.isdiffer
         g_grind.solver.iters = 1;
      end
      %g_grind.solver.opt.OutputFcn = []; %str2func('i_odespeed');
      if ~isempty(g_grind.odefile) && ~startinggrind
         try
            N0 = i_initvar;
         catch
            if g_grind.statevars.vector
               evalin('base', i_globalstr(g_grind.statevars.vectnames));
            else
               evalin('base', i_globalstr(g_grind.statevars.names));
            end
            N0 = i_initvar;
         end
         if ~isempty(g_grind.permanent)
            defpermanent('-deactivate');
         end
         if ~g_grind.solver.isimplicit&&~g_grind.checks.modelchecked
            if ~isempty(N0)
               if ~strcmp(g_grind.solver.name, 'csolvers')
                  %   disp('run parcheck');
                  feval(i_getodehandle(1,''), 1, N0);
               end
            else
               s = 'State variables not initialized';
               res = 0;
            end
         end
      end
      for i = 1:length(g_grind.externvars)
         if isempty(g_grind.externvars{i}.data)
            if ~isfield(g_grind.externvars{i}, 'nodatawarning')||~g_grind.externvars{i}.nodatawarning
               warning('GRIND:externvars:nodata','Some external variables have no data, use <a href="matlab: setdata">setdata</a> to add data');
               g_grind.externvars{i}.nodatawarning = 1;
            end
            break
         end
      end
   end
end
if ~res
   if ~nodlg
      s = ['****** ERROR in GRIND ******' LF s];
      %     errordlg(s);
      disp(' ');
      error('GRIND:parcheck:ModelNotCorrect', s);
   else
      disp('Model not correct:');
      disp(s);
   end
else
   %check if the vectorized model gives the same result as the unvectorized
   %if not (or if error), turn Vectorized off
   %vectorized is much faster in nullcline plots and in numerical jacobians
   %
   if ~g_grind.checks.modelchecked&&strcmp(g_grind.solver.opt.Vectorized, 'on')
      N0 = i_initvar;
      N0 = [N0 N0 * 1.1 N0 * 1.2];
      t1 = [1 1 1];
      if isempty(g_grind.pars)
          p0=zeros(0,3);
      else
          s = sprintf('%s(:);', g_grind.pars{:});
          p0 = evalin('base',sprintf('[%s];',s));
          p0 = [p0 p0 * 1.1 p0 * 1.2];
      end
      odefun = i_getodehandle('vectorcheck');
      try
         N1v0= odefun(t1, N0,p0(:,1));
         N1v = odefun(t1, N0, p0);
         N1nv = zeros(size(N0));
         for i = 1:3
            N1nv(:, i) = odefun(t1(i), N0(:, i), p0(:, i));
         end
         if ~all(size(N1v) == size(N1nv))||~all(size(N1v0)==size(N0))
             g_grind.solver.opt.Vectorized = 'off';
             warning('GRIND:parcheck:Vectorized','Cannot optimize speed of odefile, vectorized model gives too few results');
         else
            if ~(max(max(abs(N1v-N1nv))) < 1E-10||g_grind.solver.isstochastic)
               g_grind.solver.opt.Vectorized = 'off';
               warning('GRIND:parcheck:Vectorized','Cannot optimize speed of odefile, vectorized model results are different');
            end
         end
      catch
         g_grind.solver.opt.Vectorized = 'off';
         warning('GRIND:parcheck:Vectorized','Cannot optimize speed of odefile, error running vectorized model');
      end
   end
   g_grind.checks.modelchecked = true;
end
if nargout > 0
   Res = res;
   S = s;
end
function s1 = getline(file, line, col)
fid = fopen(file, 'r');
for i = 1:line
   s1 = fgetl(fid);
end
s1 = [s1(1:col) '#' s1(col + 1:end)];
s1 = removegX(s1);
col = strfind(s1, '#');
s1=strrep(s1,'#','');
s1=sprintf('%s\n%s^\n',s1,repmat(' ',1,col));
fclose(fid);


function s1 = removegX(s)
global g_grind;
if ~isempty(g_grind) && isfield(g_grind, 'statevars') && ~isempty(g_grind.statevars)
   s1 = s;
   for d = 1:length(g_grind.statevars.names)
      inp = sprintf('g_X1(%d)', d);
      s1 = i_changevar(s1, inp, g_grind.statevars.names{d});
      inp = sprintf('g_X2(%d,1)',d);
      s1=i_changevar(s1,inp,[g_grind.statevars.names{d} '''']);
   end
end

