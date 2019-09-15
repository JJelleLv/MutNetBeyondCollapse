function i_makemodel(mode1, comman1, inifil1, schem)
% if (nargin>=3) || (nargin==1)
%    evalin('base','initgrind');
% end

global g_grind;
if nargin >= 3
   finishgrind;
   evalin('base','initgrind');
   g_grind.model = mode1;
   g_grind.commands = comman1;
   g_grind.inifile = inifil1;
   if nargin == 4
      g_grind.scheme = schem;
   elseif ~isfield('scheme', g_grind)
      g_grind.scheme = {};
   end

   clear('mode1','comman1','inifil1');
elseif nargin == 1
   finishgrind;
   evalin('base','initgrind');
   g_grind.model = mode1.model;
   g_grind.commands = mode1.commands;
   g_grind.inifile = mode1.inifile;
   if isfield(mode1, 'scheme')
      g_grind.scheme = mode1.scheme;
   end

end

%update odefile and g_grind based on g_grind.model
if isempty(g_grind)||isempty(g_grind.model)
    error('grind:makemodel','no model defined');
end


g_grind.hfun.curr=[];
g_grind.hfun.odehandle=[];
g_grind.figopts={};
clear(g_grind.odefile);
ggrind=i_analysemodel(transpose(g_grind.model));  %analysemodel is independent
if ~isempty(ggrind.errors)
    error('grind:model','%s\n',ggrind.errors{:});
end

g_grind.model=ggrind.model;
if isfield(ggrind.solver,'dwiener')
    g_grind.solver.dwiener=ggrind.solver.dwiener;
end
if isfield(ggrind.solver,'djump')
    g_grind.solver.djump=ggrind.solver.djump;
end
g_grind.solver.nonautonomous =ggrind.solver.nonautonomous;
g_grind.solver.isstochastic =ggrind.solver.isstochastic;
g_grind.solver.eulerneeded =ggrind.solver.eulerneeded;
g_grind.solver.isdiffer =ggrind.solver.isdiffer;
g_grind.solver.isimplicit =ggrind.solver.isimplicit;
g_grind.statevars = ggrind.statevars;
g_grind.solver.haslags=ggrind.solver.haslags;
if g_grind.solver.haslags
    ggrind.solver.opt.Vectorized='off';
    g_grind.dde=ggrind.dde;
    g_grind.dde.lags=i_enterjacdlg('str2jac',g_grind.dde.lags);
    g_grind.solver.history=ggrind.solver.history;
end

g_grind.funcnames.names=ggrind.funcnames.names;
g_grind.funcs=ggrind.funcs;
g_grind.modelndx=ggrind.modelndx;
g_grind.permanent={};
for i=1:length(ggrind.permanent)
  if isfield(ggrind.permanent{i},'currentval')
     defpermanent('-i_makemodel','var',ggrind.permanent{i}.name,'value',ggrind.permanent{i}.currentval);
  else
     defpermanent('-i_makemodel','var',ggrind.permanent{i}.name);
  end

end

if ~isempty(ggrind.odefile)
    g_grind.odefile=ggrind.odefile;
end
g_grind.solver.name=ggrind.solver.name;
g_grind.solver.opt=solver('-defaultopt');
g_grind.solver.opt.Vectorized=ggrind.solver.opt.Vectorized;
if ~isempty(ggrind.permanent)||~isempty(ggrind.externvars)
    g_grind.solver.opt.Vectorized='off';
end

g_grind.boxcar.names=ggrind.boxcar.names;
g_grind.boxcar.gcycl=ggrind.boxcar.gcycl;
if isfield(ggrind,'implicdisp')
    for i=1:length(ggrind.implicdisp)
        implicitdisperse(i,ggrind.implicdisp{i}.Name,ggrind.implicdisp{i}.D);
    end

end

g_grind.pars=ggrind.pars;
out('-defaults');
if isfield(ggrind,'theodefile')
    oldcd=cd;
    cd(grindpath);
    if strcontains(ggrind.theodefile(2),'curr_ode0(')
        ggrind.theodefile{2}=strrep(ggrind.theodefile{2},'curr_ode0',g_grind.odefile);
    end
    
    if ~strcontains(g_grind.odefile,'.')
        fid=fopen(sprintf('%s.m',g_grind.odefile),'w');
    else
        fid=fopen(g_grind.odefile,'w');
    end
    
    fprintf(fid,'%s\n',ggrind.theodefile{:});
    fclose(fid);
    cd(oldcd);
end
g_grind.hfun.odehandle=str2func(g_grind.odefile);
which(g_grind.odefile);
if ~isempty(g_grind.pars)
   evalin('base', i_globalstr(g_grind.pars,[],'clear')); % clear the parameters first to avoid warnings
   evalin('base', i_globalstr(g_grind.pars));
end

g_grind.externvars={};
for i=1:length(ggrind.externvars)
  if isfield(ggrind.externvars{i},'options')
     if length(ggrind.externvars{i}.options)==2
         s1= defextern(ggrind.externvars{i}.name,ggrind.externvars{i}.default,ggrind.externvars{i}.options{:});
     else
         s1= defextern(ggrind.externvars{i}.name,ggrind.externvars{i}.default,ggrind.externvars{i}.options{1});
     end

  else
     s1= defextern(ggrind.externvars{i}.name,ggrind.externvars{i}.default);
  end

  g_grind.funcs=sprintf('%s\n%s',s1,g_grind.funcs);
end

setevent('-c')
if ~isempty(g_grind.commands)
   g_grind.commands = removepoints(g_grind.commands, 1);
   for g_l_h = 1:size(g_grind.commands, 2)
      if isempty(strfind(g_grind.commands{g_l_h}, '%'))
         g_l_c = strtrim(g_grind.commands{g_l_h});
         g_l_l = length(g_l_c);
         while (g_l_l > 1) && (g_l_c(g_l_l) == ' ')
            g_l_l = g_l_l - 1;
         end

         if (g_l_l > 0) && (g_l_c(g_l_l) ~= ';')
            g_grind.commands{g_l_h} = [g_grind.commands{g_l_h} ';'];
         end

      end

   end

end

i_modelinit;
if g_grind.solver.isdiffer && g_grind.solver.nonautonomous
   g_grind.solver.nonautonomous = 0;
end

if isempty(g_grind.statevars.names) && isempty(g_grind.statevars.vectnames)
    cd(oldcd);
   error('GRIND:model:NoStatevars','Error: no state variables detected');
else
   if g_grind.statevars.vector
      evalin('base', i_globalstr(g_grind.statevars.vectnames,[],'clear')); %avoid warnings
      evalin('base', i_globalstr(g_grind.statevars.vectnames));
   else
      evalin('base', i_globalstr(g_grind.statevars.names,[],'clear')); %avoid warnings
      evalin('base', i_globalstr(g_grind.statevars.names));
   end

   i_initvar(1); %initiates the state variables
end

if ~isempty(g_grind.commands) && iscell(g_grind.commands{1})
   g_grind.commands = g_grind.commands{1};
end

try
   g_grind.starting = 1;
   resetpars(1);
   if isfield(g_grind, 'starting')
      g_grind = rmfield(g_grind, 'starting');
   end

catch err
   if isfield(g_grind, 'starting')
      g_grind = rmfield(g_grind, 'starting');
   end

   rethrow(err);
end

if ~isempty(ggrind.definespace)
    evalin('base',sprintf('%s\n',ggrind.definespace{:}));
end

function s = stripcomments(alist)
f=transpose(strfind(alist,'%'));
s=alist;
for i=1:length(alist)
    if ~isempty(f{i})
        s{i}=alist{i}(1:f{i}(1)-1);
    end

end

function alist = removepoints(alist, docompound)
h = 1;
i = 1;
siz = length(alist);
s = stripcomments(alist);
%s=alist;
sfull = alist;
while h <= length(alist)
   alist{i} = '';
   h2=stackedfindend(s,'[',']',h);
   if isempty(h2) || (h2==h)
      h2=stackedfindend(s,'{','}',h);
      if docompound && (isempty(h2) || (h2==h))
         h2=stackedfindend(s,{'for','if','while','switch','case'},'end',h,1);
      end

   end
   eoln=sprintf('\n'); %#ok<SPRINTFN> for older versions
   %merge h2 lines
   if ~isempty(h2)
      for j = h:h2 - 1
         alist{i} = [char(alist{i}) char(sfull{j}) eoln]; 
         h = h + 1;
      end

   end

   %merge ... lines
   while (h < siz) && (length(strfind(s{h}, '...')) == 1)
      alist{i} = [char(alist{i}) char(sfull{h}) eoln]; 
      h = h + 1;
   end

   alist{i} = [char(alist{i}) char(sfull{h})];
   h = h + 1;
   i = i + 1;
end

alist = alist(1: i - 1);

function h = stackedfindend(list, sstart, send, i1,wholewords)
if nargin < 4
   i1 = 1;
end

if nargin < 5
   wholewords = 0;
end

h = i1;
stack = 0;
if isempty(findstrings(sstart, list{h}, wholewords))
   h = [];
else
   found = 0;
   while (h < length(list)) && ~found
      stack = stack + length(findstrings(sstart, list{h},wholewords)) - length(findstrings(send, list{h},wholewords));
      found=stack == 0;
      if ~found
         h = h + 1;
      end

   end

end


