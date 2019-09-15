%OPEN_MATCONT_GUI   Open MatCont with the current model
%   This function provides an interface to <a href="http://www.matcont.ugent.be/">MatCont</a>, an advanced tool for the 
%   numerical continuation of bifurcations of dynamical systems. Leaders of the project are 
%   Willy Govaerts (Gent,B) and Yuri A. Kuznetsov (Utrecht,NL). This interface works only if
%   MatCont is installed (tested with MatCont version 6p3). Currently only ODE models
%   are supported (no matrix notation).
%   Note that:
%       *Parameters/variables named I,i,c or C are renamed to I1,i1,c1 or C1.
%       *Non-autonomous models get an extra statevariable tau (tau'=1).
%       *External variables become parameters with the default value.
%       *MatcontM (for maps) is not supported.
%       
%  
%   Usage:
%   OPEN_MATCONT_GUI - the model is translated to MATCONT (as far as possible) and the 
%   system editor of MatCont is opened. Check first if the conversion succeeded.
%  
%   See also conteq, paranal
%
%   Reference page in Help browser:
%      <a href="matlab:commands('open_matcont_gui')">commands open_matcont_gui</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function open_matcont_gui(varargin)
global g_grind;
i_parseargs('','','',varargin);
matcontpath =fileparts(which('matcont.m'));
if isempty(matcontpath)
    matcontpath = findgrindfile('matcont.m');
    if isempty(matcontpath)
        matcontpath = uigetdir(matcontpath, 'Find the directory of MATCONT');
    end
end
findeq;
if g_grind.statevars.vector
   error('grind:matcont','matcont does not support vector models');
end

if g_grind.solver.isdiffer
   error('grind:matcont','matcontm (for maps) not yet supported');
end

themodel = str2cell(g_grind.funcs);
for i = 1:length(g_grind.model)
   f1=strfind(g_grind.model{i},'''');
   feq=strfind(g_grind.model{i}, '=');
   if ~isempty(f1)&&~isempty(feq)&&f1(1) < feq(1)
      themodel(end + 1, 1) = g_grind.model(i);
   end

end



pars = g_grind.pars;
statevars = g_grind.statevars.names;
allvars = [g_grind.pars, g_grind.statevars.names];
if g_grind.solver.nonautonomous
   var1 = 'tau1';
   if any(strcmp(allvars, var1))
      var1 = 'tau2';
   end

   for i = 1:length(themodel)
      fs = parsed_equation(themodel{i}).fs;
      fs(strcmp(fs, 't')) = {var1};
      themodel{i} = sprintf('%s', fs{:});
   end

   themodel{end+1}=sprintf('%s''=1',var1);
   statevars{end + 1} = var1;
end

[themodel, pars] = checkvars(themodel, pars, allvars);
[themodel, statevars] = checkvars(themodel, statevars, allvars);


gds.coordinates = cell(length(statevars), 2);
for i = 1:length(statevars)
   gds.coordinates(i, 1) = statevars(i);
   if ~any(strcmp(statevars{i}, allvars))
      gds.coordinates{i, 2} = 0;
   else
      gds.coordinates{i, 2} = evalin('base', g_grind.statevars.names{i});
   end

end

gds.parameters = cell(length(pars), 2);
for i = 1:length(pars)
   gds.parameters(i, 1) = pars(i);
   gds.parameters{i, 2} = evalin('base', g_grind.pars{i});
end

if ~isempty(g_grind.externvars)
   warning('grind:matcont','External variables are not supported, changed them to auxilliary variables')
  % npar = length(pars);
   for i = 1:length(g_grind.externvars)
      %       gds.parameters{npar+i,1}=g_grind.externvars{i}.name;
      %       gds.parameters{npar+i,2}=str2double(g_grind.externvars{i}.default);
      comm=sprintf('%s=externvar(', g_grind.externvars{i}.name);
      f = find(strncmp(themodel, comm, length(comm)));
      if length(f) == 1
         themodel{f}=sprintf('%s=%s;', g_grind.externvars{i}.name, g_grind.externvars{i}.default);
      end

   end

end
gds.time = cell(1, 2);
if strcmp(g_grind.diffto, 'Time (t)')
   tvar = 't';
else
   tvar = g_grind.diffto;
end

gds.time{1} = tvar;
gds.time{2} = evalin('base', tvar);
gds.options.InitStepsize = [];
gds.options.MinStepsize = [];
gds.options.MaxStepsize = [];
gds.options.MaxCorrIters = [];
gds.options.MaxNewtonIters = [];
gds.options.MaxTestIters = [];
gds.options.MoorePenrose = [];
gds.options.SymDerivative = [];
gds.options.SymDerivativeP = [];
gds.options.Increment = [];
gds.options.FunTolerance = [];
gds.options.VarTolerance = [];
gds.options.TestTolerance = [];
gds.options.Singularities = 1;
gds.options.MaxNumPoints = [];
gds.options.Backward = 0;
gds.options.CheckClosed = [];
gds.options.TestFunctions = [];
gds.options.WorkSpace = [];
gds.options.Locators = [];
gds.options.Adapt = [];
gds.options.IgnoreSingularity = [];
gds.options.ActiveParams = [];
gds.options.Multipliers = 1;
gds.options.Eigenvalues = 1;
gds.options.Userfunctions = [];
gds.options.UserfunctionsInfo = [];
gds.options.PRC = 0;
gds.options.dPRC = 0;
gds.options.Input = [];
gds.options.ActiveUParams = [];
gds.options.ActiveSParams = [];
gds.options.ActiveSParam = [];
[~, n] = fileparts(g_grind.inifile);
gds.system = n;
gds.curve.new = 'P_O(1)';
gds.curve.old = '';
gds.equations = sprintf('%s\n', themodel{:});
gds.equations=strrep(gds.equations,'./','/');
gds.equations=strrep(gds.equations,'.*','*');
gds.equations=strrep(gds.equations,'.^','^');
ndx=~(gds.equations==' '|gds.equations==';');
gds.equations = gds.equations(ndx);
gds.dim = g_grind.statevars.dim;
if i_hastoolbox('symbolic')
   %symbolic Jacobians are recommended
   gds.der = [0 0 0 1 1; 0 0 0 0 0; 0 0 0 0 0; 1 1 1 0 0];
else
   gds.der = [1 1 1 1 1; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
end

gds.jac = '';
gds.jacp = '';
gds.hess = {};
gds.hessp = {};
gds.tensor3 = {};
gds.tensor4 = {};
gds.tensor5 = {};
gds.point = 'P ';
gds.type = 'O ';
gds.discretization.ntst = 20;
gds.discretization.ncol = 4;
gds.period = 1;
gds.plot2 = '';
gds.plot3 = '';
gds.PRC = '';
gds.dPRC = '';
gds.open.figuur = 1;
gds.open.continuer = 0;
gds.open.numeric_fig = 0;
gds.open.D2 = 0;
gds.open.D3 = 0;
gds.open.PRC = 0;
gds.open.dPRC = 0;
gds.open.integrator = 1;
gds.integrator.method = g_grind.solver.name;
gds.integrator.options.AbsTol = g_grind.solver.opt.AbsTol;
gds.integrator.options.BDF = [];
gds.integrator.options.Events = [];
gds.integrator.options.InitialStep = [];
gds.integrator.options.Jacobian = @jacobian;
gds.integrator.options.JacobianP = [];
gds.integrator.options.Hessians = @hessians;
gds.integrator.options.HessiansP = [];
gds.integrator.options.Der3 = [];
gds.integrator.options.Der4 = [];
gds.integrator.options.Der5 = [];
gds.integrator.options.JConstant = [];
gds.integrator.options.JPattern = [];
gds.integrator.options.Mass = [];
gds.integrator.options.MassConstant = [];
gds.integrator.options.MassSingular = [];
gds.integrator.options.MaxOrder = [];
gds.integrator.options.MaxStep = [];
gds.integrator.options.NormControl = [];
gds.integrator.options.OutputFcn = @integ_prs;
gds.integrator.options.OutputSel = [];
gds.integrator.options.Refine = [];
gds.integrator.options.RelTol = g_grind.solver.opt.RelTol;
gds.integrator.options.Stats = [];
gds.integrator.options.Vectorized = [];
gds.integrator.tspan = [0; g_grind.ndays];
solvername = g_grind.solver.name;
if any(strcmp(solvername,{'euler','c.ode45','c.euler','c.rk4','rk4','back_euler'}))
   solvername = 'ode45';
end

gds.integrator.method = solvername;
% gds.numeric.O=cell(3,2);
% gds.numeric.EP=cell(5,2);
% gds.numeric.LC=cell(8,2);
% gds.numeric.PD=cell(5,2);
% gds.numeric.Hom=cell(5,2);
% gds.numeric.HSN=cell(5,2);
gds.diagram = 'diagram';
gds.T = [];
gds.eps0 = [];
gds.eps1 = [];
gds.extravec = [1 0 0];
gds.t = [];
gds.epsilon = [];
oldcd = cd;
cd([matcontpath filesep 'Systems']);
save(gds.system, 'gds');
mfile = [gds.system '.m'];
%if ~exist(mfile,'file')
copyfile('standard.m', mfile);
%end

cd ..
model('-c');
matcont;
editsystem;
set(0,'showhiddenhandles','on');
ch = get(0, 'children');
hsystems=findobj(ch,'name','Systems');
hdata = guidata(hsystems);
mfiles = get(hdata.files, 'String');
f = find(strcmp(mfile, mfiles));
if isempty(f)
   f = 1;
end

set(hdata.files, 'value', f);
drawnow;
editsystem('edit_Callback',hsystems,'edit',guidata(hsystems))
cd(oldcd);
function [themodel, vars] = checkvars(themodel, vars, allvars)
if any(strcmpi(vars,'i')|strcmpi(vars,'c'))
   no = 1;
   warning('grind:matcont','Renamed some variables ("c" and "i" are not allowed in MATCONT)');
   while any(strcmpi(allvars,sprintf('i%d',no))|strcmpi(allvars,sprintf('c%d',no)))
      no = no + 1;
   end

   vars(strcmp(vars,'i'))={sprintf('i%d',no)};
   vars(strcmp(vars,'I'))={sprintf('I%d',no)};
   vars(strcmp(vars,'c'))={sprintf('c%d',no)};
   vars(strcmp(vars,'C'))={sprintf('C%d',no)};
   for i = 1:length(themodel)
      fs = parsed_equation(themodel{i}).fs;
      fs(strcmp(fs,'i'))={sprintf('i%d',no)};
      fs(strcmp(fs,'I'))={sprintf('I%d',no)};
      fs(strcmp(fs,'c'))={sprintf('c%d',no)};
      fs(strcmp(fs,'C'))={sprintf('C%d',no)};
      themodel{i} = sprintf('%s', fs{:});
   end

end
