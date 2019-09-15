%FINDEQ   Find nearest equilibrium
%   Find the closest initial condition with zero growth. FINDEQ
%   uses the simplex optimizing method to find an initial
%   condition with zero growth.
%
%   Usage:
%   FINDEQ - use of FINDEQ command with default options.
%   FINDEQ -options - Display a dialog box with all options.
%   FINDEQ('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:      
%     'complex' [logical] - Search an equilibrium if state variable(s) can be complex numbers.
%     'display' [off | on | iter-detailed | iter | final | final-detailed | notify | notify-detailed] - Different display options of fminsearch
%     'keep' [logical] - Keep the last value
%     'mask' [logical and length(logical)==ndim or empty] - Mask for state variables [1 0 1 0] (1=optimize 0=keep fixed) (you can find a nullcline value or partial equilibrium
%    this way, leave empty for no mask)
%     'maxfunevals' [integer>0] - Max number of function evaluations
%     'maxiter' [integer>0] - Max number of iterations
%     'maxtime' [number>0|isnan(number)] - Maximum time for calculation in sec. NaN = ignore.
%     'method' [fminsearch | fminunc | fsolve | trust-region-dogleg | trust-region-reflective | levenberg-marquardt] - Algorithm of optimization
%    fminsearch or fsolve (for fsolve, trust-region or levenberg-marquardt is the optimization toolbox <a href="matlab:commands toolboxes">required</a>)
%     'n0' [number] - Initial values of state variables a vector of all state variables leave empty for current values []
%     'p0' [number] - Initial values of the parameters (vector of all values) []
%     'tolfun' [number>0] - Tolerance of function
%     'tolx' [number>0] - Tolerance of X
%     'userfun' [function_handle] - Add a user function for optimizing []
%   FINDEQ('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-default' - get default options as struct (used internally)
%     '-o' - options: display a user interface with options-
%     '-test' - test of consistency of findeq between the methods fsolve and fminsearch.
%     
%   Example:
%  [N1, isfound]=FINDEQ('Display','off') - returns the initial conditions and 
%     a flag if an equilibrium is found. The arguments ('Display','off')
%     suppresses output.
%
%   See also eigen, perturb, null, findeqs, conteq
%
%   Reference page in Help browser:
%      <a href="matlab:commands('findeq')">commands findeq</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function [N1, found] = findeq(varargin)
global g_Y g_grind;
if isfield(g_grind,'solver')
    iscomplex=g_grind.solver.iscomplex;
else
    iscomplex=false;
end
fieldnams={'complex', 'l', 'Search an equilibrium if state variable(s) can be complex numbers.', iscomplex;...
   'display', 'e[off|on|iter-detailed|iter|final|final-detailed|notify|notify-detailed]', 'Different display options of fminsearch','notify';...
   'keep', 'l', 'Keep the last value',true;...
   'mask', 'l&length(l)==ndim#E', 'Mask for state variables [1 0 1 0] (1=optimize 0=keep fixed) (you can find a nullcline value or partial equilibrium',[];...
   'maxfunevals', 'i>0', 'Max number of function evaluations',50000;...
   'maxiter', 'i>0', 'Max number of iterations',50000;...
   'maxtime', 'n>0|isnan(n)', 'Maximum time for calculation in sec. NaN = ignore.',NaN;...
   'method', 'e[fminsearch|fminunc|fs+olve|trust-region-dogleg|trust-region-reflective|l+evenberg-marquardt]', 'Algorithm of optimization','fminsearch';...
   'n0', 'n', 'Initial values of state variables a vector of all state variables leave empty for current values []',[];...
   'p0', 'n', 'Initial values of the parameters (vector of all values) []',[];...
   'tolfun', 'n>0', 'Tolerance of function', 1E-6;...
   'tolx', 'n>0', 'Tolerance of X',1E-10;...
   'userfun', 'f', 'Add a user function for optimizing []',[]}';

[args,def]=i_parseargs(fieldnams,'','-test,-o,-default',varargin);
if any(strcmp(args.opts, '-default')) %get the default settings
   %N1=def;
   N1=def;%struct('tolfun', 1E-6,'tolx',1E-10,'maxfunevals', 50000,'maxiter',50000,'maxtime',NaN,'method','fminsearch','keep',true,'complex',g_grind.solver.iscomplex,'display','notify','mask',[],'n0',[]);
   return;
end

if ~isfield(g_grind, 'findeq')
   g_grind.findeq= def;
end
%defaults per run
g_grind.findeq.keep =true;
g_grind.findeq.display  = 'notify'; % off | iter | iter - detailed | notify | notify - detailed | final | final - detailed
g_grind.findeq.mask  = [];
g_grind.findeq.n0  = [];
g_grind.findeq.userfun=[];
if isfield(args,'userfun')
    g_grind.findeq.userfun=args.userfun;
end
norun=false;
if isfield(args,'maxiter')&&args.maxiter==1
    norun=true;
    args=rmfield(args,'maxiter');
end
g_grind.findeq=mergestructs(g_grind.findeq,rmfield(args,'opts'));
if any(strcmp(args.opts, '-test'))
   %test of consistency of findeq with fsolve and fminsearch
   N0 = i_initvar;
   N1=findeq('method','fminsearch','keep',false,'display','off','n0',N0);
   N2=findeq('method','fsolve','keep',false,'display','off','n0',N0);
   if max(abs(N2-N1)) > 1E-4
      error('grind:findeq','Findeq Test failed: difference between fsolve and fminsearch is %g',max(abs(N2-N1)));
   else
       disp('findeq test ok')
   end

   if length(N0) > 1
      mask = false(size(N0));
      mask(1) = true;
      N1=findeq('Method','fminsearch','Mask',mask,'Keep',false,'Display','off');
      N2=findeq('Method','fsolve', 'Mask',mask,'Keep',false,'Display','off');
      if max(abs(N2-N1)) > 1E-4
         error('grind:findeq','Findeq test with mask failed: difference between fsolve and fminsearch is %g',max(abs(N2-N1)));
      else
          disp('findeq test (mask) ok')
      end

   end

elseif any(strcmp(args.opts, '-o'))
   args=i_parseargdlg('findeq',fieldnams,'-test,-o,-default',g_grind.findeq);
   if ~isempty(args)
       findeq(args);
   end
%    prompt={'Tolerance of function (TolFun):', 'Tolerance of X (TolX)',...
%       'Max number of function evaluations (MaxFunEvals)','Max number of iterations (MaxIter)','Max time in sec/NaN=ignore (MaxTime)',...
%       'Method of optimization (Method=fminsearch | fsolve | trust-region-reflective | levenberg-marquardt)','Display results (Display= off | on |iter)',...
%       'Mask for state variables (Mask []= no mask)','Initial values of state variables (NO []= current', 'Keep new initial conditions (Keep=yes|no)',...
%       'Allow complex numbers (complex=yes|no)'};
%    answers = {num2str(g_grind.findeq.tolfun), num2str(g_grind.findeq.tolx), num2str(g_grind.findeq.maxfunevals), num2str(g_grind.findeq.maxiter), num2str(g_grind.findeq.maxtime), ...
%       g_grind.findeq.method, g_grind.findeq.display, mat2str(g_grind.findeq.mask), mat2str(g_grind.findeq.n0),iif(g_grind.findeq.keep,'yes','no'),iif(g_grind.findeq.complex,'yes','no')};
% 
%    answers = inputdlg(prompt, 'findeq options', 1, answers);
%    if isempty(answers)
%       return;
%    end
% 
%    findeq('tolfun',answers{1},'tolx',answers{2},'maxfunevals',answers{3}, 'maxiter', answers{4}, 'maxtime', answers{5},'method',answers{6},'display', answers{7},...
%    'mask' ,answers{8},'n0', answers{9},'keep' , answers{10},'complex',answers{11});
   return
end

if g_grind.solver.nonautonomous&&strcmp(g_grind.findeq.display,'on')
   if strcmp(questdlg('Might not find correct equilibria in nonautonomous equation, continue?','Error','Yes','No','No'),'No')
      error('GRIND:findeq:NonAutonomous','Make equation autonomous by adding tau''=1 and replacing all t by tau');
   end

end
if isfield(args,'p0')
    oldp=par(struct('opts',{{'-vector'}}));
    setpars(args.p0);
else
    oldp=[];
end
dokeep = g_grind.findeq.keep;
display = ~strcmpi(g_grind.findeq.display, 'off');
if isempty(g_grind.findeq.n0)
   N0 = i_initvar;
else
   N0 = g_grind.findeq.n0(:);
end

if any(strcmp(g_grind.findeq.method , {'fsolve','trust-region-dogleg','trust-region','levenberg-marquardt','fminunc'}))&&~i_hastoolbox('optim')
   warning('grind:findeq','Optimize toolbox is needed for the method "%s", fminsearch is used instead',g_grind.findeq.method);
   g_grind.findeq.method = 'fminsearch';
end

g_grind.hfun.curr = i_getodehandle(1, ''); 
if ~isnan(g_grind.findeq.maxtime)&&~isempty(g_grind.findeq.maxtime)
    tstart=tic;
    maxtime=g_grind.findeq.maxtime;
    outf=@(x, optimValues, state)timefun(maxtime,tstart,x, optimValues, state);
else
    outf=[];
end
opts=optimset('tolfun' , g_grind.findeq.tolfun,'tolx',g_grind.findeq.tolx, 'maxfunevals' ,g_grind.findeq.maxfunevals,'maxiter',g_grind.findeq.maxiter,'display' ,g_grind.findeq.display, 'Outputfcn',outf);
if ~isempty(g_grind.findeq.mask)&&~all(logical(g_grind.findeq.mask))
   lastN0 = g_grind.findeq.n0;
   g_grind.findeq.mask = logical(g_grind.findeq.mask(:));
   g_grind.findeq.n0 = N0;
   N0 = N0(g_grind.findeq.mask);
   if norun
       NRes1=N0;
       fun=i_getodehandle(5,'');
       eqfound=fun(NRes1)<1E-4;
   elseif any(strcmp(g_grind.findeq.method , {'fsolve','trust-region-dogleg','trust-region','levenberg-marquardt'}))
       if ~strcmp(g_grind.findeq.method ,'fsolve')
           opts=optimset(opts,'Algorithm',g_grind.findeq.method);
       end
       if g_grind.findeq.complex
          [NRes1, ~, eqfound] = fminsearch(i_getodehandle('findeq_complex',''), [real(N0);imag(N0)], opts);
          NRes1=complex(NRes1(1:end/2),NRes1(end/2+1:end));
       else
          [NRes1, ~, eqfound] = fsolve(i_getodehandle(6,''), N0, opts);
       end
   else
       if ~strcmp(g_grind.findeq.method ,'fminunc')
          [NRes1, ~, eqfound] = fminsearch(i_getodehandle(5,''), N0, opts);
       else
          [NRes1, ~, eqfound] = fminunc(i_getodehandle(5,''), N0, opts);
       end
   end
   NRes = g_grind.findeq.n0;
   NRes(g_grind.findeq.mask) = NRes1;
   g_grind.findeq.n0 = lastN0;
else
   if norun
       NRes=N0;
       fun=i_getodehandle(5,'');
       res=fun(NRes);
       eqfound=res<1E-4;
   elseif any(strcmp(g_grind.findeq.method , {'fsolve','trust-region-dogleg','trust-region','levenberg-marquardt'}))
      if ~strcmp(g_grind.findeq.method ,'fsolve')
           opts=optimset(opts,'Algorithm',g_grind.findeq.method);
       end
      [NRes, res, eqfound] = fsolve(i_getodehandle(6,''), N0, opts);
   else
       if g_grind.findeq.complex
           [NRes, res, eqfound] = fminsearch(i_getodehandle('findeq_complex',''), [real(N0);imag(N0)], opts);
           NRes=complex(NRes(1:end/2),NRes(end/2+1:end));
       else
           if ~strcmp(g_grind.findeq.method ,'fminunc')
           [NRes, res, eqfound] = fminsearch(i_getodehandle(5,''), N0, opts);
           else
               warning off optim:fminunc:SwitchingMethod
           [NRes, res, eqfound] = fminunc(i_getodehandle(5,''), N0, opts);
           warning on optim:fminunc:SwitchingMethod
           end
       end
   end
   eqfound = all(res < g_grind.findeq.tolfun * 2) & eqfound; %No convert
end

if norun||eqfound
   g_Y = [];
   smallvalues = abs(NRes) < 1E-6;
   if any(smallvalues)&&~g_grind.solver.isimplicit %%&&~g_grind.solver.isdiffer
      %for border equilibria: is the result better if smallvalues are zero?
      %-is usually more precise then
      NRes1 = NRes;
      f = find(smallvalues);
      g_grind.hfun.curr = i_getodehandle(1, '');
      for i = 1:length(f)
         f1 = f(i);
         NRes1(f1) = 0;
         r1 = feval(g_grind.hfun.curr, 1, NRes);
         r2 = feval(g_grind.hfun.curr, 1, NRes1);
         if abs(r2(f1)) < abs(r1(f1))&&(sum(r2.^2) < g_grind.findeq.tolfun * 2)
            NRes = NRes1;
         else
            NRes1(f1) = NRes(f1);
         end

      end

   end

   if display
      N_old = i_initvar;
      i_keep(NRes);
      if nargout==0&&~norun&&isempty(oldp)
         findeqs('-add_curr');
      end
      disp('Equilibrium found:');
      try
         val();
         [~, eigenval] = i_eigen(isempty(g_grind.syms.Jacobian),g_grind.solver.iters);
         [stable, issaddle] =  i_stability(eigenval, g_grind.solver.isdiffer);
         if stable
            h=where('sblink','o', 'k');
         else
            if issaddle
               h=where('sblink','o',[0.8 0.8 0.8]);
            else
               h=where('sblink','o','w');
            end

         end
         set(h, 'markersize', 6)
      catch %#ok<CTCH>
         disp('Could not evaluate the stability');
         where('o')
      end
      if ~dokeep
         i_keep(N_old);
      end

   elseif dokeep
      i_keep(NRes);
      if nargout==0&&~norun&&isempty(oldp)
         findeqs('-add_curr');
      end
   end

elseif display
   i_errordlg('Could not find an equilibrium, try other initial conditions');
end
if ~isempty(oldp)
    setpars(oldp);
end
if nargout > 0
   N1 = NRes;
   found = eqfound;
end
function setpars(p0)
global g_grind
if ~isempty(p0)
    ndx=1;
    for i=1:length(g_grind.pars)
        siz=evalin('base',sprintf('size(%s);',g_grind.pars{i}));
        assignin('base',g_grind.pars{i},reshape(p0(ndx:ndx+prod(siz)-1),siz(1),siz(2)));
        ndx=ndx+prod(siz);
    end
end

function stop=timefun(maxtime,tstart,~, ~, ~) 
stop=toc(tstart)>maxtime;
    


