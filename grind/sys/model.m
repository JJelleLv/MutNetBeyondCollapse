%MODEL   Open/create/edit a model
%   Use the upper panel in the window to enter differential 
%   equation or difference equations (or iteration maps).
%   In the lower panel default values of the parameters and initial
%   values of the state variables are entered. 
%   Furthermore default commands can be entered here (for instance
%   axis limits, see <a href="matlab:help ax">ax</a>).
%   There are no restrictions for the number of parameters, functions and equations.
%
%   Note that parameters and state variables are case-sensitive.
%   Parameter names can be all alphanumeric names except the
%   following reserved words:
%   <a href="matlab:commands t">t</a>   - initial time.
%   pi  - pi=3.1416
%   inf - Infinity
%   Inf - Infinity
%   nan - NaN (not-a-number)
%   NaN - NaN (not-a-number)
%   eps - Floating point relative accuracy eps=2.2204e-016 
%
%
%   Examples:
%   Model equations:
%   cons=Z*A/(h+A)  - "function" with temporal results that can be used in the 
%   other equations (see <a href="matlab:help funcs">funcs</a>) Such function cannot take arguments.
%   N(t + 1) = N(t) * r * (1 - N(t) / K)  - Use (t+1) and (t) for difference equations 
%   Higher order equations (N(t + 2) ) are not allowed, but can be written as first
%   order equations.
%   N' = N * r * (1 - N / K) - Use ' for differential equations with respect of time (dN/dt).
%   X(1:4)'=X .* (gamma - V * X) - You may also use such matrix notation. The length 
%   of the vector of state variables must be set by use of a colon. *. is a product of 
%   arrays, * is a matrix product. See MATLAB manuals.
%   X(1:10,1:10)(t)= f(X(t-1)) - State variables can also be matrices (in this example a 10x10 
%   matrix).  There are special commands for viewing the matrices (<a href="matlab:help viewcells">viewcells</a>).
%   function res=monod(x,h);
%   res=x/(x+h);
%   return; - More complex functions with arguments. All code between function 
%   and return is copied into a function file. You may use any MATLAB statement here,
%   but matrix multiplication is not supported. Finish your function always with "return".
%
%
%   Default values of the parameters/Commands:
%   N=0.01;  - assign initial values and default parameters.
%   r=0.5;
%   ax x N [0 10];
%   (the semicolon is not required, but suppresses unnecessary
%   output)
%   setevent('simpleevent',0,'A=A+1',30);
%  
%
%   Usage:
%   MODEL - Opens with the current model
%   MODEL FILE - Opens the ini file named FILE
%   MODEL('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'commands' [cell] - cell array of strings with parameter initializations and default commands
%     'file' [string] - the name of an ini file
%     'model' [cell] - cell array of strings with model equations
%     'scheme' [cell] - cell array of strings with definition of the scheme for <a href="matlab:help vismod">vismod</a>.
%   MODEL('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' - clear: clears the current model from memory (and deletes currently temporary files)
%
%
%   See also use, savemodel, vismod, lag, rednoise, modelpanel
%
%   Reference page in Help browser:
%      <a href="matlab:commands('model')">commands model</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function model(varargin)
global g_grind;
%if the model name has spaces it can be read as separate arguments
if ~exist('i_use', 'file')
  addpath([grindpath filesep 'sys2']);
end
fieldnams={'file', 's', 'the name of an ini file','';...
   'commands', 'c', 'cell array of strings with parameter initializations and default commands',{};...
   'model', 'c', 'cell array of strings with model equations',{};...
   'scheme', 'c', 'cell array of strings with definition of the scheme for <a href="vismod.htm">vismod</a>.',{}}';
args=i_parseargs(fieldnams,'file(+)',{'-c'},varargin);
if isfield(args,'model')&&isfield(args,'commands')
    finishgrind;
    if isfield(args,'scheme')
        i_makemodel(args.model, args.commands, '', args.scheme);
    else
        i_makemodel(args.model, args.commands, '', {});
    end
    i_parcheck(1)
    return;
end
if isfield(args,'file')&&iscell(args.file)
    args.file=strtrim(sprintf('%s ', args.file{:}));
end
if ~isfield(args,'file')
    args.file=[];
end
if any(strcmp(args.opts,'-c'))
   finishgrind;
   if ~isempty(g_grind)&&isfield(g_grind,'statevars')
      if g_grind.statevars.vector
         s = sprintf('%s ', g_grind.pars{:}, g_grind.statevars.vectnames{:});
      else
         s = sprintf('%s ', g_grind.pars{:}, g_grind.statevars.names{:});
      end
      s = sprintf('clear global %s t g_*', s);
      evalin('base', s);
      rmpath([grindpath filesep 'sys2']);
      disp('Cleared model from memory and deleted temporary files');
   else
      if exist('i_use', 'file')==2
         rmpath([grindpath filesep 'sys2']);
      end
      if nargin==1
         warning('GRIND:model:noopenfile','No GRIND model is open');
      end
   end
   return;
end

if isempty(args.file)&&~isempty(g_grind)&&isfield(g_grind,'scheme')&&(length(g_grind.scheme)>10)
    warning('GRIND:model:modelhasscheme','The current model has a Forrester scheme, and should be edited with <a href="matlab:vismod">vismod</a>');
   % vismod;
   % return;
end
i_use(args.file, 1, 0);


