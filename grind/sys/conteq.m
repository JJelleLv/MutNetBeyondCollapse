%CONTEQ   Continue an equilibrium or bifurcation with COCO or MatCont
%   This command can do an 1D or 2D bifurcation analysis using the package
%   <a href="matlab:commands coco">COCO (Continuation Core and Toolboxes)</a> or <a href="matlab:commands matcont">MatCont</a>.
%   Select an equilibrium with <a href="matlab:help findeq">findeq</a>, then you can continue the stable or 
%   unstable equilibrium. It can (currently) detect a fold (=saddle-node), transcritical and a Hopf bifurcation. 
%   A plot of the equilibrium and the eigenvalues is made.
%   Limitation: implicit models, delay differential equations, nonautonomous, stochastic models and difference equations are currently not 
%   supported in COCO or MatCont. They can be continued in a more primitive way using the "grind engine" or with <a href="matlab:help paranal">paranal</a>.
%   Also higher dimensional models can better be analysed with <a href="matlab:help paranal">paranal</a>.
%
%   Usage:
%   CONTEQ -  - user is prompted for information
%   CONTEQ -out - change the default output in a dialog box and runs conteq 
%   (the plot definitions are shared with <a href="matlab:help paranal">paranal</a>).
%   CONTEQ -out plotno ['<param1>' 'funy1 funy2' funz1] [minx maxx] [miny maxy] [minz maxz] - see <a href="matlab:help paranal">paranal</a>
%   CONTEQ -list - List the plot definitions (shared with <a href="matlab:help paranal">paranal</a>)
%   CONTEQ APAR FROMPOINT CTYPE - continue APAR from the current value till TOVAL.
%   CONTEQ APAR parranges1 [FROMVAL TOVAL] - continue APAR in the range of FROMVAL to TOVAL
%   (the current value of the parameter should be included in the range).
%  
%       
%   CONTEQ('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:      
%     'activepars' [logical or cell] - indices of the active parameters (default: all 1)
%     'ctype' [cell or string] - type of curve for continuations (for instance EP)
%     'file' [string] - open a previously saved conteq session
%     'frompoint' [cell or string] - code of the point from which the continuation starts (for instance 'EP1') (default: empty)
%     'mindist' [number>0] - minimum distance between different special points (default: 1E-5)
%     'par1' [parameter] - first free parameter (default: '')
%     'par2' [parameter or empty] - second free parameter (default: '')
%     'par3' [parameter or empty] - third free parameter (default: '')
%     'parranges1' [number and length(number)==2] - range for the first parameter (default: [NaN NaN])
%     'parranges2' [number and length(number)==2] - range for the second parameter (default: [NaN NaN])
%     'stateranges' [number] - min/max value for each of the state variables, or one row if all are the same (default: [Nan NaN])
%     'symbolic' [logical] - if available use the symbolic toolbox for Jacobians (default 1)
%    
%     <a href="matlab:commands Coco">Specific options for COCO</a>
%     <a href="matlab:commands matcont">Specific options for MatCont</a> for differential equations
%     <a href="matlab:commands matcontm">Specific options for MatContM</a> (for difference equations)
%       
%   CONTEQ('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' - clear the previous run
%     '-coco' - use the COCO engine
%     '-defaults' - reset to default settings
%     '-matcont' - use the MatCont engine
%     '-l' - list the bifurcations
%     '-out' - set the output
%    
%
%     Labels of special points (in brackets if it can be continued in one or two dimensions-po=periodic orbit)
%     EP = equilibrium point (equilibria found with <a href="matlab:help findeqs">findeqs</a>) [1D].
%     H = supercritical Hopf bifurcation (negative Lyapunov coefficient) [1Dpo/2D].
%     H+ = subcritical Hopf bifurcation (positive Lyapunov coefficient) [1Dpo/2D].
%     HB = Hopf bifurcation (coco label).
%     T = transcritical bifurcation. [2D]
%     F = fold or saddle-node bifurcation. [2D]
%     BP = branch point (coincides usually with T). [1D]
%     EP = end point of continuation. [1D/2D]
%     EPS = end point of continuation if the range for state variables is exceeded [1D]
%     IO = initial orbit [1Dpo].
%     NS = neutral saddle.
%     NSA = netral saddle (coco label).
%     SN = saddle node or transcritical (coco label).
%  
%
%
%
%   See also findeq, paranal, open_matcont_gui
%
%   Reference page in Help browser:
%      <a href="matlab:commands('conteq')">commands conteq</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function conteq(varargin)
evalin('base','global g_cont;');
global g_grind g_cont;
if ~(nargin==1&&strncmp(varargin{1},'??',2))
    i_parcheck;
end
if isa(g_cont,'grind_matcont')&&~isfield(g_grind.cont,'matcont')
    g_grind.cont.matcont=g_cont;
    g_grind.cont.engine=g_cont.settings.derived.engine;
end
if isa(g_cont,'grind_coco')&&~isfield(g_grind.cont,'coco')
    g_grind.cont.coco=g_cont;
    g_grind.cont.engine='coco';
end  
docontbif=false;
if nargin>0&&strncmp(varargin{1}, '-o', 2)
    if nargin == 1
        i_paranaloutdlg('Edit Plots for conteq');
        conteq;
    else
        paranal(varargin{:});
    end
    return;
end
if nargin>0&&strncmp(varargin{1}, '-c', 2)&&~strncmp(varargin{1}, '-co', 3)
    if ~isempty(g_cont)
        disp('cleared conteq session')
        g_cont.clear;
    end
    return;
end

if  nargin>0&&strncmp(varargin{1}, '-l', 2)
    if ~isfield(g_grind, 'paranal')
        paranal('-default');
    end
    for i = 1:length(g_grind.paranal.plots)
        fprintf('conteq -out %d [''%s'' ''%s'' ''%s''] [%g %g] [%g %g] [%g %g]\n',i, strtrim(sprintf('%s ',g_grind.paranal.plots{i}.xaxis{:})),...
            strtrim(sprintf('%s ',g_grind.paranal.plots{i}.yaxis{:})),strtrim(sprintf('%s ',g_grind.paranal.plots{i}.zaxis{:})),...
            g_grind.paranal.plots{i}.xlim, g_grind.paranal.plots{i}.ylim, g_grind.paranal.plots{i}.zlim);
    end
    return;
end
if nargin>0
    ndx=cellfun(@ischar,varargin);
    f1=strcmp(varargin(ndx),'-contbif');
    if any(f1)
        docontbif=true;
    end
    f2=strcmp(varargin(ndx),'-coco');
    if any(f2)
        g_grind.cont.engine='coco';
    end
    f3= strcmp(varargin(ndx),'-matcont');
    if any(f3)
        if g_grind.solver.isdiffer
           g_grind.cont.engine='matcontm';
        else
           g_grind.cont.engine='matcont';
        end
    end
    varargin=varargin(~(f1|f2|f3));
end

if strncmp(g_grind.cont.engine,'matcont',7)&&~isa(g_cont,'grind_matcont')
    if isfield(g_grind.cont,'matcont')&&~isempty(g_cont)
        g_cont=g_grind.cont.matcont;
    else
        g_cont=grind_matcont;
        g_grind.cont.matcont=g_cont;
    end
end
if strcmp(g_grind.cont.engine,'coco')&&~isa(g_cont,'grind_coco')
    if isfield(g_grind.cont,'coco')&&~isempty(g_cont)
        g_cont=g_grind.cont.coco;
    else
        g_cont=grind_coco;
        g_grind.cont.coco=g_cont;
    end
end
args=struct('opt','');
if ~isempty(varargin)
   args=g_cont.set(varargin{:});
end
if isfield(args,'file')
    g_cont=g_cont.load_session(args.file);
    g_grind.cont.engine=g_cont.settings.derived.engine;
    if strcmp(g_cont.settings.derived.engine,'coco')
        g_grind.cont.coco=g_cont;
    elseif strcmp(g_cont.settings.derived.engine,'matcont')
        g_grind.cont.matcont=g_cont;
    end
end
if ~(isfield(args,'par1')||isfield(args,'par2'))||isempty(g_cont.settings.derived.freepars)
    if docontbif
        g_cont.gui('contbif');
    else
        g_cont.gui('conteq');
    end
else
    if numel(g_cont.points)<=1
        g_cont.add_points(findeqs('maxtime',10));
    end
    if isempty(g_cont.settings.derived.frompoint) %by default the current initial conditions
        g_cont.settings.derived.frompoint=g_cont.add_points;
    end
    g_cont.cont;
end

