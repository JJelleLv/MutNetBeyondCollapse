%MAKEMAP  Discretize a differential equation model
%   A differential equation model is changed to a discrete difference equation. THere are different
%   kinds of difference equations possible: (1) "period" = fixed period map, (2) "minima" =
%   difference equation of subsequent minima (3) "maxima" (Lorenz map) -
%   difference equation of subsequent maxima. (4) "poincare" (Poincare map) -
%   difference equation of crossings through a plane in the state space. After this you 
%   can use most normal tools like <a href="matlab:help time">time</a>, <a href="matlab:help itermap">itermap</a>, <a href="matlab:help paranal">paranal</a> etc.
%   (for minima/maxima itermap does not work properly)
%
%   Usage:
%   MAKEMAP - Select the map and parameters in a dialog box.
%   MAKEMAP 'period' TPERIOD PHASE - make a map with a fixed period of TPERIOD and a sampling 
%   time of PHASE.
%   MAKEMAP 'period' TPERIOD PHASE - make a map with a fixed period of TPERIOD and a sampling 
%   time of PHASE.
%   MAKEMAP 'minima' VAR - Map of minima of variable VAR. 
%   MAKEMAP 'maxima' VAR - Map of maxima of variable VAR. 
%   MAKEMAP 'poincare' PLANE NPOINTS - Poincare map. PLANE is the equation that defines a plane (e.g. "z=10"), 
%   NPOINTS is the number of point the model is run to find the plane (default see <a href="matlab:help simtime">simtime</a>). 
%   MAKEMAP('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'npoints' [number>0] - the period in which the model is run to find the plane of Poincare map
%     'phase' [string] - a sampling time of PHASE
%     'plane' [string] - equation defining a plane for a Poincare map
%     'tperiod' [string] - fixed period of the period map
%     'type' [normal | poincare | period | minima | maxima] - string with type of map
%     'var' [state variable] - state variable to use
%   MAKEMAP('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' - clear the map and return to the original model.
%     '-l' - list the current map.
%
%
%   See also makecartesian,lorenzmap, poincaremap, itermap  
%
%   Reference page in Help browser:
%      <a href="matlab:commands('makemap')">commands makemap</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function makemap(varargin)
%(maptype, p1, p2)
global g_grind;
if nargin == 0
    mtypes={'normal simulation','period', 'minima', 'maxima', 'poincare'};
    iinit=1;
    defper='100';
    defper2=g_grind.ndays;
    defsamp='0';
    defplane='';
    defvar=i_statevars_names(1);
    if ~isfield(g_grind.solver, 'map')
        g_grind.solver.map = [];
    elseif ~isempty(g_grind.solver.map)
        if strcmp(g_grind.solver.map.type,'period')
            defper=g_grind.solver.map.tperiod;
            defsamp=g_grind.solver.map.tsample;
        elseif  strcmp(g_grind.solver.map.type,'poincare')
            defper2=g_grind.solver.map.tperiod;
            defplane=g_grind.solver.map.plane;
        else
            defvar=i_statevars_names(g_grind.solver.map.varno);
        end
    end
    [args.type,isok] = listdlg('PromptString','Select a map:',...
        'SelectionMode','single',...
        'InitialValue',iinit,...
        'ListString', mtypes);
    if isok
        args.type = mtypes{args.type};
        switch args.type
            case 'normal simulation'
                makemap('-c');
            case 'period'
                prompt={'Enter the period (may be an expression/parameter)','Enter the sampling time (may be an expression/parameter)'};
                answer = inputdlg(prompt, 'Input for make map', 1, {defper,defsamp});
                if ~isempty(answer)
                    makemap(args.type, answer{1}, answer{2});
                end
            case 'poincare'
                prompt={'Enter an equation for the plane','Enter the period for running per step'};
                answer = inputdlg(prompt, 'Input for make map', 1, {defplane,int2str(defper2)});
                if ~isempty(answer)
                    makemap(args.type, answer{1}, answer{2});
                end
            otherwise %minima or maxima
                prompt = {sprintf('Enter the variable for %s', args.type)};
                answer = inputdlg(prompt, 'Input for make map', 1, {defvar});
                if ~isempty(answer)
                    makemap(args.type, answer{1});
                end
        end
    end
    return;
end
fieldnams={'type', 'e[n+ormal|po+incare|p+eriod|mi+nima|m+axima]', 'string with type of map','poincare';...
   'tperiod', 's', 'fixed period of the period map','';...
   'var', 'v', 'state variable to use',i_statevars_names(1);...
   'plane', 's', 'equation defining a plane for a Poincare map','';...
   'phase', 's', 'a sampling time of PHASE','0';...
   'npoints', 'n>0', 'the period in which the model is run to find the plane of Poincare map',1000}';
args=i_parseargs(fieldnams,['if(strncmpi(args{1},''po'',2)),deffields=''type,plane,npoints'';',...
    'elseif(strncmpi(args{1},''p'',1)),deffields=''type,tperiod,phase'';',...
    'elseif(strncmpi(args{1},''m'',1)),deffields=''type,var'';else,deffields=''type'';end;'],'-c,-l',varargin);

if any(strcmp(args.opts, '-c'))||(isfield(args,'type')&&strcmp(args.type,'normal'))
    domakemap('');
    makemap('-l');
    return;
elseif any(strcmp(args.opts, '-l'))
    if isfield(g_grind.solver,'map')&&~isempty(g_grind.solver.map)
        switch g_grind.solver.map.type
            case 'poincare'
                fprintf('Map type: %s\nEquation of plane: %s==0\nNumber of points: %g\n',g_grind.solver.map.type,g_grind.solver.map.plane,g_grind.solver.map.tperiod);
            case 'period'
                fprintf('Map type: %s\nPeriod: %s\nPhase: %s\n',g_grind.solver.map.type,g_grind.solver.map.tperiod,g_grind.solver.map.tsample);
            case 'maxima'
                fprintf('Map type: %s\nVariable: %s\n',g_grind.solver.map.type,i_statevars_names(g_grind.solver.map.varno));
            case 'minima'
                fprintf('Map type: %s\nVariable: %s\n',g_grind.solver.map.type,i_statevars_names(g_grind.solver.map.varno));
            otherwise
                disp('Normal simulation');
        end
    else
        disp('Normal simulation');
    end
    return;
end
if ~isfield(args,'type')
    args.type='';
end
if strncmpi(args.type, 'po', 2)
    domakemap('poincare')
    g_grind.checks.lastsettings = [];
    if ~isfield(args,'npoints')
        g_grind.solver.map.tperiod = 1000;
    else
        g_grind.solver.map.tperiod = args.npoints;
    end
    if ~isfield(args,'plane')
        prompt={'Enter an equation for the plane'};
        answer = inputdlg(prompt, 'Input for make map', 1, {''});
        if ~isempty(answer)
            args.plane=answer{1};
        end
    end
    if args.plane(end)==';'
        args.plane=args.plane(1:end-1);
    end
    fs=parsed_equation(args.plane).fs;
    eqs={'==','=>','=<','=','<','>'};
    for i=1:length(eqs)
        f=strcmp(fs,eqs{i});
        if any(f)
            f=find(f,1);
            if f==length(fs)-1&&strcmp(fs{end},'0')
                fs=fs(1:end-2);
            else
                fs=[fs(1:f-1) {'-','('} fs(f+1:end) {')'}];
            end
        end
    end
    p1=minbrackets(parsed_equation(sprintf('%s',fs{:})));
    g_grind.solver.map.plane=char(p1);
elseif strncmpi(args.type, 'p', 1)
    domakemap('period')
    g_grind.checks.lastsettings = [];
    if ~isfield(args,'phase')
        g_grind.solver.map.tsample = '0';
    else
        g_grind.solver.map.tsample = args.phase;
    end
    if ~isfield(args,'tperiod')
        prompt = {'Enter the period for the map'};
        defvar='100';
        answer = inputdlg(prompt, 'Input for make map', 1, {defvar});
        if ~isempty(answer)
            g_grind.solver.map.tperiod=answer{1};
        else
            error('makemap:grind','missing period')
        end
    else
        g_grind.solver.map.tperiod = args.tperiod;
    end
elseif strncmpi(args.type, 'm', 1)
    if g_grind.solver.nonautonomous
        error('grind:makemap:minima','A map of minima/maxima is only possible for autonomous ODE''s:\nYou can make your model suitable by adding tau''=1 and use tau in equations instead of t');
    end
    if strncmpi(args.type, 'mi', 2)
        domakemap('minima')
        g_grind.solver.map.dir = 1;
    else
        domakemap('maxima')
        g_grind.solver.map.dir = -1;
    end
    if ~isfield(args,'var')
        g_grind.solver.map.varno = 1;
    else
        ix = i_getno(args.var);
        g_grind.solver.map.varno = ix.no;
    end
    g_grind.checks.lastsettings = [];
else
    error('grind:makemap','type of map not recognized');
end
makemap('-l')

function domakemap(aname)
global g_grind;
switch aname
    case 'period'
        odename = 'i_period_map';
    case 'poincare'
        odename = 'i_poincare_map';
    case 'maxima'
        odename = 'i_minmax_map';
    case 'minima'
        odename = 'i_minmax_map';
    otherwise
        odename = '';
end
if ~isfield(g_grind.solver, 'map')
    g_grind.solver.map = [];
end
if isempty(aname)
    g_grind.checks.lastsettings = [];
    if ~isempty(g_grind.solver.map)
        if isfield(g_grind.solver.map, 'solver')
            g_grind.solver.name = g_grind.solver.map.solver;
        end
        if isfield(g_grind.solver.map, 'opt')
            g_grind.solver.opt = g_grind.solver.map.opt;
        end
        if isfield(g_grind.solver.map, 'odefile')
            g_grind.odefile = g_grind.solver.map.odefile;
        end
        if isfield(g_grind.solver.map, 'haslags')
            g_grind.solver.haslags = g_grind.solver.map.haslags;
        end
    end
    g_grind.solver.isdiffer = 0;
    g_grind.solver.map = [];
elseif isempty(g_grind.solver.map)||~isfield(g_grind.solver.map,'type')||~strcmp(aname, g_grind.solver.map.type)
    if ~isempty(g_grind.solver.map)
        domakemap('');
    end
    g_grind.solver.map.type = aname;
    g_grind.solver.map.odefile = g_grind.odefile;
    g_grind.solver.map.solver = g_grind.solver.name;
    g_grind.solver.map.opt = g_grind.solver.opt;
    g_grind.solver.map.haslags=g_grind.solver.haslags;
    g_grind.odefile = odename;
    g_grind.solver.isdiffer = 1;
    g_grind.solver.haslags =false;
    g_grind.solver.name = 'i_differ';
end


