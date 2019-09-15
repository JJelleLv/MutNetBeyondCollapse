%DEFINESPACE   Define the space of a spatial explicit model.
%   With this command you can set the size, boundary conditions and other properties of the space.
%   It works with property/value string pairs. You can use this command in the model definition, or 
%   afterwards when a model is already defined (the size can dynamically change then).
%      
%   Usage:
%   DEFINESPACE - opens dialog screens.
%   DEFINESPACE [N1,N2] - defines a space with size [N1,N2], default settings.
%   DEFINESPACE('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'bc' [periodic | neumann:0 | dirichlet:0] - sets all boundary conditions to the same value, options: 
%    periodic(default), left borders right and top borders bottom. neumann: the derivative has a fixed value (default 0), 
%    dirichlet: the border has a fixed value (default 0)
%     'bottom' [periodic | neumann:0 | dirichlet:0] - sets the bottom boundary condition, see bc.
%     'cellsize' [number>0] - size of a grid cell (not yet used)
%     'diffscheme' [center | forward | backward] - scheme of discrete partial differential equation
%    (see <a href="matlab:help Dx">Dx</a> or <a href="matlab:help Dy">Dy</a>)
%     'excludevars' [parameter+ or empty] - list of parameters that should never be resized (cell)
%     'gridsize' [integer>0 and length(integer)<=2] - size of the grid, adapts sizes of state variables and parameters
%     'left' [periodic | neumann:0 | dirichlet:0] - sets the left boundary condition, see bc.
%     'right' [periodic | neumann:0 | dirichlet:0] - sets the right boundary condition, see bc.
%     'top' [periodic | neumann:0 | dirichlet:0] - sets the top boundary condition, see bc.
%     'type' [lde | pde | ca] - type of grid, pde (pdepe, not yet implemented, lde(default)=lattice differential 
%    equation, ca = cellular automata (the only option for difference equations).
%   DEFINESPACE('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-d' - sets space to default values.
%     '-l' - Lists the current options.
%     '-listshort' - internally used for a shorter list.
%     '-u' - update the dimensions.
%  
%
%   See also Dx, Dy, Dxx, Dyy, neighborcells, model, setdimension
%
%
%   Reference page in Help browser:
%      <a href="matlab:commands('definespace')">commands definespace</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function definespace(varargin)
global g_grind;
fieldnams={'gridsize', 'i>0&length(i)<=2', 'size of the grid, adapts sizes of state variables and parameters', [1 1];...
   'type', 'e[l+de|p+de|c+a]', 'type of grid, pde (pdepe, not yet implemented, lde(default)=lattice differential','lde';...
   'bc', 'e[p+eriodic|n+eumann:0|d+irichlet:0]', 'sets all boundary conditions to the same value, options:','periodic';...
   'left', 'e[p+eriodic|n+eumann:0|d+irichlet:0]', 'sets the left boundary condition, see bc.','periodic';...
   'right', 'e[p+eriodic|n+eumann:0|d+irichlet:0]', 'sets the right boundary condition, see bc.','periodic';...
   'top', 'e[p+eriodic|n+eumann:0|d+irichlet:0]', 'sets the top boundary condition, see bc.','periodic';...
   'bottom', 'e[p+eriodic|n+eumann:0|d+irichlet:0]','sets the bottom boundary condition, see bc.','periodic';...
   'diffscheme', 'e[c+enter|f+orward|b+ackward]', 'scheme of discrete partial differential equation','center';...
   'excludevars', 'p+#E', 'list of parameters that should never be resized (cell)',[];...
   'cellsize', 'n>0', 'size of a grid cell (not yet used)',1}';
args=i_parseargs(fieldnams,'gridsize','-listshort,-u,-d,-l', varargin);
if isempty(g_grind)
    error('grind:definespace','No model defined');
end
if g_grind.solver.haslags
    error('grind:definespace','Lags are not yet supported in spatial explicit models');
end
if nargin == 0
    if ~isfield(g_grind, 'space')
        definespace('-defaults');
    end
    prompt={'Size of grid','Type of spatial model (lde/pde/ca)','Size of cell',...
        'Left boundary condition (periodic/neumann:value/dirichlet:value)','Right boundary condition (periodic/neumann:value/dirichlet:value)',...
        'Top boundary condition (periodic/neumann:value/dirichlet:value)','Bottom boundary condition (periodic/neumann:value/dirichlet:value)','Lattice scheme for single derivative (center/forward/backward)',...
        'List of vars/pars that are NOT associated to space'};
    if ~isempty(g_grind.space.excludevars)
        exvars = strtrim(sprintf('%s ', g_grind.space.excludevars{:}));
    else
        exvars = '';
    end
    defaultanswer = {mat2str(g_grind.space.gridsize), g_grind.space.type, mat2str(g_grind.space.cellsize), ...
        dispbc(1),dispbc(2),dispbc(3),dispbc(4),g_grind.space.diffscheme,exvars};
    answer = transpose(inputdlg(prompt, 'Define space for your model', 1, defaultanswer));
    res.gridsize = str2num(answer{1});  %#ok<ST2NM>
    res.type = answer{2};
    res.cellsize = str2num(answer{3});  %#ok<ST2NM>
    %  res.spacenames = answer(4:5);
    res.bc = answer(4:7);
    res.diffscheme = answer{8};
    res.excludevars = answer(9);
    definespace(res);
    return;
end
if any(strcmp(args.opts,'-d'))
    if ~isempty(g_grind)&&isfield(g_grind, 'statevars')&&g_grind.statevars.vector
        dim = [g_grind.statevars.dims{1}.dim1, g_grind.statevars.dims{1}.dim2];
    else
        dim = [1 1];
    end
    g_grind.space.gridsize = dim;
    if isfield(g_grind,'solver')&&isfield(g_grind.solver,'isdiffer')&&g_grind.solver.isdiffer
        g_grind.space.type = 'ca'; %ca, pde (only for 1D: pdepe)
    else
        g_grind.space.type = 'lde'; %ca, pde (only for 1D: pdepe)
    end
    
    g_grind.space.cellsize = [1 1];
    g_grind.space.bc = {'periodic','periodic','periodic','periodic'}; %left right top bottom ...
    g_grind.space.bcvalue=[NaN, NaN, NaN, NaN]; %left right top bottom ...
    g_grind.space.diffscheme = 'center';
    g_grind.space.excludevars = {}; %exclude parameters / variables that may not change in dimensions
    return;
end
if any(strcmp(args.opts,'-u'))  %update
    pars = setdiff(g_grind.pars, g_grind.space.excludevars);
    if g_grind.statevars.vector
        vars = setdiff(g_grind.statevars.vectnames, g_grind.space.excludevars);
    else
        vars = setdiff(g_grind.statevars.names, g_grind.space.excludevars);
    end
    if ~isempty(vars)
        if ~isempty(pars)
            setdimension(pars{:}, vars{1});
        end
        for i = 1:length(vars)
            setdimension(vars{i}, g_grind.space.gridsize);
        end
        if prod(g_grind.space.gridsize) == 1
            setdimension('-vector','off');
        end
    end
    return;
end
if any(strcmp(args.opts,'-listshort'))
    if ~isfield(g_grind,'space')
        return;
    end
    fprintf('<a href="matlab:definespace(''-l'')">Space</a> (%s) defined, size = <%dx%d>\n', g_grind.space.type,g_grind.space.gridsize);
    if all(strcmp(g_grind.space.bc{1}, g_grind.space.bc))&&(all(g_grind.space.bcvalue(1)==g_grind.space.bcvalue)||all(isnan(g_grind.space.bcvalue)))
        fprintf('All boundary conditions = ''%s''\n',dispbc(1));
    else
        fprintf('Boundary conditions = (''%s'',''%s'',''%s'',''%s'')\n',dispbc(1),dispbc(2),dispbc(3),dispbc(4));
    end
    if strcmp(g_grind.space.type, 'lde')
        fprintf('Scheme for single differentiation (see Dx or Dy) = ''%s''\n', g_grind.space.diffscheme);
    end
    return
elseif any(strcmp(args.opts,'-l'))
    if ~isfield(g_grind, 'space')||isempty(g_grind.space)
        disp('No space defined');
    else
        fprintf('Type of space [type]:                   ''%s''\n',g_grind.space.type);
        fprintf('Grid size [gridsize]:                   <%dx%d>\n', g_grind.space.gridsize);
        fprintf('Size of gridcell [cellsize]:            <%gx%g>\n', g_grind.space.cellsize);
        %           fprintf('Name of x direction in space [xname]:   ''%s''\n',g_grind.space.spacenames{1});
        %           fprintf('Name of y direction in space [yname]:   ''%s''\n',g_grind.space.spacenames{2});
        if all(strcmp(g_grind.space.bc{1}, g_grind.space.bc))&&(all(g_grind.space.bcvalue(1)==g_grind.space.bcvalue)||all(isnan(g_grind.space.bcvalue)))
            fprintf('Boundary conditions [bc]:               ''%s''\n',dispbc(1));
        else
            fprintf('Left boundary condition [left]:         ''%s''\n',dispbc(1));
            fprintf('Right boundary condition [right]:       ''%s''\n',dispbc(2));
            fprintf('Top boundary conditions [top]:          ''%s''\n',dispbc(3));
            fprintf('Bottom boundary conditions [bottom]:    ''%s''\n',dispbc(4));
        end
        if strcmp(g_grind.space.type, 'lde')
            fprintf('Scheme for [du/dx] [diffscheme]:       ''%s''\n', g_grind.space.diffscheme);
        end
        fprintf('Pars not related to space[excludevars]: {%s}\n',sprintf('%s ',g_grind.space.excludevars{:}));
    end
    return;
end

if ~isfield(g_grind, 'space')
    definespace(struct('opts',{'-d'}));
end
oldspace=g_grind.space;
s = rmfield(args,'opts');
f1 = fieldnames(s);
for i = 1:length(f1)
    if strcmpi(f1{i}, 'bc') %shortcut if all the same boundary conditions
        if ischar(s.bc)
            s.bc={s.bc,s.bc,s.bc,s.bc};
        end
        g_grind.space.bcvalue(1:length(s.bc))=getbcvalue(s.bc);
        g_grind.space.bc(1:length(s.bc)) = getbc(s.bc);
        if length(g_grind.space.bcvalue)==1
            g_grind.space.bcvalue=repmat(g_grind.space.bcvalue,1,4);
        end
    elseif strcmpi(f1{i}, 'diffscheme')
        if strncmpi(s.(f1{i}), 'c', 1)
            g_grind.space.(f1{i}) = 'center';
        elseif strncmpi(s.(f1{i}), 'f', 1)
            g_grind.space.(f1{i}) = 'forward';
        elseif strncmpi(s.(f1{i}), 'b', 1)
            g_grind.space.(f1{i}) = 'backward';
        else
            error('grind:definespace','Unknown scheme for single differentiation');
        end
    elseif strcmpi(f1{i},'left')
        g_grind.space.bcvalue(1)=getbcvalue(s.(f1{i}));
        g_grind.space.bc{1} = getbc(s.(f1{i}));
    elseif strcmpi(f1{i},'right')
        g_grind.space.bcvalue(2)=getbcvalue(s.(f1{i}));
        g_grind.space.bc{2} = getbc(s.(f1{i}));
    elseif strcmpi(f1{i},'top')
        g_grind.space.bcvalue(3)=getbcvalue(s.(f1{i}));
        g_grind.space.bc{3} = getbc(s.(f1{i}));
    elseif strcmpi(f1{i},'bottom')
        g_grind.space.bcvalue(4)=getbcvalue(s.(f1{i}));
        g_grind.space.bc{4} = getbc(s.(f1{i}));
    elseif any(strcmpi(f1{i}, {'type'}))
        if any(strcmpi(s.(f1{i}),{'pde','lde','ca'}))
            g_grind.space.(f1{i}) = s.(f1{i});
        else
            error('grind:definespace','Unknown type for space ("%s")',s.(f1{i}));
        end
    elseif any(strcmpi(f1{i},{'cellsize','gridsize'}))
        g_grind.space.(f1{i}) = i_checkstr(s.(f1{i}));
        if numel(g_grind.space.(f1{i})) == 1
            g_grind.space.(f1{i}) = [g_grind.space.(f1{i}), g_grind.space.(f1{i})];
        end
        if strcmpi(f1{i}, 'gridsize')
            definespace(struct('opts',{'-u'}));
        end
    else
        error('grind:definespace','Unknown property "%s" for space',f1{i});
    end
end
if ~struccmp(g_grind.space,oldspace)
    disp('Space settings changed to:');
    definespace(struct('opts',{'-l'}))
    g_grind.checks.lastsettings=[];
end

function v=getbcvalue(s)
if iscell(s)
    v=zeros(size(s));
    for i=1:length(s)
        v(i)=getbcvalue(s{i});
    end
    return;
end
if strncmpi(s, 'p', 1)
    v=nan;
    return;
else
    v=0;
end
f=strfind(s,':');
if ~isempty(f)
    v=str2double(s(f+1:end));
end
function s1=getbc(s)
if iscell(s)
    s1=cell(size(s));
    for i=1:length(s)
        s1{i}=getbc(s{i});
    end
    return;
end
if strncmpi(s, 'p', 1)
    s1 = 'periodic';
elseif strncmpi(s, 'n', 1)
    s1 = 'neumann';
elseif strncmpi(s, 'd', 1)||strncmpi(s, 'f', 1) %fixed
    s1 = 'dirichlet';
else
    error('grind:definespace','Boundary condition ''%s'' not recognized',s);
end
function s=dispbc(i)
global g_grind;
s=g_grind.space.bc{i};
if ~strcmp(g_grind.space.bc{i},'periodic')
    s=sprintf('%s:%g',s,g_grind.space.bcvalue(i));
end


