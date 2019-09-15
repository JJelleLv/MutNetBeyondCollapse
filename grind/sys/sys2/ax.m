%AX   Define axes of phase plane (variable and range) 
%   Define an axis for the phase plane.
%
%   Usage:
%   AX - show current axis settings.
%   AX AXNAME - clear the axis AXNAME (AXNAME can be X, Y or Z). 
%   AX AXNAME [] [LIM(1) LIM(2)] - set the limits of axis AXNAME.
%   AX AXNAME VAR [LIM(1) LIM(2)] - set the function/variable of
%   AXNAME and set the limits. VAR1 may be a state variable, a 
%   parameter, an auxiliary variable (see <a href="matlab:help funcs">funcs</a>), or an MATLAB 
%   expression with these variables.
%   AX AXNAME1<>AXNAME2 - exchange two axes.
%   AX('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'axname' [x | y | z] - the name of the axis (can be X, Y, or Z)
%     'lim' [number and length(number)<=2] - limits or minimum limit of the axis.
%     'limmax' [number] - maximum limit of the axis.
%     'var' [string] - state variable or function on axis.
%   AX('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-def' - default axis.
%     '-?' - list the current axes.
%        
%
%
%   Examples: 
%   ax x A  - sets the state variable 'A' on the x axis.
%   ax y K [0 10] - sets the parameter 'K' on the y axis with a
%   range of 0-10
%   ax z K*A - sets the expression 'K*A' on the z axis.
%   ax z [] [0 100] - sets the range of the z axis to 0-100
%
%   ax z - clears the z axis.
%   ax x<>y - exchanges x and y axis.
%
%   See also null, null3, phas, ru, paranal, outf
%
%   Reference page in Help browser:
%      <a href="matlab:commands('ax')">commands ax</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
function ax(varargin)
%axname, axvar, axlim, axmax)
global g_grind g_paranal;
fieldnams={'axname', 'e[x|y|z]', 'the name of the axis (can be X, Y, or Z)','x';...
   'var', 's', 'state variable or function on axis.',g_grind.xaxis.var;...
   'lim', 'n&length(n)<=2', 'limits or minimum limit of the axis.',[0 10];...
   'limmax', 'n', 'maximum limit of the axis.',[]}';
args=i_parseargs(fieldnams,'if ~(argtype(1,''e[x|y|z]'')),deffields=''var,lim,limmax'';elseif(argtype(2,''d'')),deffields=''axname,lim,limmax'';else,deffields=''axname,var,lim,limmax'';end;','-def,-?',varargin);
if isfield(args,'var')
    if ~isempty(strfind(args.var,'['))
        args.lim=i_checkstr(args.var);
        args.var='';
    elseif ~isempty(strfind(args.var,'<>'))
        args.axname=args.var;
        args =rmfield(args,'var');
    end
end
if isfield(args,'limmax')
    args.lim=[args.lim args.limmax];
end
if ~isfield(args,'var')&&~isfield(args,'lim')
    args.var=' ';
end

if any(strcmp(args.opts,'-def'))
    irhs = g_grind.statevars.dim;
    g_grind.xaxis.var ='';
    g_grind.yaxis.var ='';
    g_grind.zaxis.var ='';
    if (irhs > 0)
        g_grind.xaxis.var = char(i_statevars_names(1));
    end
    if (irhs > 1)
        g_grind.yaxis.var = char(i_statevars_names(2));
    end
    if (irhs > 2)
        g_grind.zaxis.var = char(i_statevars_names(3));
    end
    if nargin == 1
        ax('-?');
    end
    return;
end
% if strcmp(args.var, '[]')
%    args.var = [];
% end
if isfield(g_paranal,'nulldata')
    g_paranal.nulldata=[];
end
if ~isfield(args,'lim')
    args.lim=[];
end
if isfield(args,'var')&&~isfield(args,'axname')
    if strcmp(g_grind.xaxis.var,args.var)
        args.axname='x';
    elseif strcmp(g_grind.yaxis.var,args.var)
        args.axname='y';
    elseif strcmp(g_grind.zaxis.var,args.var)
        args.axname='z';
    end
end

%    f1=strfind(args.var,'[');
%    f2=strfind(args.var,']');
%    if ~isempty(f1)&&~isempty(f2)
%       args.lim = i_checkstr(args.var);
%       args.var = [];
%    end
% end
if any(strcmp(args.opts,'-?'))||~isfield(args,'axname') || strcmp(args.axname, '?')
    
    if ~isfield(g_grind, 'xaxis')||isempty(g_grind.xaxis) || isempty(g_grind.yaxis) || isempty(g_grind.zaxis)
        disp('currently no axis defined');
    else
        if nargin == 0
            i_axdlg
        end
        disp(' ');
        disp ('Current axis settings:');
        printax('x', g_grind.xaxis.var , g_grind.xaxis.lim);
        printax('y', g_grind.yaxis.var , g_grind.yaxis.lim);
        printax('z', g_grind.zaxis.var , g_grind.zaxis.lim);
    end
else
    switch args.axname
        case '1'
            args.axname = 'x';
        case '2'
            args.axname = 'y';
        case '3'
            args.axname = 'z';
    end
    if isfield(args,'var')&&~isempty(args.var)&&~isnan(str2double(args.var))
        error('GRIND:ax:VarUnknown','Axis "%s" cannot be "%s"',args.axname,args.var);
    end
    if strcmpi(args.axname, 'x')
        if isfield(args,'var')&&~isempty(args.var)
            g_grind.xaxis.var = strtrim(args.var);
        end
        g_grind.xaxis.lim = i_checklim(g_grind.xaxis.lim,args.lim);
    elseif strcmpi(args.axname, 'y')
        if isfield(args,'var')&&~isempty(args.var)
            g_grind.yaxis.var = strtrim(args.var);
        end
        g_grind.yaxis.lim = i_checklim(g_grind.yaxis.lim,args.lim);
    elseif strcmpi(args.axname, 'z')
        if isfield(args,'var')&&~isempty(args.var)
            g_grind.zaxis.var = strtrim(args.var);
        end
        g_grind.zaxis.lim = i_checklim(g_grind.zaxis.lim,args.lim);
    elseif strcmpi(args.axname,'x<>y')||strcmpi(args.axname,'y<>x')
        hlp = g_grind.yaxis;
        g_grind.yaxis = g_grind.xaxis;
        g_grind.xaxis = hlp;
    elseif strcmpi(args.axname,'x<>z')||strcmpi(args.axname,'z<>x')
        hlp = g_grind.zaxis;
        g_grind.zaxis = g_grind.xaxis;
        g_grind.xaxis = hlp;
    elseif strcmpi(args.axname,'y<>z')||strcmpi(args.axname,'z<>y')
        hlp = g_grind.zaxis;
        g_grind.zaxis = g_grind.yaxis;
        g_grind.yaxis = hlp;
    else
        if isempty(args.lim)
            args.lim = str2num(args.var); %#ok<ST2NM>  not str2double
        end
        if (~isempty(args.lim)) && (length(args.lim) > 1)
            if strcmp(args.axname, g_grind.xaxis.var)
                g_grind.xaxis.lim = i_checklim(g_grind.xaxis.lim,args.lim);
            elseif strcmp(args.axname, g_grind.yaxis.var)
                g_grind.yaxis.lim = i_checklim(g_grind.yaxis.lim,args.lim);
            elseif strcmp(args.axname, g_grind.zaxis.var)
                g_grind.zaxis.lim = i_checklim(g_grind.zaxis.lim,args.lim);
            else
                error('GRIND:ax:VarUnknown','Axis "%s" does not exist',args.axname);
            end
        end
    end
end
if isfield(g_grind,'xaxis')
    g_grind.xaxis.var = outf('changeshortcut', g_grind.xaxis.var);
    g_grind.yaxis.var = outf('changeshortcut', g_grind.yaxis.var);
    g_grind.zaxis.var = outf('changeshortcut', g_grind.zaxis.var);
end
function printax(axname, var, lim)
varno = i_getno(var,true);
var = outf('changeshortcut', var);
if ~isempty(varno.no)
    fprintf('ax %s %s [%g %g]\n',char(axname), char(var) ,lim)
else
    fprintf('ax %s "%s" is not a valid axis\n',char(axname), char(var));
end




