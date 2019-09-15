%MAKECARTESIAN  Change the coordinates of a model in polar coordinates
%   A differential equation can be defined in polar coordinates. This command can
%   automatically change the coordinates to normal cartesian coordinates.
%   Example of a model in polar coordinates:
%   r'=-a*r
%   theta'=b
%   Where r is the radius and theta is the angle. With makecartesian this model is changed
%   to cartesian coordinates.
%
% 
%   Usage:
%   MAKECARTESIAN - make cartesian coordinates with default settings. (RADIUSVAR = r (or the first
%   variable, THETAVAR = thetat (or the second variable))
%   MAKECARTESIAN RADIUSVAR THETAVAR - make cartesian coordinate of RADIUSVAR and THETAVAR.
%   MAKECARTESIAN -C - Clear cartesian coordinates and revert to polar coordinates.
%      
%
%   MAKECARTESIAN('argname',argvalue,...) - Valid argument name-value pairs [with type]:
%     'analytic' [logical] - analytic means that the whole odefile is replaced. This is faster
%  than the simpler non-analytic method (default true).
%     'radiusvar' [identifier] - the name of the radius variable (default 'r')
%     'thetavar' [identifier] - the name of the angle variable (default 'theta')
%     'xvar' [identifier] - the name of x variable in cartesian coordinates (default 'x')
%     'yvar' [identifier] - the name of y variable in cartesian coordinates (default 'y')
%   MAKECARTESIAN('-opt1','-opt2',...) - Valid command line options:
%     '-?' - show if the coordinates are cartesian or polar
%     '-c' - go back to polar coordinates.
%     '-s' - silent mode, suppresses output in command window
%     
%See also makemap
%
%   Reference page in Help browser:
%      <a href="matlab:commands('makecartesian')">commands makecartesian</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function makecartesian(varargin)
global g_grind
fieldnams={'radiusvar', 'U1', 'the name of the radius variable (default ''r'')','r';...
   'thetavar', 'U1', 'the name of the angle variable (default ''theta'')','theta';...
   'xvar', 'U1', 'the name of x variable in cartesian coordinates (default ''x'')','x';...
   'yvar', 'U1', 'the name of y variable in cartesian coordinates (default ''y'')','y';...
   'analytic', 'l', 'analytic means that the whole odefile is replaced. This is faster',true}';
args=i_parseargs(fieldnams,'radiusvar,thetavar', '-c,-?,-s',varargin,false,{@i_isid});
if g_grind.statevars.dim<2
    error('grind:makecartesian:dim','The model has too few statevariables for this option (at least 2 required)')
end
if g_grind.statevars.vector
    error('grind:makecartesian:dim','Vector notation is not supported for makecartesian')
end
if (nargin==1) &&any(strcmp(args.opts, '-?'))
    if isfield(g_grind.solver,'polar')
        fprintf('Polar coordinates (radius=%s, angle=%s) replaced by cartesian (%s, %s)\n',...
            g_grind.solver.polar.polarnames{g_grind.solver.polar.radius},...
            g_grind.solver.polar.polarnames{g_grind.solver.polar.theta},...
            g_grind.solver.polar.cartnames{g_grind.solver.polar.radius},...
            g_grind.solver.polar.cartnames{g_grind.solver.polar.theta});
    else
        disp('Unchanged coordinates')
    end
    return
end

if ~isfield(args,'analytic')
    args.analytic=true;
end
silent=any(strcmp(args.opts, '-s'));
if isfield(g_grind.solver,'polar')
    if ~g_grind.solver.polar.analytic
        g_grind.odefile=g_grind.solver.polar.odefile;
        g_grind.statevars.names=g_grind.solver.polar.polarnames;
    end
elseif ~any(strcmp(args.opts, '-c'))&& any(strcmp('%model adapted by grind',g_grind.model))
    if ~silent
        disp('Already cartesian coordinates')
    end
    return;
end
if ~isfield(args,'radiusvar')
    if any(strcmp(g_grind.statevars.names,'r'))&&any(strcmp(g_grind.statevars.names,'theta'))
        args.radiusvar='r';
    else
        args.radiusvar=g_grind.statevars.names{1};
    end
end
if ~isfield(args,'thetavar')
    if any(strcmp(g_grind.statevars.names,'r'))&&any(strcmp(g_grind.statevars.names,'theta'))
        args.thetavar='theta';
    else
        args.thetavar=g_grind.statevars.names{2};
    end
end
if ~isfield(args,'xvar')
    args.xvar='x';
end
if ~isfield(args,'yvar')
    args.yvar='y';
end
if any(strcmp(args.xvar,g_grind.statevars.names))
    args.xvar=[args.xvar '1'];
end
if any(strcmp(args.yvar,g_grind.statevars.names))
    args.yvar=[args.yvar '1'];
end
if strcmp(args.xvar,args.yvar)
    error('grind:makecartesion','xvar and yvar should be different names')
end
if strcmp(args.radiusvar,args.thetavar)
    error('grind:makecartesion','radiusvar and thetavar should be different state variables')
end
if any(strcmp(args.opts, '-c'))
    if isfield(g_grind.solver,'polar')||any(strcmp('%model adapted by grind',g_grind.model))
        if isfield(g_grind.solver,'polar')&&~g_grind.solver.polar.analytic
            args.xvar=g_grind.solver.polar.cartnames{g_grind.solver.polar.radius};
            args.yvar=g_grind.solver.polar.cartnames{g_grind.solver.polar.theta};
            args.radiusvar=g_grind.solver.polar.polarnames{g_grind.solver.polar.radius};
            args.thetavar=g_grind.solver.polar.polarnames{g_grind.solver.polar.theta};
            if ~silent
                fprintf('Cartesian coordinates (%s, %s) replaced by polar (radius=%s,angle=%s)\n',args.xvar,args.yvar,args.radiusvar,args.thetavar)
            end
            evalin('base',sprintf('%s=sqrt(%s^2+%s^2);',args.radiusvar,args.xvar,args.yvar));
            evalin('base',sprintf('%s=atan2(%s,%s);',args.thetavar,args.yvar,args.xvar));
            evalin('base',sprintf('clear(''%s'',''%s'');',args.xvar,args.yvar))
        else
            if ~isfield(g_grind.solver,'polar')||~isfield(g_grind.solver.polar,'polarmodel')
                polarmodel=makepolar(g_grind.model);
            else
                polarmodel=g_grind.solver.polar.polarmodel;
            end
            finishgrind;
            i_makemodel(polarmodel,g_grind.commands,g_grind.inifile)
        end
        %keep the same initial conditions
        if isfield(g_grind.solver,'polar')
            g_grind.solver=rmfield(g_grind.solver,'polar');
            g_grind.checks.lastsettings=[];
            out('-defaults','-silent')
            ax('-defaults','-silent')
        end
    elseif ~silent
        disp('Already polar coordinates')
    end
else
    %
    wascart=isfield(g_grind.solver,'polar');
    if wascart&&~args.analytic
        oldids=[g_grind.solver.polar.radius,g_grind.solver.polar.theta];
        newids=[find(strcmp(g_grind.statevars.names,args.radiusvar)), find(strcmp(g_grind.statevars.names,args.thetavar))];
        if any(oldids~=newids)
            makecartesian('-clear')
            wascart=false;
        end
    end
    g_grind.solver.polar.polarnames=g_grind.statevars.names;
    g_grind.solver.polar.analytic=args.analytic;
    g_grind.solver.polar.radius=find(strcmp(g_grind.statevars.names,args.radiusvar));
    g_grind.solver.polar.theta=find(strcmp(g_grind.statevars.names,args.thetavar));
    
    if args.analytic
        g_grind.solver.polar.polarmodel=g_grind.model;
        cartmodel=makecartmod(g_grind.model,args.radiusvar,args.thetavar,args.xvar,args.yvar);
        polar=g_grind.solver.polar;
        finishgrind
        i_makemodel(cartmodel,g_grind.commands,g_grind.inifile)
        g_grind.solver.polar=polar;
    else
        g_grind.solver.polar.odehandle=str2func(g_grind.odefile);
        g_grind.solver.polar.odefile=g_grind.odefile;
        %       g_grind.solver.polar.polarhandle=i_getodehandle('polar2cart');
        g_grind.odefile='g_grind.solver.polar.polarhandle';
        g_grind.odefile='i_polar2cart';
        g_grind.statevars.names(g_grind.solver.polar.radius)={args.xvar};
        g_grind.statevars.names(g_grind.solver.polar.theta)={args.yvar};
    end
    %keep the same initial conditions
    if ~wascart
        evalin('base',sprintf('%s=%s*cos(%s);',args.xvar,args.radiusvar,args.thetavar));
        evalin('base',sprintf('%s=%s*sin(%s);',args.yvar,args.radiusvar,args.thetavar));
        evalin('base',sprintf('clear(''%s'',''%s'');',args.radiusvar,args.thetavar))
        g_grind.solver.polar.cartnames=g_grind.statevars.names;
        if ~silent
            makecartesian('-?')
        end
    else
        if any(~strcmp(g_grind.solver.polar.cartnames,g_grind.statevars.names))
            oldx=g_grind.solver.polar.cartnames{g_grind.solver.polar.radius};
            oldy=g_grind.solver.polar.cartnames{g_grind.solver.polar.theta};
            if ~strcmp(oldx,args.xvar)
                evalin('base',sprintf('%s=%s;clear(''%s'')',args.xvar,oldx,oldx));
            end
            if ~strcmp(oldy,args.yvar)
                evalin('base',sprintf('%s=%s;clear(''%s'')',args.yvar,oldy,oldy));
            end
            if ~silent
                fprintf('Updating cartesian coordinates\nNames of state variables changed:\n%s -> %s\n%s -> %s\n',oldx,args.xvar,oldy,args.yvar)
            end
        elseif ~silent
            fprintf('Already cartesian coordinates\n')
        end
    end
    g_grind.solver.polar.cartnames=g_grind.statevars.names;
    g_grind.checks.lastsettings=[];
    out -defaults -silent;
    ax -defaults -silent
end

function mod1=makepolar(mod)
mod1=mod;
ndx=true(size(mod));
if any(strcmp('%model adapted by grind',mod1))
    f=find(strcmp('%model adapted by grind',mod1));
    ndx(f)=false;
    f=f+1;
    while f<length(mod1)&&isempty(strfind(mod1{f},'=sqrt('))
        f=f+1;
    end
    ndx(f)=false;
    f1=strfind(mod1{f},'=');
    radius=mod1{f}(1:f1-1);
    %         vars= regexp(mod1(f),'(?<=[\(\+])[a-zA-Z0-9_]*','match');
    %         xvar=vars{1}{1};
    %         yvar=vars{1}{2};
    while f<length(mod1)&&isempty(strfind(mod1{f},'=atan2('))
        f=f+1;
    end
    ndx(f)=false;
    f1=strfind(mod1{f},'=');
    theta=mod1{f}(1:f1-1);
end
dr_dt_=sprintf('d%s_dt_',radius);
f=find(strncmp(mod1,dr_dt_,length(dr_dt_)));
if ~isempty(f)
    s=mod1{f};
    mod1{f}=sprintf('%s''%s',radius,s(length(dr_dt_)+1:end));
end
dtheta_dt_=sprintf('d%s_dt_',theta);
f=find(strncmp(mod1,dtheta_dt_,length(dtheta_dt_)));
if ~isempty(f)
    s=mod1{f};
    mod1{f}=sprintf('%s''%s',theta,s(length(dtheta_dt_)+1:end));
end
f=strcmp('%original polar model:',mod1);
if any(f)
    ndx(f)=false;
end
f=strcmp('%new cartesian model:',mod1);
if any(f)
    ndx(f)=false;
end
f=strfind(mod1,dtheta_dt_);
f=~cellfun('isempty',f);
ndx(f)=false;
mod1=mod1(ndx);



function mod1=makecartmod(mod,radius,theta,xvar,yvar)
mod1=mod;
if ~any(strcmp('%model adapted by grind',mod))
    i1=strncmp([radius ''''],mod1,length(radius)+1);
    if any(i1)
        dr_dt_=sprintf('d%s_dt_',radius);
        mod1{i1}=[dr_dt_ mod1{i1}(length(radius)+2:end)];
    end
    i1=strncmp([theta ''''],mod1,length(theta)+1);
    if any(i1)
        dtheta_dt_=sprintf('d%s_dt_',theta);
        mod1{i1}=[dtheta_dt_ mod1{i1}(length(theta)+2:end)];
    end
    iscomment=strncmp('%',mod1,1);
    i=find(~iscomment,1);
    mod1=[mod1(1:i-1) ,...
        '%model adapted by grind',...
        {sprintf('%s=sqrt(%s^2+%s^2);',radius,xvar,yvar),...
        sprintf('%s=atan2(%s,%s);',theta,yvar,xvar),...
        '%original polar model:'},mod1(i:end),...
        {'%new cartesian model:',...
        sprintf('%s''=%s/%s*%s-%s*%s;',xvar,xvar,radius,dr_dt_,yvar,dtheta_dt_),...
        sprintf('%s''=%s/%s*%s+%s*%s;',yvar,yvar,radius,dr_dt_,xvar,dtheta_dt_)}];
    %fprintf('%s\n',mod1{:})
end
