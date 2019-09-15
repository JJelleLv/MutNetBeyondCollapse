%ENTERJAC   Enter equations of the Jacobian matrix
%   Enter the analytical Jacobian matrixes, consisting of the derivatives
%   of the right-hand sides of all differential/difference equations with respect to each state
%   variable, parameters and sensitivity.
%   This is used to calculate the linearized system in equilibria (stability) and for 
%   continuation of bifurcations. This way you can study the stability of equilibria.
%
%   Usage:
%   ENTERJAC - the user is prompted to give the elements of the Jacobian
%   ENTERJAC({ J11,J12;J21,J22}) - the whole matrix is entered in one message 
%   (J11 is the first string with an equation.
%   ENTERJAC 1 1 J11 - each element is entered separately
%   ENTERJAC('jacobian',@myjac) - use user supplied function (myjac):
%   this function should look like this:
%   jac=function myjac(t,x,par1,par2,par3...) t=time,x=vector of statevars, par1-parn the names of the parameters.
%   Please note that they need to be in the same order as g_grind.pars (will not be checked)!
%
%   ENTERJAC('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'colnr' [integer>0] - column number for input
%     'hessian' [function_handle or cell] - the Hessian 3d matrix as cell of strings or function handle.
%     'implicit' [logical] - the Jacobian for a implicit model, to be used for ODE15i
%     'jacelement' [string] - single element of Jacobian
%     'jacobian' [function_handle or cell] - the whole Jacobian as cell of strings or function handle.
%     'jacobianp' [function_handle or cell] - the Jacobian of parameters (columns) as cell of strings or function handle.
%     'jacobian_y' [function_handle or cell] - the Jacobian of an implicit model.
%     'jacobian_yp' [function_handle or cell] - the Jacobian with respect to the derivatives of an implicit model.
%     'maxorder' [integer>=0 and integer<=5] - limit order of derivatives to maxorder (default=5)
%     'maxtime' [number] - maximum time for deriving the derivatives (default NaN)
%     'rownr' [integer>0] - row number for input
%   ENTERJAC('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-?' - list all Jacobians.
%     '-?der3' - list the 3rd derivatives
%     '-?der4' - list the 4th derivatives
%     '-?der5' - list the 5th derivatives
%     '-?hess' - list the Hessian matrix
%     '-?Hessianp' - list the Hessian with respect to parameters
%     '-?jac' - list the Jacobian only.
%     '-?jacp' - list the Jacobian with respect to parameters
%     '-c' - clear all symbolic Jacobians
%     '-s' - the Jacobian is calculated using the Symbolic toolbox.
%     '-t' - test if the analytical Jacobian is different from the numerical one.
%     '-update' - update the handles for the implicit Jacobian.
%
%   See also eigen, findeq, lyapspect
%
%   Reference page in Help browser:
%      <a href="matlab:commands('enterjac')">commands enterjac</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function [res, maxdif] = enterjac(varargin)
global g_grind;
% if ~isempty(g_grind)
%     n = g_grind.statevars.dim;
% else
%     n = 0;
% end
fieldnams={'jacobian', 'f#c', 'the whole Jacobian as cell of strings or function handle','';...
   'jacobianp', 'f#c', 'the Jacobian of parameters (columns) as cell of strings or function handle.','';...
   'hessian', 'f#c', 'the Hessian 3d matrix as cell of strings or function handle.','';...
   'jacobian_y', 'f#c', 'the Jacobian of an implicit model.','';...
   'jacobian_yp', 'f#c', 'the Jacobian with respect to the derivatives of an implicit model.', '';...
   'maxtime', 'n', 'maximum time for deriving the derivatives (default NaN)',NaN;...
   'maxorder', 'i>=0&i<=5', 'limit order of derivatives to maxorder (default=5)',5;...
   'implicit', 'l', 'the Jacobian for a implicit model, to be used for ODE15i',false;...
   'rownr', 'i>0', 'row number for input',1;...
   'colnr', 'i>0', 'column number for input',1;...
   'jacelement', 's', 'single element of Jacobian',''}';
args=i_parseargs(fieldnams,['if(hasoption(''-update'')),deffields=''implicit'';',...
    'elseif(nargs==1&&iscell(args{1}))',...
    'deffields=''jacobian'';',....
    'else,deffields=''rownr,colnr,jacelement'';end;'],'-?jacp,-?Hessianp,-?hess,-?der3,-?der4,-?der5,-?jac,-?,-c,-update,-t,-s',varargin);
i_parcheck;
if strcmp(args.opts, '-update')
    %use analytical Jacobian for some ode solvers, update the function
    %handles
    g_grind.solver.opt.Jacobian=i_getodehandle('Jacobian');
    return;
end
if isfield(args,'jacobian')
    if all(size(args.jacobian)==g_grind.statevars.dim)||isempty(args.jacobian)
        g_grind.syms.Jacobian = i_enterjacdlg('str2jac',args.jacobian);
    elseif isa(args.jacobian,'function_handle')
        g_grind.syms.Jacobian=args.jacobian;
    else
        error('grind:enterjac:dimension','Error entering the Jacobian matrix, the number of elements does not match the model');
    end
end
if (isfield(args,'jacobian_y')||isfield(args,'jacobian_yp'))&& ~g_grind.solver.isimplicit
    error('grind:enterjac','The fields "jacobian_y" and "jacobian_yp" can only be used for implicit models');
end
if isfield(args,'jacobian_y')
    if all(size(args.jacobian_y)==g_grind.statevars.dim)||isempty(args.jacobian_y)
        g_grind.syms.Jacobian_y = i_enterjacdlg('str2jac',args.jacobian_y);
    elseif isa(args.jacobian_y,'function_handle')
        g_grind.syms.Jacobian_y=args.jacobian_y;
    else
        error('grind:enterjac:dimension','Error entering the Jacobian_y matrix, the number of elements does not match the model');
    end
end
if isfield(args,'jacobian_yp')
    if all(size(args.jacobian_yp)==g_grind.statevars.dim)||isempty(args.jacobian_yp)
        g_grind.syms.Jacobian_yp = i_enterjacdlg('str2jac',args.jacobian_yp);
    elseif isa(args.jacobian_yp,'function_handle')
        g_grind.syms.Jacobian_yp=args.jacobian_yp;
    else
        error('grind:enterjac:dimension','Error entering the Jacobian_yp matrix, the number of elements does not match the model');
    end
end
if isfield(args,'jacobianp')
    if size(args.jacobianp,1)==g_grind.statevars.dim&&size(args.jacobianp,2)==numel(g_grind.pars)
        g_grind.syms.Jacobianp = i_enterjacdlg('str2jac',args.jacobianp);
    elseif isa(args.jacobianp,'function_handle')
        g_grind.syms.Jacobianp=args.jacobianp;
    else
        error('grind:enterjac:dimension','Error entering the Jacobianp matrix, the number of elements does not match the model');
    end
end
if isfield(args,'hessian')
    if all(size(args.hessian)==g_grind.statevars.dim)
        g_grind.syms.Hessian = i_enterjacdlg('str2jac',args.hessian);
    elseif isa(args.hessian,'function_handle')
        g_grind.syms.Hessian=args.hessian;
    else
        error('grind:enterjac:dimension','Error entering the Hessian matrix, the number of elements does not match the model');
    end
end
if isfield(args,'jacelement')
    g_grind.syms.Jacobian{args.rownr,args.colnr}=args.jacelement;
end
if any(strcmp(args.opts, '-c'))
    g_grind.syms.fulleq = {};
    g_grind.syms.Jacobian = {};
    g_grind.syms.Hessian = {};
    g_grind.syms.der3 = {};
    g_grind.syms.der4 = {};
    g_grind.syms.der5 = {};
    g_grind.syms.Jacobianp = {};
    g_grind.syms.Hessianp = {};
    g_grind.syms.Sensitivp = {};
    g_grind.syms.errormsg = {};
    if g_grind.solver.isimplicit
        g_grind.syms.Jacobian_y={};
        g_grind.syms.Jacobian_yp={};
    end
end
if any(strcmp(args.opts, '-?jac'))
    states=i_statevars_names;
    fstates=  regexp(sprintf('f_%d,',1:g_grind.statevars.dim),',','split');
    fstates=fstates(1:end-1);
    plotjac('Jacobian',g_grind.syms.Jacobian,fstates,states);
    disp('');
    if g_grind.solver.isimplicit
        fstates=  regexp(sprintf('eq_%d,',1:g_grind.statevars.dim),',','split');
        fstates=fstates(1:end-1);
        plotjac('Jacobian_y',g_grind.syms.Jacobian_y,fstates,states);
        dstates=regexp(sprintf('%s'',',states{:}),',','split');
        dstates=dstates(1:end-1);
        plotjac('Jacobian_yp',g_grind.syms.Jacobian_yp,fstates,dstates);
    end
end
if any(strcmp(args.opts, '-?hess'))
    states=i_statevars_names;
    fstates=  regexp(sprintf('f_%d,',1:g_grind.statevars.dim),',','split');
    if iscell(g_grind.syms.Hessian)&&numel(g_grind.syms.Hessian)>1
        for j=1:length(states)
            var3=states{j};
            fprintf('Hessian(:,:,%d) derivatives of Jacobian with respect to %s\n',j,var3)
            plotjac('Hessian',g_grind.syms.Hessian(:,:,j),fstates(1:end-1),states,var3);
            disp('');
        end
    else
        plotjac('Hessian',g_grind.syms.Hessian,fstates(1:end-1),states,states{1});
    end
end
if any(strcmp(args.opts, '-?der3'))
    plotjac('der3',g_grind.syms.der3);
end
if any(strcmp(args.opts, '-?der4'))
    plotjac('der4',g_grind.syms.der4);
end
if any(strcmp(args.opts, '-?der5'))
    plotjac('der5',g_grind.syms.der5);
end
if any(strcmp(args.opts, '-?jacp'))
    fstates=  regexp(sprintf('f_%d,',1:g_grind.statevars.dim),',','split');
    plotjac('Jacobianp',g_grind.syms.Jacobianp,fstates(1:end-1),g_grind.pars);
    disp('');
end
if any(strcmp(args.opts, '-?Hessianp'))
    %hessian to parameters of the Jacobian matrix
    plotjac('Hessianp',g_grind.syms.Hessianp);
end
if any(strcmp(args.opts, '-?'))
    if ~g_grind.solver.isdiffer&&~g_grind.solver.isimplicit&&isfield(g_grind.syms,'fulleq')
        states=i_statevars_names;
        plotjac('Single line equations',g_grind.syms.fulleq(:),states,{'t'});
    end
    if iscell(g_grind.syms.Jacobian)&&~isempty(g_grind.syms.Jacobian)
        if g_grind.solver.isimplicit
            fprintf('Symbolic Jacobians defined  <a href="matlab:enterjac(''-?jac'')">3x{%dx%d}</a>\n',size(g_grind.syms.Jacobian,1),size(g_grind.syms.Jacobian,2));
        else
            fprintf('Symbolic Jacobian defined  <a href="matlab:enterjac(''-?jac'')">{%dx%d}</a>\n',size(g_grind.syms.Jacobian,1),size(g_grind.syms.Jacobian,2));
        end
    else
        plotjac('Jacobian',g_grind.syms.Jacobian,[],[]);
    end
    if iscell(g_grind.syms.Hessian)&&~isempty(g_grind.syms.Hessian)||isstruct(g_grind.syms.Hessian)
        fprintf('Symbolic Hessian defined   <a href="matlab:enterjac(''-?hess'')">{%dx%dx%d}</a>\n',g_grind.syms.Hessian.size);
    else
        plotjac('Hessian',g_grind.syms.Hessian,[],[]);
    end
    if iscell(g_grind.syms.der3)&&~isempty(g_grind.syms.der3)||isstruct(g_grind.syms.der3)
        fprintf('Symbolic der3 defined      <a href="matlab:enterjac(''-?der3'')">{%dx%dx%dx%d}</a>\n',g_grind.syms.der3.size);
    else
        plotjac('der3',g_grind.syms.der3,[],[]);
    end
    if iscell(g_grind.syms.der4)&&~isempty(g_grind.syms.der4)||isstruct(g_grind.syms.der4)
        fprintf('Symbolic der4 defined      <a href="matlab:enterjac(''-?der4'')">{%dx%dx%dx%dx%d}</a>\n',g_grind.syms.der4.size);
    else
        plotjac('der4',g_grind.syms.der4,[],[]);
    end
    if iscell(g_grind.syms.der5)&&~isempty(g_grind.syms.der5)||isstruct(g_grind.syms.der5)
        fprintf('Symbolic der5 defined      <a href="matlab:enterjac(''-?der5'')">{%dx%dx%dx%dx%dx%d}</a>\n',g_grind.syms.der5.size);
    else
        plotjac('der5',g_grind.syms.der5,[],[]);
    end
    if iscell(g_grind.syms.Jacobianp)&&~isempty(g_grind.syms.Jacobianp)
        fprintf('Symbolic Jacobianp defined <a href="matlab:enterjac(''-?jacp'')">{%dx%d}</a>\n',size(g_grind.syms.Jacobianp,1),size(g_grind.syms.Jacobianp,2));
    else
        plotjac('Jacobianp',g_grind.syms.Jacobianp,[],[]);
    end
    if iscell(g_grind.syms.Hessianp)&&~isempty(g_grind.syms.Hessianp)||isstruct(g_grind.syms.Hessianp)
        fprintf('Symbolic Hessianp defined  <a href="matlab:enterjac(''-?Hessianp'')">{%dx%dx%d}</a>\n',g_grind.syms.Hessianp.size);
    else
        plotjac('Hessianp',g_grind.syms.Hessianp,[],[]);
    end
end
if any(strcmp(args.opts, '-t'))
    if isempty(g_grind)||~isfield(g_grind,'syms')||isempty(g_grind.syms.Jacobian)
        if i_hastoolbox('symbolic')
            disp('Cannot test analytical Jacobian,  use <a href="matlab:enterjac(''-sym'')">enterjac -sym</a> to generate analytical Jacobian');
        else
            disp('Cannot test analytical Jacobian, enter the Jacobian first with <a href="matlab:enterjac">enterjac</a> ');
        end
        if nargout==0
            res1 = false;
            maxdif1 = nan;
        end
    elseif g_grind.solver.isimplicit
        N0 = i_initvar;
        [~,~] = i_implicit_Jac(0,N0,N0);
        res1=true;
        maxdif1=0;
    else
        N0 = i_initvar;
        try
            i_calcjac(1, g_grind.solver.iters, N0);
            Jnum = i_calcjac(1, g_grind.solver.iters, N0); %run twice to use improved factors
            Janal = i_calcjac(0, g_grind.solver.iters, N0);
            maxdif1 = max(max(abs(Jnum - Janal)));
            if nargout  > 0
                res1 = maxdif1 < 1E-4;
            else
                if maxdif1 < 1E-4
                    fprintf('OK: max difference numeric and analytical Jacobian = %g\n', maxdif1);
                else
                    fprintf('Error: max difference numeric and analytical Jacobian = %g\n', maxdif1);
                end
            end
        catch err
            disp('cannot run analytical Jacobian');
            enterjac('-?')
            rethrow(err);
        end
    end
    if nargout>0
        res = res1;
        maxdif = maxdif1;
    end
    return;
end
if any(strcmp(args.opts, '-s'))
    %Help from symbolic toolbox
    if ~isfield(args,'maxtime')
        args.maxtime=30;
    end
    if ~isfield(args,'maxorder')
        %default values for the max order of the derivatives
        %der3-der5 get very slow with increasing dimensiona
        %more than exponential increasing
        if g_grind.statevars.dim<=5
            args.maxorder=5; %der5
            %5^2*3^4=2025
        elseif g_grind.statevars.dim<=7
            args.maxorder=4; %der4
            %7^2*4^3=3136
        elseif g_grind.statevars.dim<=10
            args.maxorder=3; %der3
            %10^2*6^2=3600
        elseif g_grind.statevars.dim<=18
            args.maxorder=2;
            %18^2*10 = 3240
        else
            args.maxorder=1;
        end
    end
    runwhich=ones(1,7);
    % runwhich makes it possible to run part of the Jacobians (no check for dependency)
    % vector with 7 elements 1=yes 0=no  order: Jacobian, Jacobianp Hessian Hessianp der3 der4 der5
    switch args.maxorder
        case 4
            runwhich(7)=0;
        case 3
            runwhich(6:7)=0;
        case 2
            runwhich(5:7)=0;
        case 1
            runwhich(3:7)=0;
        case 0
            runwhich=zeros(1,7);
    end
    if isempty(g_grind.pars)
        runwhich(2)=0;
        runwhich(4)=0;
    end
    
    if g_grind.statevars.vector
        error('grind:enterjac:vector','Model with vector/matrix notation cannot be differentiated')
    end
    if g_grind.solver.haslags
        error('grind:enterjac:lags','Model with time lags cannot be differentiated')
    end
    i_parcheck;
    i_model2mupad(false,g_grind.syms,args.maxtime,runwhich);
    
    if ~isempty(g_grind.syms.errormsg)
        error('grind:enterjac','%s\n%s',g_grind.syms.errormsg{:});
    end
    %check if they are runnable
    
    f={'Jacobian','Jacobianp','Hessian','Hessianp','der3','der4','der5'};
    N0=i_initvar;
    pars=struct2cell(par('-v'));
    pars=pars(2:end);
    if isempty(pars)
        parlist='';
    else
        parlist=sprintf(',%s',g_grind.pars{:});
    end
    for i=1:length(f)
        if ~isempty(g_grind.syms.(f{i}))
            han=i_getodehandle(f{i},parlist);
            try
                if g_grind.solver.isimplicit
                    han(0,N0,N0);
                else
                    han(0,N0,pars{:});
                end
            catch
                fprintf('Cannot run analytical %s\n',f{i});
                g_grind.syms.(f{i})=[];
            end
        end
    end
    enterjac('-?')
    try
        [isok, differ] = enterjac('-t');
    catch
        isok=false;
        differ=nan;
    end
    if ~isok
        if ~isnan(differ)
            warning('grind:enterjac','Analytical and numerical Jacobians differ, max difference=%g',differ)
        else
            %   g_grind.syms.Jacobian = {};
            if ~isa(g_grind.syms.Jacobian,'function_handle')
                % g_grind.syms.Jacobian = {};
            end
            if ~isa(g_grind.syms.Jacobianp,'function_handle')
                %  g_grind.syms.Jacobianp = {};
            end
            if ~isa(g_grind.syms.Hessian,'function_handle')
                %  g_grind.syms.Hessian = {};
            end
            g_grind.syms.errormsg{end + 1} = 'Cannot run analytical Jacobian';
            error('grind:enterjac','Cannot run analytical Jacobian');
        end
    else
        p = solver('-properties');
        if p.usesJac
            enterjac('-update'); %use the jacobian for some implicit solvers
        end
    end
    return;
end

if nargin == 0
    i_enterjacdlg; %dialog box where you can enter and view all equations
end
%enterjac('-?')
function s=index2sub(siz,index)
a=cell(size(siz));
[a{1:end}]=ind2sub(siz,index);
s=sprintf('%d,',a{:});
s=s(1:end-1);
function jac=plotjac(name,jac,dx,dy,dhess)
if nargin<5
    dhess='';
end
if isempty(jac)
    fprintf('No %s defined\n',name)
elseif isa(jac,'function_handle')
    fprintf('Symbolic %s defined   @%s\n',name,func2str(jac))
elseif isstruct(jac)
    s1=sprintf('%d,',jac.size);
    fprintf('%s = zeros(%s)\n',name,s1(1:end-1))
    if isfield(jac,'unique')
        for i=1:numel(jac.unique)
            u=jac.unique(i);
            for j=1:size(u.indices)
                c=i_enterjacdlg('jac2str',{u.equation});
                fprintf('%s(%s) = %s\n',name,index2sub(jac.size,u.indices(j)),c{1})
            end
        end
    end
else
    jac=i_enterjacdlg('jac2str',jac);
    if isempty(dhess)
        fprintf('%s:\n',name)
    end
    for i = 1:length(dx)
        for j = 1:length(dy)
            if ~isempty(dhess)
                fprintf('d2(%s)/d(%s)/d(%s) = %s\n',dx{i},dy{j},dhess,jac{i,j});
            else
                fprintf('d(%s)/d(%s) = %s\n',dx{i},dy{j},jac{i,j});
            end
        end
    end
end