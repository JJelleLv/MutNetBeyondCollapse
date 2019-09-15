%FORWARD_STABIL   Analyse the forward basin stability of a bistable system
%   Analyse by simulation which parameter changes are "forward basin unstable" i.e. they can potentially 
%   result in rate tipping. The result of the last run <a href="matlab:help paranal">paranal</a> is analysed (can be 2 dimensional).
%   Paranal should be run in two directions (if this is not done, paranal -1 is run automatically).
%   The user should either give the parameter values of a target point or the initial point. All
%   other points of paranal are tested for forward basin stability by simulation.
%  
%   Citation for the concept of forward basin stability:
%   Ashwin, P., C. Perryman, and S. Wieczorek. 2017. Parameter shifts for nonautonomous systems in low 
%   dimension: Bifurcation- and rate-induced tipping. Nonlinearity 30:2185-2210 doi: <a href="http://dx.doi.org/10.1088/1361-6544/aa675b">10.1088/1361-6544/aa675b</a>.
%
%  
%
%   Usage:
%   RES=FORWARD_STABIL - opens an inputbox to enter all settings. RES is a structure with the
%   results:
%     RES.istarget - logical value that is true if the point is a target and false if it is an initial condition.
%     RES.parvalue - values of the parameters of the point.
%     RES.forward - consider the forward or backward run of <a href="matlab:help paranal">paranal</a>.
%     RES.pars - list with parameter names.
%     RES.initial or RES.target - structure with the initial point or the target.
%     RES.points - structure with the stability of parameter values.
%        RES.points.parvalues - all unique parameter values.
%        RES.points.Y - all equilibrium values.
%        RES.points.unstable -logical value that is true if the parameter combination is forward unstable
%        RES=FORWARD_STABIL(parvalue,forward) - order of default arguments: parvalue(s) and forward(true/false).
%   FORWARD_STABIL('argname',argvalue,...) - Valid argument name-value pairs [with type]:
%     'forward' [logical] - true if the forward run of paranal is considered.
%     'istarget' [logical] - true if the point is the target, false if it is the initial condition
%     'parvalue' [number] - the values of the parameters of the point
%     'runno' [integer] - plot previous run number RUNNO (default the last run)
%   FORWARD_STABIL('-opt1','-opt2',...) - Valid command line options:
%     '-p' - plot the last run.
%     
%         
%   See also paranal, paranal2d
%
%   Reference page in Help browser:
%      <a href="matlab:commands('forward_stabil')">commands forward_stabil</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function res=forward_stabil(varargin)
global g_paranal g_grind g_Y
fieldnams={'parvalue', 'n', 'the values of the parameters of the point',0;...
   'forward', 'l', 'true if the forward run of paranal is considered.',true;...
   'istarget', 'l', 'true if the point is the target, false if it is the initial condition',true;...
   'runno', 'i', 'plot previous run number RUNNO (default the last run)',[]}';
args=i_parseargs(fieldnams,'parvalue,forward','-p',varargin);
if any(strcmp(args.opts,'-p'))
    %plot last run
    try
        if isfield(args,'runno')
            res1=g_paranal.forward_stabil{args.runno};
        else
            res1=g_paranal.forward_stabil{end};
        end
    catch err
        if strcmp(err.identifier,'MATLAB:nonExistentField')
            error('grind:forward_stabil','Cannot plot the results as <a href="matlab:forward_stabil">forward_stabil</a> has not been run before');
        else
            rethrow(err);
        end
    end
    if length(res1.parvalue)==2
        %statevars=outfun('-p','fun',g_grind.statevars.names,'runfield',res.points.runfield);
        pars=outfun('-p','fun',g_paranal.run.pars,'runfield',res1.points.runfield);
        X=reshape(pars(res1.points.index,1),res1.points.size);
        X=X(:,1);
        Y=reshape(pars(res1.points.index,2),res1.points.size);
        Y=Y(1,:)';
        Z=double(reshape(res1.points.unstable==0,res1.points.size)).';
        if Y(end)<Y(1)
 %           Z=fliplr(Z);
            Y=flipud(Y);
        end
        if X(1)>X(end)
 %           Z=flipud(Z);
            X=flipud(X);
        end
        h=imagesc(X,Y,Z);
        set(gca,'ydir','normal')
        hold on;
        shading flat;
        i_plotdefaults();
        hax=get(h,'parent');
        hfig=get(hax,'parent');
        if res1.istarget
         set(hfig,'colormap',[0.7 0.7 0.7;1 1 1]);
           h=plot(hax,res1.target.parvalues(1),res1.target.parvalues(2),'ro');
            set(h,'MarkerFaceColor','r');
        else
         set(hfig,'colormap',[0.8 0.8 0.8;1 1 1]);
           h=plot(hax,res1.initial.parvalues(1),res1.initial.parvalues(2),'ko');
            set(h,'MarkerFaceColor','k');
        end
    elseif length(res1.parvalue)==1
        Xs=outfun('-p','fun', g_grind.paranal.plots{1}.xaxis,'runfield',res1.points.runfield);
        X=Xs(res1.points.index(res1.points.unstable==0),1);
        Ys=outfun('-p','fun', g_grind.paranal.plots{1}.yaxis,'runfield',res1.points.runfield);
        Y=Ys(res1.points.index(res1.points.unstable==0),1);
        hfig=i_makefig('paranal');
        hold on
        plot(X,Y,'ro');
    end
    return
end
if isempty(g_paranal)||isempty(g_paranal.run)
    error('grind:forward_stabil','First run <a href="matlab:paranal">paranal</a> or <a href="matlab:paranal2d">paranal2d</a> before you can analyse forward basin stability')
end
if ~isfield(args,'istarget')
    args.istarget=true;
end

if ~isfield(args,'parvalue')
    if length(g_paranal.run.pars)==1
        answer= {'','true','true'};
    else
        answer= {'','','true','true'};
    end
    if isfield(g_paranal,'forward_stabil')
        for i=1:length(g_paranal.forward_stabil{end}.parvalue)
            answer{i}=num2str(g_paranal.forward_stabil{end}.parvalue(i));
        end
        if g_paranal.forward_stabil{end}.forward
            answer{length(answer)-1}='true';
        else
            answer{length(answer)-1}='false';
        end
        if g_paranal.forward_stabil{end}.istarget
            answer{length(answer)}='true';
        else
            answer{length(answer)}='false';
        end
    end
    if length(g_paranal.run.pars)==2
        answer=inputdlg({sprintf('Reference point: value of "%s"',g_paranal.run.pars{1}),sprintf('Reference point: value of "%s"',g_paranal.run.pars{2}),'Reference attractor: forward run? (true/false)','Reference point is target (true) or initial condtion (false)'},'forward basin stability',1,answer);
        args.parvalue(2)=str2num(answer{2}); %#ok<ST2NM>
        args.parvalue(1)=str2num(answer{1}); %#ok<ST2NM>
        args.forward=eval(answer{3});
        args.istarget=eval(answer{4});
    else
        answer=inputdlg({sprintf('Reference point: value of "%s"',g_paranal.run.pars{1}),'Reference attractor: forward run? (true/false)','Reference point is target (true) or initial condtion (false)'},'forward basin stability',1,answer);
        args.parvalue(1)=str2num(answer{1}); %#ok<ST2NM>
        args.forward=eval(answer{2});
        args.istarget=eval(answer{3});
    end
end

if numel(args.parvalue)~=numel(g_paranal.run.pars)
    error('grind:forward_stability','The number of parameter values does not match the last paranal run');
end
if isempty(g_paranal.prevrun)
    paranal('-1')
end
if ~isfield(args,'forward')
    args.forward=true;
end
args.parvalue=args.parvalue(:).';

[forward,backward]=anal_paranal;
f=[find(any(diff(forward.parvalues)~=0,2));size(forward.parvalues,1)];
[~,currndx]=min(sum((forward.parvalues-repmat(args.parvalue,size(forward.parvalues,1),1)).^2,2));
parameter_state=forward.parvalues(currndx,:);
if args.forward
    YYcurr=forward.Y(currndx,:);
    YYalt=backward.Y(currndx,:);
    curdir=forward;
    curalt=backward;
else
    YYcurr=backward.Y(currndx,:);
    YYalt=forward.Y(currndx,:);
    curdir=backward;
    curalt=forward;
end
isunstable=double(false(size(f)));
res=rmfield(args,'opts');
res.pars=g_paranal.run.pars;
if args.istarget
    %the point is the target
    for i=1:length(parameter_state)
        assignin('base',g_paranal.run.pars{i},parameter_state(i));
    end
    for j=1:length(f)
        i=f(j);
        if curdir.ndx(i)&&i~=currndx && ~any(isnan(curdir.Y(i,:)))
            N0=curdir.Y(i,:).';
            i_keep(N0);
            stabil('-s')
            Nend=g_Y(end,:).';
            if sqrt(sum((Nend-YYcurr).^2))>sqrt(sum((Nend-YYalt).^2))
                isunstable(j)= true;
            end
            %fprintf('%g %g - %g %g   %d\n',curdir.parvalues(i,1),curdir.parvalues(i,2),N0,Nend,isunstable(j))
        end
    end
    res.target.parvalues=parameter_state;
    res.target.Y=YYcurr;
    res.target.forward=args.forward;
else
    %the point is the source
    for j=1:length(f)
        i=f(j);
        if curdir.ndx(i)&&curalt.ndx(i)&&i~=currndx && ~any(isnan(curdir.Y(i,:)))&& ~any(isnan(curalt.Y(i,:)))
            for k=1:length(parameter_state)
                assignin('base',g_paranal.run.pars{k},curdir.parvalues(i,k));
            end
            %N0=curdir.Y(i,:).';
            i_keep(YYcurr);
            stabil('-s')
            Nend=g_Y(end,:).';
            if sqrt(sum((Nend-curdir.Y(i,:).').^2))<sqrt(sum((Nend-curalt.Y(i,:).').^2))
                isunstable(j)= true;
                disp(curdir.parvalues(f(j),:));
            end
        end
    end
    res.initial.parvalues=parameter_state;
    res.initial.Y=YYcurr;
end
if ~args.forward
    res.points.runfield='prevrun';
    res.points.index=flipud(f);
else
    res.points.runfield='run';
    res.points.index=f;
end
%res.points.parvalues=curdir.parvalues(f,:);
%res.points.Y=curdir.Y(f,:);
res.points.unstable=isunstable;
if length(res.pars)==1
    res.points.size= [g_grind.paranal.dlg.steps(1)+1,1];
else
    res.points.size= g_grind.paranal.dlg.steps+1;
end

if ~isfield(g_paranal,'forward_stabil')||isempty(g_paranal.forward_stabil)
    g_paranal.forward_stabil{1}=res;
else
    i=1;
    while (i<=length(g_paranal.forward_stabil)) && ((g_paranal.forward_stabil{i}.forward~=args.forward) ...
            || (g_paranal.forward_stabil{i}.istarget~=args.istarget) || any(g_paranal.forward_stabil{i}.parvalue~=res.parvalue))
       i=i+1;
    end 
    if i<=length(g_paranal.forward_stabil)
        g_paranal.forward_stabil(i)=[];
    end
    g_paranal.forward_stabil{end}=res;
end

function [forward,backward]=anal_paranal
%analyse the result of 1D or 2D paranal run. It is needed that they are run
%forward and backward (order is not important)
global g_paranal g_grind;
statevars=g_grind.statevars.names;
%needs to be adapted for vector models
forward.parvalues=outfun('-p',g_paranal.run.pars);
forward.Y=outfun('-p',statevars);
backward.parvalues=outfun('-p',g_paranal.run.pars,'runfield','prevrun');
backward.Y=outfun('-p',statevars,'runfield','prevrun');

if forward.parvalues(end,1)<forward.parvalues(1,1)
    %exchange to match forward and backward run
    h=backward;
    backward=forward;
    forward=h;
end
backward.parvalues=flipud(backward.parvalues);
backward.Y=flipud(backward.Y);
if size(forward.parvalues,1)>1
    breaks=[0;find(diff(forward.parvalues(:,1))>0);size(forward.Y,1)];
else
    breaks=[0;size(forward.Y,1)];
end
diff_fb=sqrt(sum((forward.Y-backward.Y).^2,2));
hasASS=double(diff_fb>max(diff_fb)*0.01);
for j=2:length(breaks)
    ndx=breaks(j-1)+1:breaks(j);
    h=hasASS(ndx);
    res=h;
    fw=forward.Y(ndx,:);
    bw=backward.Y(ndx,:);
    f=[find(diff(h)~=0);length(h)];
    l=1;
    for k=1:length(f)-1
        if h(f(k))
            %from ASS to no ASS
            varnoASS=(fw(f(k)+1,:)+bw(f(k)+1,:))/2;
            if norm(fw(f(k),:)-varnoASS)<norm(bw(f(k),:)-varnoASS)
                res(f(k)+1:f(k+1))=2; %forward
            else
                res(f(k)+1:f(k+1))=3;
            end
        else
            %fromd no ASS to ASS
            varnoASS=(fw(f(k),:)+bw(f(k),:))/2;
            if norm(fw(f(k)+1,:)-varnoASS)<norm(bw(f(k)+1,:)-varnoASS)
                res(l:f(k))=2; %forward
            else
                res(l:f(k))=3;
            end
            
        end
        l=f(k);
    end
    hasASS(ndx)=res;
end
forward.ndx=hasASS==1|hasASS==2;
backward.ndx=hasASS==1|hasASS==3;



