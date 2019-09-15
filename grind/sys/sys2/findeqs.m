%FINDEQS   Find all equilibria
%   Try to find all equilibria (may fail or find double equilibria). 
%   It uses the simplex optimizing method and tries various random initial conditions.
%   Any equlibria found with <a href="matlab:help findeq">findeq</a> are added to the list.
%
%   Usage:
%   FINDEQS - find equilibria select one.
%   FINDEQS -r - refreshes the list of equilibria.
%   FINDEQS -a - try to add to the existing list of equilibria.
%   FINDEQS(ntrial,scale) - ntrial is the number of trials (default=100), scale is the 
%   scale of the initial conditions (default=50).
%   eqlist=FINDEQS - saves the list of equilibria to a structure (fields 'id' (IC1-n), 'label','IC','msg',message, 
%   'N0', vector with state variables).
%   FINDEQS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'complex' [logical] - Search an equilibrium if state variable(s) can be complex numbers.
%     'display' [off | on | iter-detailed | iter | final | final-detailed | notify | notify-detailed] - Different display options 
%    of fminsearch
%     'id' [integer>0] - Select the equilibrium with this id number
%     'keep' [logical] - Keep the last value
%     'mask' [logical and length(logical)==ndim or empty] - Mask for state variables [1 0 1 0] (1=optimize 0=keep fixed) (you can find a nullcline value or partial
%    equilibrium this way, leave empty for no mask)
%     'maxfunevals' [integer>0] - Max number of function evaluations
%     'maxiter' [integer>0] - Max number of iterations
%     'maxtime' [number>0|isnan(number)] - Maximum time for calculation in sec. NaN = ignore.
%     'method' [fminsearch | fsolve | trust-region-dogleg | trust-region | levenberg-marquardt] - Algorithm of optimization 
%    fminsearch or fsolve (for fsolve, trust-region or levenberg-marquardt is the optimization toolbox <a href="matlab:commands toolboxes">required</a>)
%     'n0' [number] - Initial values of state variables a vector of all state variables leave empty for current values []
%     'ntrial' [integer>0] - Number of attempts to find equilibria (default 100)
%     'p0' [number] - Initial values of the parameters (vector of all values) []
%     'scale' [number>0] - Scale of the initial conditions (default=50)
%     'toleq' [number>0] - Tolerance for distinguishing different equilibria
%     'tolfun' [number>0] - Tolerance of function
%     'tolx' [number>0] - Tolerance of X
%     'userfun' [function_handle] - Add a user function for optimizing []
%   FINDEQS('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-a' - add new equilibria to the list of equilibria
%     '-add_curr' - add the current state to the list of equilibria (used internally)
%     '-c' - clear the current list of equilibria
%     '-l' - list the equilibria that are found.
%     '-o' options: display a user interface with options-
%     '-r' - rerun FINDEQS
%
%
%   See also findeq, eigen
%
%   Reference page in Help browser:
%      <a href="matlab:commands('findeqs')">commands findeqs</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function [eqlist,eqs] = findeqs(varargin)
%(ntrial, scal)
global g_grind;
if isfield(g_grind,'solver')
    iscomplex=g_grind.solver.iscomplex;
else
    iscomplex=false;
end
fieldnams={'complex', 'l', 'Search an equilibrium if state variable(s) can be complex numbers.',iscomplex;...
   'display', 'e[off|on|iter-detailed|iter|final|final-detailed|notify|notify-detailed]', 'Different display options','notify';...
   'id', 'i>0', 'Select the equilibrium with this id number',1;...
   'keep', 'l', 'Keep the last value',true;...
   'mask', 'l&length(l)==ndim#E', 'Mask for state variables [1 0 1 0] (1=optimize 0=keep fixed) (you can find a nullcline value or partial',[];...
   'maxfunevals', 'i>0', 'Max number of function evaluations',50000;...
   'maxiter', 'i>0', 'Max number of iterations',50000;...
   'maxtime', 'n>0|isnan(n)', 'Maximum time for calculation in sec. NaN = ignore.',NaN;...
   'method', 'e[fm+insearch|fs+olve|trust-region-dogleg|trust-region|l+evenberg-marquardt]', 'Algorithm of optimization','fminsearch';...
   'n0', 'n', 'Initial values of state variables a vector of all state variables leave empty for current values []',[];...
   'ntrial', 'i>0', 'Number of attempts to find equilibria (default 100)',100;...
   'p0', 'n', 'Initial values of the parameters (vector of all values) []',[];...
   'scale', 'n>0', 'Scale of the initial conditions (default=50)',50;...
   'toleq', 'n>0', 'Tolerance for distinguishing different equilibria',1E-7;...
   'tolfun', 'n>0', 'Tolerance of function',1E-6;...
   'tolx', 'n>0', 'Tolerance of X',1E-10;...
   'userfun', 'f', 'Add a user function for optimizing []',[]}';

args=i_parseargs(fieldnams,'ntrial,scale','-c,-r,-add_curr,-a,-l,-o',varargin);
if isfield(g_grind,'findeq')
    findeq_args=mergestructs(g_grind.findeq,rmfield(args,'opts'));
else
    findeq_args=mergestructs(findeq('-default'),rmfield(args,'opts'));
end
if ~isfield(args, 'toleq')
    if isfield(g_grind,'findeqs')&&isfield(g_grind.findeqs,'toleq')
        args.toleq= g_grind.findeqs.toleq;
    else
        args.toleq=1E-7;
    end
else
    findeq_args=rmfield(findeq_args,'toleq');
end
if isfield(args,'id')
    [eqlist1,eqs]=update_eqlst(args.toleq);
    if args.id>size(eqs,1)
        args.id=[];
    end
    eqs=eqs(args.id,:);
    if nargout>0
        eqlist=eqlist1(args.id);
    elseif ~isempty(eqs)
        i_keep(transpose(eqs));
        findeq('maxiter',1)
    end
    return;
end
if any(strcmp(args.opts, '-l'))
    if isfield(g_grind,'findeqs')
        settings=par(struct('opts','-vector'));
        changed=(length(settings)~=length(g_grind.findeqs.settings)) || (sum(settings~=g_grind.findeqs.settings) > 0);
        if changed
            disp('Parameter(s) changed since last findeqs run')
            return;
        end
%    else
%        changed=true;
    end
    if nargout>0
        [eqlist,eqs]=update_eqlst(args.toleq);
    else
        if ~isfield(g_grind,'findeqs')
            disp('No previous findeqs run')
        elseif changed
            disp('Parameters have changed')
        else
            eqlist1=update_eqlst(args.toleq);
            fprintf('%d equilibria found in %d trials: ',length(eqlist1),g_grind.findeqs.ntrials)
            fprintf('[<a href="matlab:findeqs(''-a'')">try more</a>|<a href="matlab:findeqs(''-o'',''-r'')">change settings</a>]\n');
            for i=1:length(eqlist1)
                fprintf('<a href="matlab:findeqs(''id'',%d)">EP%d</a> - %s\n',i,i,eqlist1(i).msg);
            end
            
        end
    end
    return;
end
if any(strcmp(args.opts, '-c'))
    if isfield(g_grind,'findeqs')
        g_grind=rmfield(g_grind,'findeqs');
        disp('Cleared the list of equilibria');
    else
        disp('There is no list of equilibria');
    end
    return;
end
if ~isfield(findeq_args,'display')
    findeq_args.display='off';
end
if ~isfield(args, 'maxtime')||isempty(args.maxtime)
    if isfield(g_grind,'findeq')&&isfield(g_grind.findeq,'maxtime')
        findeq_args.maxtime=g_grind.findeq.maxtime;
    else
        findeq_args.maxtime = NaN;
    end
end

if ~isfield(args, 'ntrial')
    if isfield(g_grind,'findeqs')&&isfield(g_grind.findeqs,'ntrial')
        args.ntrial= g_grind.findeqs.ntrial;
    else
        args.ntrial=100;
    end
else
    findeq_args=rmfield(findeq_args,'ntrial');
end
if ~isfield(args, 'scale')
    if isfield(g_grind,'scale')&&isfield(g_grind.findeqs,'scale')
        args.scale= g_grind.findeqs.scale;
    else
        args.scale=50;
    end
else
    findeq_args=rmfield(findeq_args,'scale');
end
if any(strcmp(args.opts, '-o'))
    %findeq_args=mergestructs(findeq_args, findeq('-default'));
    prompt={ 'Number of trials to find equilibria (ntrial):',...
        'Scale of the perturbations in state variables (Scale):', ...
        'Tolerance for distinguishing equilibria (TolEq):', ...
        'Max number of function evaluations findeq (MaxFunEvals)',...
        'Max number of iterations findeq (MaxIter)',...
        'Max time in sec/NaN=ignore (MaxTime)',...
        'Method of optimization (Method=fminsearch | fsolve | trust-region | levenberg-marquardt)' ,...
        'Tolerance of function findeq (TolFun):', ...
        'Tolerance of X findeq (TolX)'};
    answers={int2str(args.ntrial),num2str(args.scale),num2str(args.toleq),int2str(findeq_args.maxfunevals),int2str(findeq_args.maxiter),num2str(findeq_args.maxtime),...
        findeq_args.method,num2str(findeq_args.tolfun),num2str(findeq_args.tolx)};
    answers = inputdlg(prompt, 'findeqs options', 1, answers);
    if isempty(answers)
        return;
    end
    args.ntrial=str2double(answers{1});
    args.scale=str2double(answers{2});
    args.toleq=str2double(answers{3});
    findeq_args.maxfunevals=str2double(answers{4});
    findeq_args.maxiter=str2double(answers{5});
    findeq_args.maxtime=str2double(answers{6});
    findeq_args.method=answers{7};
    findeq_args.tolfun=str2double(answers{8});
    findeq_args.tolx=str2double(answers{9});
end
do_add_curr=any(strcmp(args.opts, '-add_curr'));
do_add=any(strcmp(args.opts, '-a'));
do_rerun=any(strcmp(args.opts, '-r'));
neqs_before=0;
ntrials=args.ntrial;
settings = par(struct('opts','-vector'));
if isfield(g_grind, 'findeqs')
    parchanged=(length(settings)~=length(g_grind.findeqs.settings)) || (sum(settings~=g_grind.findeqs.settings) > 0);
else
    g_grind.findeqs=struct('ntrials',0,'eqs',[],'settings',settings);
    parchanged=true;
end
g_grind.findeqs.toleq=args.toleq;
if parchanged  %never save equilibria if parameters have been changed
    eqs=[];
    g_grind.findeqs.ntrials=0;
    g_grind.findeqs.settings=settings;
elseif do_add
    eqs=g_grind.findeqs.eqs;
    neqs_before=size(clusteqs(eqs,g_grind.findeqs.toleq, g_grind.solver.iscomplex),1);
elseif do_rerun
    eqs=[];
    g_grind.findeqs.ntrials=0;
else
    eqs=g_grind.findeqs.eqs;
    neqs_before=size(clusteqs(eqs,g_grind.findeqs.toleq,g_grind.solver.iscomplex),1);
    if isfield(g_grind.findeqs,'maxtimereached')&&g_grind.findeqs.maxtimereached
        ntrials=0;
    else
        ntrials=ntrials-g_grind.findeqs.ntrials;
    end
end
g_grind.findeqs.ntrial=args.ntrial;
g_grind.findeqs.scale=args.scale;
if ~isfield(g_grind.findeqs,'maxtimereached')&&~isnan(findeq_args.maxtime)||parchanged||do_rerun||do_add
    g_grind.findeqs.maxtimereached=false;
end
N0 = i_initvar;
if do_add_curr
    findeq_args.maxiter=1;
    %add equilibrium only if almost in equilibrium
    eqs = addeq(eqs,N0,findeq_args);
    if neqs_before<size(eqs,1)
        [~, index] = sort(min(eqs, [], 2) + sum(eqs, 2) / 1000);
        eqs = eqs(flipud(index), :); %sort in an order that the trivial / negative equilibria are last
        g_grind.findeqs.eqs = eqs;
        g_grind.findeqs.settings=settings;
        g_grind.findeqs.ntrials=g_grind.findeqs.ntrials+1;
    end
    return
end
oldN0 = N0;
try
    if ntrials<=0
        eqs = g_grind.findeqs.eqs;
    else
        %  eqs = [];
        tstart=tic;
        [eqs] = addeq(eqs, N0, findeq_args);
        if args.ntrial > 0 && (isnan(findeq_args.maxtime)||toc(tstart)<findeq_args.maxtime)
            hwait=i_waitbar(0, ntrials, 'findeqs','Searching equilibria',0.5,true);
            for i = 1:ntrials
                f = rand(size(N0)) < 0.1;
                f2 = rand(size(N0)) < 0.1;
                f3 = rand(size(N0)) < 0.05;
                N0 = rand(size(N0)) * args.scale;
                N0(f) = 0.001 * N0(f);
                N0(f2) = 100 * N0(f2);
                N0(f3) = -N0(f3);
                [eqs] = addeq(eqs, N0, findeq_args);
                g_grind.findeqs.ntrials=g_grind.findeqs.ntrials+1;
                if ~isnan(findeq_args.maxtime)&&toc(tstart)>findeq_args.maxtime
                    disp('Maximum calculation time reached');
                    g_grind.findeqs.maxtimereached=true;
                    break;
                end
                if mod(i,10) %update the title only 1 time in 10 runs
                    neqs=size(clusteqs(eqs,g_grind.findeqs.toleq,g_grind.solver.iscomplex),1);
                    if neqs_before==0
                        if neqs==1
                            atitle='Found 1 equilibrium';
                        elseif neqs==0
                            atitle='Searching equilibria';
                        else
                            atitle=sprintf('Found %d equilibria',neqs);
                        end
                    else
                        n=neqs-neqs_before;
                        switch n
                            case 0
                                atitle=sprintf('Found no new equilibria (%d)',neqs);
                            case 1
                                atitle=sprintf('Found %d new equilibrium (%d)',n,neqs);
                            otherwise
                                atitle=sprintf('Found %d new equilibria (%d)',n,neqs);
                        end
                    end
                end
                i_waitbar(1,atitle);
                if ~isempty(hwait)&&ishandle(hwait)&&get(hwait,'userdata')==1
                    disp('Cancelled by user');
                    g_grind.findeqs.maxtimereached=true;
                    break;
                end
            end
            i_waitbar([]);
        end
    end
    if ntrials > 0
        g_grind.findeqs.eqs = eqs;
        g_grind.findeqs.settings = settings;
    end
    i_keep(oldN0);
catch err
    %   err=lasterror;
    i_keep(oldN0);
    rethrow(err);
end
[eqlist1,eqs] = update_eqlst;
if nargout ==0
    findeqs('-l')
    %     % fprintf('Found %d equilibria (might not be all!):\n', size(eqlst, 1));
    %     [Ans, ok] = listdlg( 'ListString', eqlist1 , ...
    %         'SelectionMode' ,'multiple' ,'ListSize' , [400 300],...
    %         'Name','findeqs',       'PromptString'  ,...
    %         sprintf('Found %d equilibria (might not be all!), select one:', length(eqlist1)));
    %     if ok
    %         for i = 1:length(Ans)
    %             N1= g_grind.findeqs.eqs(Ans(i),:)';
    %             i_keep(N1);
    %             N1=findeq('maxiter',1); %#ok<NASGU>
    %         end
    %     end
else
    eqlist=eqlist1;
end
function [eqlst,eqs]=update_eqlst(toleq)
global g_grind
if ~isfield(g_grind,'findeqs')||~isfield(g_grind.findeqs,'eqs')||isempty(g_grind.findeqs.eqs)
    eqlst=struct('id',{},'label',{},'msg','','data',[],'N0',[]);
    eqs=[];
    return;
end
if nargin<1
    toleq=g_grind.findeqs.toleq;
end
eqs=clusteqs(g_grind.findeqs.eqs,toleq,g_grind.solver.iscomplex);
if ~isempty(eqs)
    [~, index] = sort(min(eqs, [], 2) + sum(eqs, 2) / 1000);
    eqs = eqs(flipud(index), :); %sort in an order that the trivial / negative equilibria are last
end
eqlst = struct('id','EP','label','EP','msg','','data',[],'N0', num2cell(transpose(eqs),1));
for i = 1:size(eqs, 1)
    N0 = eqlst(i).N0;
    [lab,eigenval]=i_get_equil_label(N0);
    eqlst(i).msg = lab;
    eqlst(i).data.eigenval=eigenval;
    eqlst(i).id=['EP',num2str(i)];
end


function [eqs] = addeq(eqs, N0, findeq_args)
findeq_args.n0=N0;
[N1, found] = findeq(findeq_args,'display','off');
if found
    eqs(end + 1, :) = transpose(N1(:));
end

%se(j + 1) = ntrial;
%settings(j + 2) = scal;
%j=j+2

function eqs=clusteqs(eqs,toleq,iscomplex)
%add a random number to each column to see the difference between [0 10 0] and [0 0 10]
if isempty(eqs)
    return;
end
if iscomplex
    eqs=[real(eqs),imag(eqs)];
end
eqs1=bsxfun(@times,eqs,1:size(eqs,2));
%determine the distance with the first element
eqs1(eqs1<-0.1)=eqs1(eqs1<-0.1)*10; %some models are symmetric in zero
dist=sum((bsxfun(@minus,eqs1,rand(1,size(eqs1,2)))).^2,2);
[~,ndx]=sort(dist,'descend');
eqs1=eqs(ndx,:);
%distance between subsequent points
d=sum((eqs1(1:end-1,:)-eqs1(2:end,:)).^2,2);
if length(d)<3
    crit=toleq;
else
    d1=sort(d,'descend');
    i=2;
    while var(d1(i:end))>toleq
        i=i+1;
    end
    crit=d1(i)*1.5;
end
borders=[1; find(d>crit)+1;size(eqs1,1)+1];
eqs=zeros(length(borders)-1,size(eqs,2));
for i=2:length(borders)
    eqs(i-1,:)=mean(eqs1(borders(i-1):borders(i)-1,:),1);
end
if iscomplex
    eqs=complex(eqs(:,1:end/2),eqs(:,end/2+1:end));
end



