%EIGEN   Calculate eigenvalues of currently selected equilibrium
%   Display eigenvalues and eigenvectors of the current
%   initial conditions. The eigenvalues are plotted in
%   the complex plane.
%   Additionally the command EIGEN can do an extended analysis of the reactivity, see:
%   
%   Neubert, M. G. and Caswell, H. (1997), Alternatives to resilience for measuring the responses of ecological systems to
%   perturbations. Ecology, 78: 653-665.
%
%   Usage:
%   EIGEN - Calculate eigenvalues with analytical Jacobian
%   (see <a href="matlab:help enterjac">enterjac</a>). If there is no Jacobian entered, the Jacobian
%   matrix is approximated numerically.
%   Eigs=EIGEN - calculates all eigenvalues during the last run.
%   Eigs=EIGEN('-paranal') - calculates all eigenvalues from the last paranal run.
%   EIGEN('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'nangles' [integer>=0] - Number of angles for simulation of reactivity
%     'ndays' [number] - Number of time steps for simulating of reactivity (NaN=automatic)
%     'npoints' [integer>0] - Number of steps for simulating of reactivity
%     'sizedist' [number>0] - Size of the perturbation of perturbation
%     'statevars' [state variable and length(state variable)==2] - Which two variables should be used for the perturbations of reactivity
%   EIGEN('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-a' - analytical: Jacobian matrix is calculated numerically, only if an analytical Jacobian is entered.
%     '-n' - numerical: Jacobian matrix is approximated numerically, even if an analytical Jacobian is entered.
%     '-p' - Eigs=EIGEN('-paranal') - calculates all eigenvalues from the last paranal run.
%     '-r' - -reactivity calculates also the reactivity (=max eigenvalue of the Hermitian matrix (=(Jacobian+Jacobian')/2)).
%
%
%   See also findeq, enterjac, paranal, trdet
%
%   Reference page in Help browser:
%      <a href="matlab:commands('eigen')">commands eigen</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function Eigs = eigen(varargin)
global g_grind g_Y t g_t;
fieldnams={'sizedist', 'n>0', 'Size of the perturbation of perturbation',0.00001;...
    'nangles', 'i>=0', 'Number of angles for simulation of reactivity',100;...
    'npoints', 'i>0', 'Number of steps for simulating of reactivity',200;...
    'ndays', 'n', 'Number of time steps for simulating of reactivity (NaN=automatic)',NaN;...
    'statevars', 'v&length(v)==2', 'Which two variables should be used for the perturbations of reactivity',{g_grind.xaxis.var g_grind.yaxis.var}}';
args=i_parseargs(fieldnams,'sizedist,statevars','-n,-a,-r,-p',varargin);
i_parcheck;
doparanal = any(strcmp(args.opts, '-p'));
doreactivity = any(strcmp(args.opts, '-r'));
donumerical= isempty(g_grind.syms.Jacobian);
if any(strcmp(args.opts, '-n'))
    donumerical = true;
end
if any(strcmp(args.opts, '-a'))
    donumerical = false;
end
niters = g_grind.solver.iters;
% ud = get(get(0,'currentfigure'), 'userdata');
% if isfield(ud, 'iters')
%     niters = ud.iters;
% end
if nargout > 0
    if ~doparanal
        N0 = i_initvar;
        if isempty(g_Y)
            i_ru(t, g_grind.ndays, N0, 1);
        end
        if g_grind.solver.haslags
            Eigs = repmat(zeros(size(g_Y)), 1, 10);
        else
            Eigs = zeros(size(g_Y));
        end
        for i = 1:size(g_Y, 1)
            N0 = transpose(g_Y(i, :));
            [~, eigenvalues] = i_eigen(donumerical, niters, N0);
            Eigs(i, :) = transpose(eigenvalues);
        end
    else
        Eigs=outfun('eigen()','-p');
    end
    return;
elseif donumerical
    disp('Jacobian numerically approximated');
else
    disp('Uses Jacobian entered by user');
end

N0 = i_initvar;
[Jacobian, eigenvalues, eigenvect] = i_eigen(donumerical,niters,N0);
if g_grind.solver.haslags
    disp('This model has infinite number of eigenvalues, only the first 10 are shown');
end
disp('Jacobian:'); disp(Jacobian);
%eigenvalues = diag(eigenvalues);
disp('Eigen vectors:'); disp(eigenvect);
disp('Eigenvalues:'); disp(eigenvalues);
if doreactivity
    % disp('Reactivity:');
    if g_grind.statevars.dim==1
        warning('grind:reactivity:oneD','1D systems are never reactive: the reactivity is always the same as the eigenvalue');
    end
    %Hermitian matrix:
    H = (Jacobian + ctranspose(Jacobian)) / 2;  % should be Conjugate transpose!
    %    eigH = eig(H);
    [vectH,eigH] = eig(H);
    eigH=diag(eigH);
    disp('Eigen vectors Hermitian:'); disp(vectH);
    disp('Eigenvalues Hermitian:'); disp(eigH); %seems to have no
    %relation with the most sensitive direction
    isstable=i_stability(eigenvalues,g_grind.solver.isdiffer,g_grind.solver.haslags);
    fprintf('Total reactivity = %g\n',max(eigH))
    % if isstable
    %we can now calculate the maximum amplification (Neubert and
    %Caswell, 1997)
    %initial guess
    if max(eigH)>0&&isstable
        tt=1:0.1:100;
        rho=zeros(size(tt));
        for i=1:length(tt)
            rho(i)=norm(expm(Jacobian*tt(i)));
        end
        tmax=tt(rho==max(rho));
        tmax=tmax(1);
        han=@(t)(norm(expm(Jacobian*(t+1E-6)))-norm(expm(Jacobian*(t-1E-6))))/2E-6;
        tmax1=fzero(han,tmax);
        if ~isnan(tmax1)
            tmax=tmax1;
        end
        fprintf('Time of maximum amplification = %g\n',tmax);
        fprintf('Max amplification = %g\n',norm(expm(Jacobian*tmax)));
    else
        tmax=30;
    end
    if g_grind.statevars.dim>1
        %parameters 
        if ~isfield(args,'sizedist')
            args.sizedist=0.00001;
        end
        if ~isfield(args,'nangles')
            args.nangles=100;
        end
        if ~isfield(args,'npoints')
            args.npoints=200;
        end
        if ~isfield(args,'ndays')
            args.ndays=NaN; %number of days to simulate (if NaN is as determined by tmax or 30);
        end
        if ~isfield(args,'statevars')
            args.statevars={g_grind.xaxis.var g_grind.yaxis.var};
        end
        oldsolver=[];
        varnos=[i_varno(args.statevars{1}) i_varno(args.statevars{2})];
        if any(strcmp(solver('name'), {'c.ode45','ode45'}))
            if g_grind.solver.opt.AbsTol > 1e-10 || g_grind.solver.opt.RelTol > 1e-10
                % disp('timesens has set solver tolerances to values 1E-11');
                oldsolver.name=solver('name');
                oldsolver.RelTol=solver('RelTol');
                oldsolver.AbsTol=solver('AbsTol');
                solver(oldsolver.name, 1E-11, 1E-11);
            end
        end
        if args.nangles>0
            %get 100 random perturbations in 100 directions, by taking normal
            %distributed points and normalize them (known method)
            %circ=randn(100,g_grind.statevars.dim);
            %dcirc=sqrt(sum(circ.^2,2));
            %circ=circ./dcirc(:,ones(1,g_grind.statevars.dim));
            theta=linspace(0,2*pi,args.nangles).';
            circ=zeros(args.nangles,g_grind.statevars.dim);
            circ(:,varnos(1))=cos(theta);
            circ(:,varnos(2))=sin(theta);
            circ=circ*args.sizedist;
            try
                close(i_figno('eigen')+1);
            catch err
                if ~strcmp(err.identifier,'MATLAB:close:InvalidFigureHandle')
                    rethrow(err);
                end
            end
            figure(i_figno('eigen')+1)
            i_plotdefaults;
            xlabel(i_statevars_names(varnos(1)));
            ylabel(i_statevars_names(varnos(2)));
            hold on
            N0=i_initvar;
            eqn=N0';
            eqn=eqn(ones(1,args.npoints+1),:);
            if isnan(args.ndays)
                args.ndays=5*tmax;
            end
            oldndays=g_grind.ndays;
            oldtstep=g_grind.tstep;
            simtime(0,args.ndays,args.npoints);
            dists=zeros(args.npoints+1,args.nangles);
            maxdist=zeros(args.nangles,1);
            initdist=zeros(args.nangles,1);
            Ys=zeros(args.npoints+1,g_grind.statevars.dim,args.nangles);
            g_grind.solver.opt.OutputFcn = [];
            for i=1:args.nangles
                i_keep(N0+circ(i,:)');
                time('-s');
                if size(g_Y,1)<args.npoints+1
                    g_Y(end+1:args.npoints+1,:)=NaN;
                    g_t=linspace(0,g_grind.ndays,args.npoints+1);
                end
                Ys(:,:,i)=g_Y;
                plot(g_Y(:,varnos(1)),g_Y(:,varnos(2)));
                dists(:,i)=sqrt(sum((g_Y-eqn).^2,2));
                maxdist(i)=max(dists(:,i));
                initdist(i)=dists(2,i);
                if maxdist(i)<=1
                    amax=[diff(sign(diff(dists(:,i).^2)))<0; false];
                    if any(amax)
                        maxdist(i)=max(dists(amax,i));
                    end
                end
            end
            if ~isempty(oldsolver)
                solver(oldsolver);
            end
            %   axis tight
            hax=gca;
            [~,Hndx]=sort(real(eigH));
            addvector(hax,vectH(:,Hndx(1)),N0,varnos,'r:','1st reactivity direction')
            addvector(hax,vectH(:,Hndx(2)),N0,varnos,'b:','2nd reactivity direction')
            [~,Endx]=sort(real(eigenvalues));
            addvector(hax,eigenvect(:,Endx(1)),N0,varnos,'r-','1st eigen direction')
            addvector(hax,eigenvect(:,Endx(2)),N0,varnos,'b-','2nd eigen direction')
            hold off
            figure(i_figno('eigen')+2)
            plot(g_t,dists/args.sizedist);
            i_plotdefaults;
            xlabel('Time (t)')
            ylabel('Distance from equilibrium')
            figure(i_figno('eigen')+3);
            if exist('polarplot','file')
                polarplot(theta,ones(size(theta)),'b-')
                hold on
                polarplot(theta,maxdist/args.sizedist,'r-')
                polarplot(theta,initdist/args.sizedist,'g-');
                set(gca,'ThetaAxisUnits','radians')
                set(gca,'ThetaTick',[0 pi/2 pi 3/2*pi])
                rtick= get(gca,'Rtick');
                set(gca,'Rtick',rtick(2:end));
                i_plotdefaults
                hold off
            else
                %if polarplot doesn't exist we can draw our own polar plot
                %(the old function polar is not good)
                plot(cos(theta),sin(theta),'b-');
                i_plotdefaults;
                xlabel(sprintf('Perturbation size in %s',i_statevars_names(varnos(1))));
                ylabel(sprintf('Perturbation size in %s',i_statevars_names(varnos(2))));
                hold on;
                plot(cos(theta).*initdist/args.sizedist,sin(theta).*initdist/args.sizedist,'g-');
                plot(cos(theta).*maxdist/args.sizedist,sin(theta).*maxdist/args.sizedist,'r-')
                %                 hax=gca;
                %                 addvector(hax,vectH(:,Hndx(1)),zeros(size(N0)),varnos,'r:')
                %                 addvector(hax,vectH(:,Hndx(2)),zeros(size(N0)),varnos,'b:')
                %                 addvector(hax,eigenvect(:,Endx(1)),zeros(size(N0)),varnos,'r-')
                %                 addvector(hax,eigenvect(:,Endx(2)),zeros(size(N0)),varnos,'b-')
                hold off
            end
            
            [~,ndx]=max(maxdist);
            fprintf('Direction of max reactivity:\n    %g\n    %g\n',cos(theta(ndx)),sin(theta(ndx)));
            [~,ndx]=max(initdist);
            fprintf('Direction of initial amplification:\n    %g\n    %g\n',cos(theta(ndx)),sin(theta(ndx)));
            g_grind.ndays=oldndays;
            g_grind.tstep=oldtstep;
            i_keep(N0);
        end
    end
    %end
end
%eigenvect
%eigenvalues
%N0 = i_initvar;
Nres1 = i_runsinglestep(1, transpose(N0), true);
if g_grind.solver.isdiffer
    Nres1 = Nres1 - transpose(N0);
end
hfig = i_makefig('eigen');
plot(real(eigenvalues),imag(eigenvalues),'rx','MarkerSize',12);
ud = get(gca, 'userdata');
ud.meta=struct('func','eigen','xname','real(eig)','yname','imag(eig)','zname','');
set(gca, 'userdata', ud);

i_plotdefaults(hfig)
if sum(Nres1) > 0.001
    warning('GRIND:eigen:noequilibrium', 'The model is not in an equilibrium\nUse <a href="matlab:findeq">findeq</a> to find the equilibrium');
end
daspect([1 1 1]);
oldhold = ishold;
hold on;
if g_grind.solver.isdiffer
    x = -pi:.1:pi + 0.1;
    plot(sin(x), cos(x), 'k');
end
set(hfig,'Name','Plot of eigenvalues in the complex plane');
xlabel('Real part of eigenvalue');
ylabel('Imaginary part of eigenvalue');
H1 = get(hfig, 'CurrentAxes');
if g_grind.solver.isdiffer
    max1=max(abs([get(H1, 'Xlim') get(H1, 'Ylim')]));
    lim = [ - max1, max1];
else
    max1 = max(abs([real(eigenvalues); imag(eigenvalues)]));
    if max1 < 0.18
        lim = [ - 0.2 0.2];
    elseif max1 < 0.35
        lim = [ - 0.4 0.4];
    elseif max1 < 0.45
        lim = [ - 0.5 0.5];
    elseif max1 < 0.95
        lim = [ - 1 1];
    elseif max1 < 1.4
        lim = [ - 1.5, 1.5];
    elseif max1 < 10
        lim = [ - (round(max1) + 1), round(max1) + 1];
    else
        max1=max(abs([get(H1, 'Xlim') get(H1, 'Ylim')]));
        lim = [ - max1, max1];
    end
end
set(H1, 'Xlim', lim);
set(H1, 'Ylim', lim);
plot([complex(0, lim(1)), complex(0, lim(2))], 'k');
plot([complex(lim(1), 0), complex(lim(2), 0)], 'k');
if ~oldhold
    hold off;
end
if g_grind.solver.nonautonomous
    error('GRIND:eigen:NonAutonomous','Eigenvalues are probably not correct for this non-autonomous model. Consider to make equation autonomous by adding tau''=1 and replacing all t by tau');
end
function addvector(hax,vectmax,N0,varnos,color,legendtext)
xrange=xlim(hax);
yrange=ylim(hax);
a=vectmax(varnos(2))/vectmax(varnos(1));
if ~isreal(a)
    warning('grind:eigen:complex','Cannot plot complex eigenvectors');
    return;
end
plot(xrange,N0(varnos(2))+(xrange-N0(varnos(1)))*a,color,'DisplayName',legendtext);
xlim(hax,xrange);
ylim(hax,yrange);
