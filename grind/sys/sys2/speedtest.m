%test the speed of solvers (uses a max time for slow solvers)
function [speed,solverlist]=speedtest(asolver,maxt,ndays)
global g_grind g_t g_Y;
if nargin<2
   maxt=3;
end
if nargin<3
    ndays=g_grind.ndays;
end

if nargin>=0&&~isempty(asolver)
    if any(strcmp(asolver,{'-all','-noc'}))
        solverlist=solver('list');
        oldsolver=solver('name');
        if strncmp(asolver,'-noc',4)
            solverlist=solverlist(~strncmp(solverlist,'c.',2));
        end

        if nargout==0
            for i=1:length(solverlist)
                try
                    speedtest(solverlist{i},maxt,ndays);
                catch
                    fprintf('Error running %s\n',solverlist{i});
                end
            end

        else
            speed=zeros(size(solverlist))+NaN; 
            for i=1:length(solverlist)
                try
                    speed(i)= speedtest(solverlist{i},maxt,ndays);
                catch
                    fprintf('Error running %s\n',solverlist{i});
                end
            end

        end

        solver(oldsolver);
        return;
    end

    solver(asolver);
else
    asolver=g_grind.solver.name;
end
loc_odespeedtest(maxt,[],'setup');
oldtstep=g_grind.tstep;
simtime(0,g_grind.ndays,nan);
tic
g_grind.solver.opt.OutputFcn=@loc_odespeedtest;
time(ndays,'-s','-r');
g_grind.solver.opt.OutputFcn=[];
if isnan(g_Y(end,1))&&size(g_Y,1)>1
    ndx=size(g_Y,1)-1;
else
    ndx=size(g_Y,1);
end
speed1=toc/(g_t(ndx)-g_t(1));
g_grind.tstep=oldtstep;
p=solver('-properties');
if nargout>0
    speed=speed1*1000;
    solverlist={asolver};
else
    if p.hasfixedstep
       fprintf('Speed of %s (fixed step=%g) = %g ms/step\n',asolver,g_grind.solver.opt.StepSize,speed1*1000);
    else
       fprintf('Speed of %s = %g ms/step\n',asolver,speed1*1000);
    end 
end


%for a speed test set a maximum time to a simulation
%use this function as OutputFcn in solvers
%initialise the maximum time with i_odespeedtest(maxtime,[],'setup');
%
%use tic before running
function status = loc_odespeedtest(t, ~, flag)
persistent maxtime;
if strcmp(flag,'setup')
    maxtime=t;
    return;
end

status = toc>maxtime; 


