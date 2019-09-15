function [x0,v0]= init_P_EP(odefile, x, p, ap, ndays)
global g_Y;
if nargin<5
    ndays=1000;
end
N0=x;
[N0,isequil]=findeq('display','off','maxiter',1,'tolfun',1E-3,'n0',N0,'keep',false);
%if you start in equilibrium don't do anything
if ~isequil
    stabil('ndays',ndays,'keep',false);
    N0_1=g_Y(end,:)';
    [N0_1,iseq]=findeq('display','off','maxiter',1,'tolfun',1E-3,'n0',N0_1,'keep',false);
    if iseq
        N0=N0_1;
        fprintf('Initial equilibrium found by simulation:\n%s\n',i_get_equil_label(N0))
    else
        stabil('backwards','true','ndays',ndays,'keep',false);
        N0_2=g_Y(1,:)';
        i=1;
        while isnan(N0_2(1))
            i=i+1;
            N0_2=g_Y(i,:)';
        end
        [N0_2,isequil]=findeq('display','off','maxiter',1,'tolfun',1E-3,'n0',N0_2,'keep',false);
        if isequil
            N0=N0_2;
            fprintf('Initial equilibrium found by simulation:\n%s\n',i_get_equil_label(N0))
        else
            fprintf('Failed to find equilibrium:\n%s\n',i_get_equil_label(N0))
        end
    end
else
    fprintf('Initial value in or close to equilibrium:\n%s\n',i_get_equil_label(N0))
end
[x0,v0]= init_EP_EP(odefile, N0, p, ap);
