function [x,v] = init_P_LC(odefile,x,p,ap,ntst,ncol,ndays,cycletol)
if nargin<8
    cycletol=1E-2;
end
if nargin<7
    ndays=1000;
end
i_keep(x);
stabil(ndays);
simtime(0,ndays,2000)
time('-s');
nm=nmaxima;
%if nm.cyclic&&~nm.nonperiodic
    [x,v]=initOrbLC(odefile, nm.tcycle,nm.Ycycle, p, ap, ntst, ncol,cycletol);
%else
%    error('matcont:init_P_LC','Did not find a limit cycle');
%end
