%SIMPLEEVENT Event for many purposes
%
%   Usage:
%   Setevent('simpleevent',tstart,commands,deltat) - if tstart or deltat is a 
%   parameter then use string
%
%   Example:
%   setevent('simpleevent','thatch','Z=Z+Zin','yearlength');
%
%   See also:
%   SETEVENT
function next=simpleevent(t,comm,nextt)
comm=[comm ';'];
evalin('base',comm);
if nargin>=3
    next=t+i_checkstr(nextt);
else
    next=inf;
end