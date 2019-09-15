%MULTASSIGNIN Assign variable in workspace.
%    MULTASSIGNIN(WS,'name',V) assigns the variable 'name' in the
%    workspace WS the value V.  WS can be one of 'caller' or 'base'. 
%    The variable 'name' may also include vector or matrix elements, e.g. 
% 
%    See also:
%    ASSIGNIN, EVALIN
%
function multassignin(ws, name, V)
if ~strcontains(char(name),'(')
   assignin(ws, name, V);
else
   assignin(ws,'g_l_xxxx',V);
   evalin(ws, sprintf('%s=g_l_xxxx;clear(''g_l_xxxx'');',name));
end
