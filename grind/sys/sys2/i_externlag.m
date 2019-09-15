%EXTERNLAG   External variable with time lag.
%
%   Reference page in Help browser:
%      <a href="matlab:commands('lag')">commands lag</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 09-Apr-2019 22:50:09 $
function res=i_externlag(var,timelag,t)
global g_grind;
i=1;
while i<=length(g_grind.externvars)
    if strcmp(var,g_grind.externvars{i}.name)
        res=externvar(i,str2num(g_grind.externvars{i}.default),t-timelag); %#ok<ST2NM>
        return
    end;
    i=i+1;
end
res=evalin('base',var);

