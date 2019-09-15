%EXTERNVAR   Get the value of an external variable
%   This function can evaluate an external variable. Use <a href="matlab:help defextern">defextern</a> to 
%   define external variables. The function interpolates the external (file) data. 
%   The used method is linear interpolation. 
%
%
%   Usage:
%   EXTERNVAR(varnr,default,t) - the function has three parameters: (1) variable number (see
%   g_grind.externalvars. 2) default, a value which is used outside the range of the data
%   (3) t, current time or a vector with time.
%   EXTERNVAR(varnr,t) - if no default value is given, the function returns NaN outside
%   the range of the data.
%   EXTERNVAR(name,t) - You can also use the name of the variable.
%
%
%   See also defextern, setdata, loaddata, model, rednoise, externlag
%
%   Reference page in Help browser:
%      <a href="matlab:commands('externvar')">commands externvar</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function V = externvar(varnr, default, at)
global g_grind g_t;
if nargin == 2
    at = default;
    default = NaN;
elseif (nargin==1)&& ischar(varnr)
    i=1;
    while i<=length(g_grind.externvars)&&~strcmp(varnr,g_grind.externvars{i}.name)
        i=i+1;
    end
    varnr=i;
    if i<=length(g_grind.externvars)
        V=externvar(varnr, str2num(g_grind.externvars{varnr}.default),g_t); %#ok<ST2NM>
        return;
    else
        error('GRIND:externvar:UnknownName','Unknown name');
    end
elseif nargin==0
    error('GRIND:externvar:TooFewArgs','EXTERNVAR: too few arguments of the function');
end
options=g_grind.externvars{varnr}.options;
data=g_grind.externvars{varnr}.data;
if ~options.active
    V= evalin('base',g_grind.externvars{varnr}.name);
    if numel(at)>1
        V=transpose(V(:));
        V=repmat(V,numel(at),1);
    end
    return;
end
s = size(data);
% if ((s(1) == 2) || (s(1) == 1)) && (s(2) > 2)
%    data = data';
%    s = size(data);
% end
if options.tofloor
    at=floor(at);
end
if s(2) == 1
    if options.cycle
        at = mod(at, s(1));
    end
    if length(at) > 1
        V = interp1(transpose(1:s(1)), data, at);
    else
        if (at + 1 > s(1)) || (at < 0)
            V = NaN;
        else
            remain = rem(at + 1, 1);
            V = data(floor(at + 1)) * (1 - remain) + data(ceil(at + 1)) * remain;
        end
    end
elseif s(2) >= 2&& s(1)>0
    if options.cycle
        t1=data(1,1); tend=data(end,1);
        at = mod(at, tend - t1)+t1;
    end
    if isempty(at)
        V=[];
    elseif length(at) > 1
        if size(at,2)>1
            V = interp1(data(:, 1), data(:, 2:end), transpose(at));
        else
            V = interp1(data(:, 1), data(:, 2:end), at);
        end
    else
        V = parlookup(data, at,0,0);
    end
elseif s(1) == 0
    V = ones(length(at),numel(default)) * NaN;
else
    error('GRIND:externvar:dataerror','EXTERNVAR: data matrix not correct')
end
if (length(at)==1)
    if numel(V)==1&&isnan(V)
        V=default;
    else
        ndx=isnan(V);
        V(ndx) = default(ndx);
        V=reshape(V,g_grind.externvars{varnr}.dim1,g_grind.externvars{varnr}.dim2);
    end
else
    ndx=find(all(isnan(V),2));
    V(ndx,:) = repmat(transpose(default(:)),length(ndx),1);
end
if size(at,2)>0&&strcmp(g_grind.solver.opt.Vectorized,'on')
    V=transpose(V);
end
