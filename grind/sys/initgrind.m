%INITGRIND   Initiate grind
%   initiate global variables. If you use <a href="matlab:help use">use</a> or <a href="matlab:help model">model</a>, you never
%   need to call this function directly. Only used in combination
%   with <a href="matlab:help setodefile">setodefile</a>.
%
%   See also use, model, setodefile
%
%   Reference page in Help browser:
%      <a href="matlab:commands('initgrind')">commands initgrind</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
%start
global t g_Y g_t g_data g_grind g_paranal g_cont;
addpath([grindpath filesep 'sys2']);
g_grind=i_init_g_grind;
if ishandle(findobj(0,'tag','i_conteq_cont_dlg'))
    delete(findobj(0,'tag','i_conteq_cont_dlg'));
end
if ~isempty(g_cont)
    g_cont.close;
end
clear global g_cont
if ~isoctave&&verLessThan('matlab','7.12')
    rand('twister',sum(100*clock)); %#ok<RAND>
elseif ~isoctave
    rng('shuffle','twister');
end
warning off backtrace
%warning('on','backtrace');
clear l_r l_asys;
t = 0;
%g_noise=[];
g_Y = [];
g_t = [];
setdefaults(1);
g_data=[];
g_paranal.run=[];

