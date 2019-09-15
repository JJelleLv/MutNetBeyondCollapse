function res = i_time_options(args)
%-add
%-no check
%-p1
%-p2
%-r
%-s
% number of days
global g_grind ;
res=args;
res.OldY = [];
res.Oldt = [];
%res.opt=args.opts;
res.nocheck = any(strcmp(args.opts,'-n'));
res.adding = any(strcmp(args.opts,'-a'));
if res.adding
    g_grind.solver.addmode = true;
end

res.silent = any(strcmp(args.opts,'-s'));
res.rerun = any(strcmp(args.opts,'-r'));
if res.rerun
    timesens(struct('opts','-u'));
end

%out1 = [];
if isfield(args,'ndays')&&~isempty(args.ndays)
    res.ndays = args.ndays;
else
    res.ndays= g_grind.ndays;
end

if any(strcmp(args.opts,'-p1'))
    [res.OldY, res.Oldt] = setparanal(1);
    res.nocheck = true;
end

if any(strcmp(args.opts,'-p2'))
    [res.OldY, res.Oldt] = setparanal(2);
    res.nocheck = true;
end

res.N0 = i_initvar;
if res.rerun || (~res.nocheck && i_settingschanged(res.N0,i_checkstr(res.ndays)))
 %   disp('running');
    res.rerun=true;
end



function [OldY, Oldt] = setparanal(iparanal)
global g_Y g_t g_paranal;
if iparanal > 0
    OldY = g_Y;
    Oldt = g_t;
    g_Y = g_paranal.run.Y;
    g_t = g_paranal.run.t;
    if iparanal == 2&&isfield(g_paranal,'prevrun')&&~isempty(g_paranal.prevrun)
        g_Y = cat(3,g_paranal.prevrun.Y, g_Y);
        g_t = cat(3,g_t, g_paranal.prevrun.t + g_t(end));
    end

    g_Y = transpose(reshape(permute(g_Y, [2, 1, 3]), size(g_Y, 2), size(g_Y, 1) * size(g_Y, 3)));
    g_t = g_t(:);
else
    OldY = [];
    Oldt = [];
end

