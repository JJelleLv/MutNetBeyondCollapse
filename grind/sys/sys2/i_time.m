function [outmat, leg] = i_time(varargin)
global g_Y g_t t g_grind g_data g_paranal;
if isfield(g_grind,'timevars')&&~isempty(g_grind.timevars)
    s=sprintf('%s ',g_grind.timevars{1}{:});
else
    s='';
end
fieldnams={'ndays', 'n>0', 'runs the model for NDAYS time units.', g_grind.ndays;...
   'timevars', 's', 'list of equations, use in combination with ''-out''',s}';
res = i_time_options(i_parseargs(fieldnams,'if hasoption(''-o''),deffields=''timevars(+)'';else,deffields=''ndays'';end;',...
    {'-p','(-[1-9]+[0-9]*)','-h','-f','-o','-r','-a','-n','-s','-p2|-paranal2','-p1|-pa'},varargin));
i_parcheck;
hol = 0;
filename='';
if any(strcmp(res.opts,'-o'))
    %out1 = g_grind.timevars;
    if isfield(res,'timevars')
        v = res.timevars;
    else
        v=[];
    end
    numopts=regexp(res.opts,'[\-][1-9]+[0-9]*','once','match');
    ndx=cellfun('isempty',numopts);
    v=[numopts(~ndx) v];
    if isempty(v)
        uiwait(i_outdlg);
    else
        out(v{:});
        %return;
    end

end

if any(strcmp(res.opts,'-h'))
    hol = 1;
end

f1=strncmpi(res.opts, '-f', 2);
if any(f1)
    f=strfind(res.opts{f1}, '=');
    if ~isempty(f)
        filename = res.opts{f1}(f(1) + 1:end);
        pathname = pwd;
    else
        [filename, pathname] = uiputfile('*.txt', 'Save As');
        if filename == 0
            filename = [];
        end

    end

end

if any(strcmp(res.opts,'-p'))
    if isempty(g_paranal.run)
        error('grind:time','Paranal should be run before using this option');
    end

    if isempty(g_Y)
        g_Y=g_paranal.run.Y(:,:,1);
        g_t=g_paranal.run.t(:,:,1);
    end

    time('-n');
    addreplaydata;
    return;
end


if res.rerun
    i_ru(t, res.ndays, res.N0, 1);
end

hastimevars = ~isempty(g_grind.timevars);
if hastimevars
    hastimevars = 0;
    for i = 1:length(g_grind.timevars)
        if ~isempty(g_grind.timevars{i})
            hastimevars = 1;
            break;
        end

    end

end

if ~hastimevars
    i_warningdlg('GRIND:time:novars','No time variables defined, use <a href="matlab:out">out</a> to define output');
end

if (nargout > 0) || ~isempty(filename)
    outmat1 = cell(1,size(g_grind.timevars, 2));
    leg = cell(1,size(g_grind.timevars, 2));
    for No = 1:size(g_grind.timevars, 2)
        if ~isempty(g_grind.timevars{No})
            leg{No} = g_grind.outt(No);
            outdata = outfun(g_grind.outt{No});
            outdata2 = [];
            for i = 1:length(g_grind.timevars{No})
                avar=char(g_grind.timevars{No}{i});
                leg{No}{i + 1} = i_disptext(avar);
                aa = outfun(avar);
                if ~isempty(g_data) && strncmpi(avar, 'Observed ',9)
                    outdata2 = [outdata2, aa];
                else
                    outdata = [outdata, aa];
                end

            end

            if ~isempty(outdata2)
                if ~strcmp(g_grind.outt{No},'t')
                    error('grind:time','cannot get "%s" for observed data',g_grind.outt{No});
                end

                if ~isempty(outdata)
                    [tr, datar]=i_concatdata(outdata(:,1),outdata(:,2:end),g_data.t,outdata2);
                    outdata=[tr,datar];
                else
                    outdata=[g_data.t,outdata2];
                end

            end

            outmat1{No} = outdata;
        end
    end

    
    if ~isempty(filename)
        if ~isempty(pathname)
            if pathname(end) ~= filesep
                pathname = [pathname, filesep];
            end

            filename = [pathname, filename];
        end
        fid = fopen(filename , 'wt');
        for No = 1:size(g_grind.timevars, 2)
            if isempty(leg{No})
               fprintf(fid, '\t');
            else
               fprintf(fid, '%s\t', leg{No}{:});
            end

            fprintf(fid, '\n');
            for i = 1:size(outmat1{No}, 1)
                s = sprintf('%g\t', outmat1{No}(i,:));
                s = s(1:end - 1);
                s=strrep(s,'NaN','');
                fprintf(fid, '%s\n',s);
            end

            fprintf(fid, '\n');
            fprintf(fid, '\n');
        end

        fclose(fid);
        fprintf('Data written to %s\n', filename);
    end

    if length(outmat1) == 1
        outmat1 = outmat1{1};
    end

    if nargout > 0
        outmat = outmat1;
    end

end

if ~res.silent
    i_timeplot(hol);
end

if res.adding
    g_grind.solver.addmode = false;
end



%if ~isempty(out1)
%   g_grind.timevars = out1;
%end

if ~isempty(res.OldY)
    g_Y = res.OldY;
    g_t = res.Oldt;
end


function addreplaydata
global g_paranal g_grind;
if isempty(g_paranal.run)
    error('grind:time:paranal','<a href="matlab:paranal">Paranal</a> should be run before using this option');
end

for No = 1:size(g_grind.timevars, 2)
    if ~isempty(g_grind.timevars{No})
        [h] = i_makefig('time', No - 1);
        hax=findobj(h,'type','axes');
        tags = get(hax, 'tag');
        hax=hax( ~(strcmp(tags,'Colorbar')|strcmp(tags,'legend')));
        ud = get(hax, 'userdata');
        ud.replay.callback = @replaycallback;
        ud.replay.onstart = @replaystart;
        ud.replay.onend = @replayend;
        ud.replay.onturn = @i_replayparanalturn;
        ud.replay.pars = g_paranal.run.parvalues;
        tstep = (ud.replay.pars(end) - ud.replay.pars(1)) / length(ud.replay.pars);
        ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',length(ud.replay.pars)*10,'ndata',numel(g_paranal.run.t));
        ud.replay.opt = No;
        set(hax,'userdata', ud);
    end

end

replayall('value','start','variable', g_paranal.run.pars{1});

function replaystart(hax)
global g_grind g_paranal;
if ishandle(hax)
    i_figure(get(hax, 'parent'));
    % set(hax,'ylimmode','manual');
    % set(hax,'zlimmode','manual');
    ud = get(hax, 'userdata');
    if isfield(ud, 'replay')
        ud.replay.pars = g_paranal.run.parvalues;
        tstep = (ud.replay.pars(end) - ud.replay.pars(1)) / length(ud.replay.pars);
        ud.replay.settings=struct('tvar',g_paranal.run.pars{1},'tlim',[g_paranal.run.parvalues(1,1) g_paranal.run.parvalues(end,1)+tstep-1e-9],'tstep',tstep,'numt',length(ud.replay.pars)*10,'ndata',numel(g_paranal.run.t));
        No = ud.replay.opt;
        ser = get(hax, 'children');
        ud.replay.ydata = cell(size(ser));
        ud.replay.xdata = cell(size(ser));
        if No > length(g_grind.outt)
            No = length(g_grind.outt);
        end

        tvar = g_grind.outt{No};
        ud.replay.tvar=tvar;
        k = size(g_grind.timevars{No}, 2);
        maxy = -inf;
        mint=inf;
        maxt=-inf;
        for i = 1:size(g_grind.timevars{No}, 2)
            yvar = g_grind.timevars{No}{i};
            if ~(strncmpi(yvar,'observed(',9)||strncmpi(yvar,'observed ',9)&&strcmp(tvar, 't'))
                tt1 = outfun(tvar, '-p');
                if any(~isreal(tt1))
                    warning('grind:time','Imaginary parts of complex values ignored');
                    tt1=real(tt1);
                end
                res = outfun(yvar, '-p');
                if any(~isreal(res))
                    warning('grind:time','Imaginary parts of complex values ignored');
                    res=real(res);
                end
                ud.replay.xdata{k} = tt1;
                ud.replay.ydata{k} = res;
                if maxt<max(tt1)
                    maxt=max(tt1);
                end

                if mint>min(tt1)
                    mint=min(tt1);
                end

                if maxy < max(res)
                    maxy = max(res);
                end

            else
                ud.replay.xdata{k}=NaN;
                ud.replay.ydata{k}=NaN;
            end

            k = k - 1;
        end

        ud.replay.p=outfun(ud.replay.settings.tvar,'-p');
    end

    set(hax, 'userdata', ud);
    %   set(hax,'drawmode','fast');
    if ~isempty(g_paranal.run.t)
        if strcmp(tvar,'t')
            mint=0;
            maxt=max(g_paranal.run.t(:,1,1))-min(g_paranal.run.t(:,1,1));
        end

        xlim(hax, [mint maxt]);
        if maxy==0
            maxy=1;
        end

        ylim(hax, sort([0 maxy]));
        set(hax,'ylimmode','manual');
        set(hax,'xlimmode','manual');
    end

end


function replayend(hax,closedlg)
try
    if closedlg
        replaycallback(hax, '', 1);
    end

catch
end

if ishandle(hax)
    i_figure(get(hax, 'parent'));
    set(hax,'xlimmode','auto');
    set(hax,'ylimmode','auto');
    %  set(hax,'zlimmode','auto');
    ud = get(hax, 'userdata');
    ud.replay.xdata = {};
    ud.replay.ydata = {};
    ud.replay.p=[];
    set(hax, 'userdata', ud);
end

function p = replaycallback(hax,avar, relt)
global g_paranal;
p = [];
if ishandle(hax)&&isempty(avar)||strcmp(avar, g_paranal.run.pars{1})
    ud = get(hax, 'userdata');
    if isfield(ud, 'replay')&& isfield(ud.replay,'p')
        ndx1 = floor(relt * (numel(g_paranal.run.t) - 1)) + 1;
        p = ud.replay.p(ndx1);
        ndx=(p == ud.replay.p);
        i=1;
        lengthxdata=length(ud.replay.xdata{i});
        while lengthxdata==1&&i<length(ud.replay.xdata)
            i=i+1;
            lengthxdata=length(ud.replay.xdata{i});
        end

        if lengthxdata ~= length(ndx)
            replaystart(hax);
            ud = get(hax, 'userdata');
        end

        %sum(ndx)
        ser = get(hax, 'children');
        if any(ndx)
            for i = 1:min(length(ser), length(ud.replay.xdata))
                if length(ndx)==length(ud.replay.xdata{i})
                    tt = ud.replay.xdata{i}(ndx);
                    ydata = ud.replay.ydata{i}(ndx);
                    if sum(isnan(ydata))>length(ydata)/10
                        set(ser(i),'Marker','.');
                    else
                        set(ser(i),'Marker','none');
                    end

                    if strcmp(ud.replay.tvar,'t')
                        tt = tt - tt(1);
                    end

                    if ~isempty(tt)
                        set(ser(i), 'XData', tt, 'YData', ydata);
                    end

                end

            end

        end

    end

end


