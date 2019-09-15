% flexible legend (including settings)
%    flag= -1: %remove item
%    flag=0: add text, do not save settings
%    flag=1: %add item, save settings
%    flag=2: %add item, save settings (only paraeters)
%    flag= 10: %update legend (indicate if settings changed)
%    flag= 11: %update and show legend (indicate if settings changed)
%    flag= 12: %remove settings from legend
function pout=i_grindlegend(flag, h, par2, xtrapar)
if nargin==0
    flag=11;
end

switch flag
    case  -1 %remove item
        if ishandle(h)
            for i=1:numel(h)
                set(get(get(h(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude series in legend
            end

        end

    case 0 %add item, do not save settings
        if ishandle(h)
            set(h, 'displayname',i_disptext(par2));
            set(h, 'userdata',struct('nr',-1))
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); % Include series in legend
        end

    case 1 %add item, save settings
        if nargin<4
            xtrapar={};
        end

        addleg(h,par2,1,xtrapar)
    case 2 %update legend, save settings (only parameters)
        if nargin<4
            xtrapar={};
        end

        addleg(h,par2,0,xtrapar)
    case 3 %get the series in the legend
        if nargin < 2
            h = get(get(0,'currentfigure'),'currentaxes');
        end

        pout = getseries(h);
    case 4 %are parameters changed?
        pout=1;
        if nargin < 2
            h = get(get(0,'currentfigure'),'currentaxes');
        end

        if ishandle(h)
            ud = get(h, 'userdata');
            if ~isempty(ud)&&isfield(ud,'settings')
                pout=~struccmp(ud.settings{length(ud.settings)}, par('-v',0));
            end

        end

    case 10 %update legend
        if nargin < 2
            h = get(get(0,'currentfigure'),'currentaxes');
        end

        if ishandle(h)
            ud = get(h, 'userdata');
            if ~isempty(ud)&&isfield(ud,'settings')
                s1 = cell(size(ud.settings));
            else
                return;
            end

            dfields = {};
            for k=1:length(ud.settings)-1
                for i = k+1:length(ud.settings)%can be improved
                    [cmp, d] = struccmp(ud.settings{k}, ud.settings{i});
                    if ~cmp
                        if size(d,1)>1
                            d=d.';
                        end
                        dfields = [dfields d];
                    end

                end

            end

            dfields = unique(dfields);
            if length(dfields) <= 3
                for i = 1:length(ud.settings)
                    s1{i} = '';
                    for j = 1:length(dfields)
                        if isfield(ud.settings{i},dfields{j})
                            v = ud.settings{i}.(dfields{j});
                            if strcmp(dfields{j},'g_increasingpar')
                                if v==1
                                    s1{i}=sprintf('forward %s', s1{i});
                                elseif v==2
                                    s1{i}=sprintf('unstable node%s', s1{i});
                                elseif v==3
                                    s1{i}=sprintf('stable node%s', s1{i});
                                elseif v==4
                                    s1{i}=sprintf('saddle %s', s1{i});
                                elseif v==5
                                    s1{i}=sprintf('Hopf %s', s1{i});
                                elseif v==6
                                    s1{i}=sprintf('Fold %s', s1{i});
                                elseif v==7
                                    s1{i}=sprintf('Transcritical %s', s1{i});
                                else
                                    s1{i}=sprintf('backward %s', s1{i});
                                end

                            elseif numel(v) > 1
                                s1{i}=sprintf('%s mean %s=%s', s1{i}, dfields{j}, strcompact( mean(mean(v))));
                            else
                                s1{i}=sprintf('%s %s=%s', s1{i}, dfields{j}, strcompact(v));
                            end

                        end

                    end

                    s1{i} = strtrim(s1{i});
                end

            else
                for i = 1:length(ud.settings)
                    s1{i} = sprintf('run %d', i);
                end

            end

            [ser, uds] = getseries(h);
            if length(ser) <= 1
                legend(h,'hide');
            else
                s = get(ser, 'DisplayName');
                for i = 1:length(s)
                    f = strfind(s{i}, ' [');
                    if ~isempty(f)
                        s{i} = s{i}(1:f(1) - 1);
                    end

                    if (length(uds)>=i)&&isfield(uds{i},'nr')&&(uds{i}.nr<=length(s1))
                        txt = s1{uds{i}.nr};
                    else
                        txt='';
                    end

                    if isempty(txt)
                        set(ser(i), 'DisplayName',s{i});
                    else
                        set(ser(i), 'DisplayName', sprintf('%s [%s]',s{i},  txt));
                    end

                end

            end

        end
    case 11
        if nargin < 2
            h = get(get(0,'currentfigure'),'currentaxes');
        end

        i_grindlegend(10, h);
        [ser] = getseries(h);
        if length(ser)>1
            s = get(ser, 'DisplayName');
            showlegend(h,s);
        elseif ~isempty(ser)
            s = get(ser(1), 'DisplayName');
            s1= get(get(h,'ylabel'), 'string');
            if isempty(s1)
                ylabel(s);
            end

        end
    case 12 %strip the settings
        if nargin < 2
            h = get(get(0,'currentfigure'),'currentaxes');
        end

        [ser] = getseries(h);
        if length(ser) <= 1
            legend(h,'hide');
        else
            s = get(ser, 'DisplayName');
            for i = 1:length(s)
                f = strfind(s{i}, ' [');
                if ~isempty(f)
                    s{i} = s{i}(1:f(1) - 1);
                end

                set(ser(i), 'DisplayName',s{i});
            end

            showlegend(h,s);
        end
        set(h,'userdata',[]);
        set(ser,'userdata',[]);
    case 13 % special for combined figures
        if nargin < 2
            hfig=get(0,'currentfigure');
            hs = findall(hfig,'type','axes');
        else
            hs=h;
            hfig= get(h(1),'parent');
        end

        for i=1:length(hs)
            i_grindlegend(11,hs(i));
        end

        hleg=findobj(hfig,'tag','legend');
        set(hleg,'fontsize',8);
    otherwise
        error('grind:legend','wrong flag');
end

function [ser, uds] = getseries(h)
ser = get(h, 'children');
if ~isempty(ser)
    ndx=~strcmp(get(ser,'type'),'text');
    ser=ser(ndx);
end

hasprop=isprop(ser,'UserData')&isprop(ser,'DisplayName');
ser = ser(hasprop);
uds = get(ser, 'UserData');
if ~iscell(uds)%&&hasfield(uds,'nr')
    uds = {uds};
end

hasprop = true(size(uds));
for i = 1:length(hasprop)
    hasprop(i) = hasprop(i) & ~isempty(uds{i});
end

ser = ser(hasprop);
uds = uds(hasprop);
function addleg(h,par2,addvar,xtrapar)
if ishandle(h)
    for i=1:length(h)
        set(get(get(h(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); % Include series in legend
    end
    
    hax = get(h(1), 'parent');
    ud = get(hax, 'userdata');
    if isempty(ud)||~isfield(ud, 'settings')
        ud.settings = {};
    end

    if isstruct(xtrapar)
        psett=xtrapar;
    else
        psett=par(struct('opts', {'-v'},'statevars',addvar));
        for i=1:2:length(xtrapar)
            psett.(xtrapar{i})=xtrapar{i+1};
        end

    end

    [ud.settings, no] = findsettings(ud.settings, psett);
    set(hax, 'userdata', ud);
    if iscell(par2)
        for i=1:length(h)
            set(h(i), 'userdata', struct('nr',no), 'displayname',i_disptext(par2{i}));
        end

    else
        set(h, 'userdata', struct('nr',no), 'displayname',i_disptext(par2));
    end

end

function hleg=showlegend(h,s)
if ~all(strcmp(s,' '))
    maxlen=0;
    for i=1:length(s)
        lens=length(s)-length(strfind(s,'\partial'))*7;
        if lens>maxlen
            maxlen=lens;
        end

    end

    legend(h,flipud(s),'Location','NorthEast');
    hleg=findobj(get(h,'parent'),'tag','legend');
    if (length(s)>10)||(maxlen>15)
        set(hleg,'fontsize',8);
        if length(s)>18
            legend(h,'boxoff');
        end

    else
        set(hleg,'fontsize',14);
    end

end

function [settings, no] = findsettings(settings, asett)
for i = 1:length(settings)
    if struccmp(settings{i}, asett)
        no = i;
        return;
    end

end

%add the settings if none is already there
settings = [settings {asett}];
no = length(settings);
function [s]=strcompact(d)
d=full(d);
s=sprintf('%.4G',d);
if strcontains(s,'E')
    s=strrep(s,'E+0','E+');
    s=strrep(s,'E+0','E+');
    s=strrep(s,'E+','E');
    s=strrep(s,'E-0','E-');
    s=strrep(s,'E-0','E-');
end



