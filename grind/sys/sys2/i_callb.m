function i_callb(flag, hobject,avar, ud)
global  g_grind t;
hfig=getparentfig(hobject);
switch flag 
    case 'mdown'
        hax = get(hfig, 'CurrentAxes');
        ud=get(hax,'userdata');
        ud.initpt = get(hax, 'CurrentPoint');
        if isfield(ud,'meta')
            if ~isempty(ud.meta.zname)
                ud.initpt = i_depthslider(hfig,ud.initpt, ud.meta);
                if isempty(ud.initpt)
                    return;
                end
            else
                ud.initpt = ud.initpt(1, [1, 2]);
            end
            if ~isfield(ud.meta,'ylim')
                ylim1=get(hax,'ylim');
            else
                ylim1=ud.meta.ylim;
            end
            if ~isfield(ud.meta,'xlim')
                ud.meta.xlim=get(hax,'xlim');
            end
            if ~(ud.initpt(1)<ud.meta.xlim(1)||ud.initpt(1)>ud.meta.xlim(2)||ud.initpt(2)<ylim1(1)||ud.initpt(2)>ylim1(2))
                set(hax,'userdata',ud);
                i_mymenu(ud)
            end
        end
%     case 'mdown'
%         %  hax = get(hfig, 'CurrentAxes');
%         N0=i_initvar;
%         oldNO=N0;
%         hax = get(hfig, 'CurrentAxes');
%         atag  = get(hfig, 'tag');
%         initpt = get(hax, 'CurrentPoint');
%         XLim = get(hax, 'Xlim');
%         YLim = get(hax, 'Ylim');
%         Ylabel= get(get(hax,'ylabel'),'string');
%         Zlabel= get(get(hax,'zlabel'),'string');
%         
%         labels={get(get(hax,'xlabel'),'string'),get(get(hax,'ylabel'),'string'),get(get(hax,'zlabel'),'string')};
%         if ~isempty(labels{3})
%             initpt = i_depthslider(initpt, labels);
%             if isempty(initpt)
%                 return;
%             end
% 
%         else
%             initpt = initpt(1, [1, 2]);
%         end
% 
%         p2 = [];
%         ndx = 1:g_grind.statevars.dim;
%         if ~((initpt(1) < XLim(1)) || (initpt(1) > XLim(2)) || (initpt(2) < YLim(1)) || (initpt(2) > YLim(2)))
%             %g_Y = transpose(N0);
%             iy = i_getno(g_grind.yaxis.var);
%             if ~strcontains(char(Ylabel),'''')
%                 p2 = i_getparvar(iy);
%                 N0 = i_setparvar(iy, N0, initpt(2));
%                 ndx=ndx(ndx ~= iy.no);
%             end
%             ud=get(hax,'userdata');
%             if isfield(ud,'var')
%                 ix = i_getno(ud.var);
%             else
%                 ix = i_getno(g_grind.xaxis.var);
%             end
% 
%             if strcmp(atag, 'plotdiff')&&~ix.isvar
%                 ix = i_getno(g_grind.yaxis.var);
%             end
% 
%             p1 = i_getparvar(ix);
%             N0 = i_setparvar(ix, N0, initpt(1));
%             ndx=ndx(ndx ~= ix.no);
%             if length(initpt) > 2
%                 iz = i_getno(g_grind.zaxis.var);
%                 if ~strcontains(char(Zlabel),'''')
%                     p2 = i_getparvar(iz);
%                     N0 = i_setparvar(iz, N0, initpt(3));
%                     ndx=ndx(ndx ~= iz.no);
%                 end
%             end
% 
%             mess = cell(1, g_grind.statevars.dim + 1);
%             mess{1} = 'Mouse cursor at:';
%             for i = 1:g_grind.statevars.dim
%                 mess{i + 1}=[i_statevars_names(i) ' = ' num2str(N0(i)) ];
%             end
% 
%             i = i + 1;
%             if ix.ispar
%                 i = i + 1;
%                 mess{i}=[g_grind.pars{ix.no} ' = ' num2str(initpt(1))];
%             end
% 
%             if iy.ispar
%                 i = i + 1;
%                 mess{i}=[g_grind.pars{iy.no} ' = ' num2str(initpt(2))];
%             end
% 
%             ud = get(hfig, 'userdata');
%             if isfield(ud, 'stop')
%                 ud.stop = 0;
%                 set(hfig, 'UserData', ud);
%                 set(findobj(hfig,'Tag','stop'),'Visible','off');
%                 drawnow;
%             end
% 
%             %set(hfig, 'WindowButtonDown', '');
%             %set(hfig, 'WindowButtonMotionFcn', '');
%             %oldpoint = get(hfig, 'Pointer');
%             %set(hfig,'Pointer','arrow');
%             try
%                 i_mymenu(mess,struct('ndx',ndx,'N0',N0,'ix',ix,'iy',iy,'p1',p1,'p2',p2,'oldNO',oldNO));  %dialog box for selecting an action
%             catch
%             %    resetpoint(hfig, oldpoint);
%             end
% 
%            % resetpoint(hfig, oldpoint);
%             %       if isfield(ud, 'iters')
%             %          g_grind.solver.iters = g_grind.solver.iters * ud.iters;
%             %       end
% 
%             % try
%         end
     case 'mymenucallback'
         opt=avar;
         switch opt
             case 1
                 if isfield(g_grind, 'tilman')
                     BtnName=questdlg('Do you want to set the resource supply parameters?','Yes');
                     if strcmp(BtnName, 'Yes')
                         oldhold = ishold;
                         try
                             evalin('base', g_grind.tilman.setS);
                             hold on;
                             where('sblink','o',[0.8 0.8 0.8])
                             if ~oldhold
                                 hold off;
                             end
                             
                         catch err
                             if ~oldhold
                                 hold off;
                             end
                             disp('error in evaluating g_grind.tilman.setS');
                             rethrow(err);
                         end
                     end
                     setinitialvalues(ud);
                     ru;
                 else
                     setinitialvalues(ud);
                     ru;
                 end           
             case 2
                 setinitialvalues(ud);
                 findeq;
             case 3
                 setinitialvalues(ud);
             case 4
                 plotedit on;
                 if ishandle(hfig)
                     set(hfig,'Pointer','arrow');
                 end             
         end
         %       if isfield(ud, 'iters')
         %          g_grind.solver.iters = g_grind.solver.iters / ud.iters;
         %       end
         
         %  catch
         %     if isfield(ud, 'iters')
         %        g_grind.solver.iters = g_grind.solver.iters / ud.iters;
         %    end
         
         % end
%     case 'mymenucallback'
%         opt=avar;
%         ndx=ud.ndx;
%         N0=ud.N0;
%         H=hfig;
%         hax=get(hfig, 'CurrentAxes');
%         XLim = get(hax, 'Xlim');
%         YLim = get(hax, 'Ylim');
%  %       Ylabel= get(get(hax,'ylabel'),'string');
% %        Zlabel= get(get(hax,'zlabel'),'string');
%         
%         if ~any(opt==[4 5])
%             i_keep(transpose(N0))
%         end
%         old_Y = g_Y;
%         oldN0 = ud.oldNO;
%         if ~isempty(ndx) && ((opt==1) || (opt==2))
%             
%             ud2 = get(hax, 'userdata');
%             if isfield(ud2, 'quasiss')&&ud2.quasiss
%                 %if quasi steady state interpolate eq value from userdata
%                 mask = ud2.mask;
%                 siz1 = size(ud2.qssvalues);
%                 if length(siz1) == 2
%                     siz1(3) = 1;
%                 end
% 
%                 qq = zeros(siz1(3), 1);
%                 initpt=N0(~mask);
%                 for i = 1:siz1(3)
%                     qq(i) = interp2(linspace(XLim(1), XLim(2), siz1(1)), linspace(YLim(1), YLim(2), siz1(2)), ud2.qssvalues(:, :, i), initpt(1), initpt(2));
%                 end
% 
%                 N0(mask) = qq;
%             else
%             vals =  cell(1, length(ndx));
%             for i = 1:length(ndx)
%                 vals{i} = num2str(N0(ndx(i)));
%             end
% 
%             answer = inputdlg(i_statevars_names(ndx), 'Not all initial conditions are on the axes, if needed edit these', 1, vals);
%             if ~isempty(answer)
%                 for i = 1:length(answer)
%                     N0(ndx(i)) = str2double(answer{i});
%                 end
% 
%             end
% 
%             end
% 
%             i_keep(transpose(N0));
%         end
% 
%         switch opt
%             case 1
%                 if isfield(g_grind, 'tilman')
%                     BtnName=questdlg('Do you want to set the resource supply parameters?','Yes');
%                     if strcmp(BtnName, 'Yes')
%                         oldhold = ishold;
%                         try
%                             evalin('base', g_grind.tilman.setS);
%                             hold on;
%                             where('sblink','go')
%                             if ~oldhold
%                                 hold off;
%                             end
% 
%                         catch err
%                             %                err=lasterror;
%                             if ~oldhold
%                                 hold off;
%                             end
%                             disp('error in evaluating g_grind.tilman.setS');
%                             rethrow(err);
%                         end
% 
%                     end
% 
%                     ru;
%                 else
%                     ru;
%                 end
% 
%             case 2
%                 findeq;
%             case 4
%                 docancel(ud.ix, ud.iy, oldN0, ud.p1, ud.p2, old_Y);
%                 plotedit on;
%                 if ishandle(H)
%                     set(H,'Pointer','arrow');
%                 end
% 
%             case 5
%                 docancel(ud.ix, ud.iy, oldN0, ud.p1, ud.p2, old_Y);
%         end
%         %       if isfield(ud, 'iters')
%         %          g_grind.solver.iters = g_grind.solver.iters / ud.iters;
%         %       end
% 
%         %  catch
%         %     if isfield(ud, 'iters')
%         %        g_grind.solver.iters = g_grind.solver.iters / ud.iters;
%         %    end
% 
%         % end
        
    case 'mmove'
        hax = get(hfig, 'CurrentAxes');
        if ishandle(hax)
            initpt = get(hax, 'CurrentPoint');
            Lims = zeros(3, 2);
            Lims(1,:) = get(hax, 'Xlim');
            Lims(2,:) = get(hax, 'Ylim');
            Lims(3,:) = get(hax, 'Zlim');
            axlabels={'xlabel','ylabel','zlabel'};
            fax=find(abs(initpt(1, :) - initpt(2, :)) <= 0.05);
            nval = 0;
            if length(fax) > 1
                s = sprintf(' | ');
                for i = 1:length(fax)
                    if ~((initpt(1, fax(i)) < Lims(fax(i), 1)) || (initpt(1, fax(i)) > Lims(fax(i), 2)))
                        lab = get(get(hax, axlabels{fax(i)}), 'string');
                        lab = clearcode(lab);
                        if isempty(lab) || (length(lab) > 8)
                            lab = [axlabels{fax(i)}(1) '-axis'];
                        end

                        nval = nval + 1;
                        s=sprintf('%s  %s = %0.3g', s, lab, initpt(1, fax(i)));
                    end

                end

            else
                s = sprintf(' | ');
                fax = 1:3;
                for i = 1:length(fax)
                    if ~((initpt(1, fax(i)) < Lims(fax(i), 1)) || (initpt(1, fax(i)) > Lims(fax(i), 2)))
                        lab = get(get(hax, axlabels{fax(i)}), 'string');
                        lab = clearcode(lab);
                        if isempty(lab) || (length(lab) > 8)
                            lab = [axlabels{fax(i)}(1) '-axis'];
                        end

                        nval = nval + 1;
                        s=sprintf('%s  %s = %0.3g', s, lab, initpt(1, fax(i)));
                    end

                end

            end

            if nval > 1
                set(hfig,'Pointer','crosshair')
            else
                s = '';
                set(hfig,'Pointer','arrow');
            end

            t1 = get(hfig, 'name');
            f1 = strfind(t1, ' | ');
            if ~isempty(f1)
                t1 = t1(1:f1(1) - 1);
            end

            set(hfig, 'name', [t1 s]);
        end

    case 'mdown2'
        %  global g_t g_Y g_grind;
        N0 = i_initvar;
        hax = get(hfig, 'CurrentAxes');
        initpt = get(hax, 'CurrentPoint');
        initpt = initpt(1, [1, 2]);
        XLim = get(hax, 'Xlim');
        YLim = get(hax, 'Ylim');
        if ~((initpt(1) < XLim(1)) || (initpt(1) > XLim(2)) || (initpt(2) < YLim(1)) || (initpt(2) > YLim(2)))
            oldt = t;
            t = initpt(1); 
            iy = i_getno(avar);
            if ~isempty(iy.no)
                N0(iy.no) = initpt(2);
            end

            i_keep(transpose(N0));
            ru
            t = oldt;
            if g_grind.statevars.dim > 1
                htitle=title(['Other intial conditions ' i_othervars(N0, iy.no)]);
                set(htitle,'fontweight','normal');
            end

            nextpen;
        end

    case 'keypressed'
        key = double(get(hfig, 'CurrentCharacter'));
        if isempty(key)
            return
        end
        switch key
            case 1 %Ctrl - A
                checkdlg(hfig);
                hax = get(gcf, 'CurrentAxes');
                initpt = get(hax, 'CurrentPoint');
                if ~isempty(initpt)
                    i_draw2Darrow(initpt(1, 1), initpt(1, 2));
                end

            case 65 %Shift - A
                checkdlg(hfig);
                hax = get(hfig, 'CurrentAxes');
                initpt = get(hax, 'CurrentPoint');
                if ~isempty(initpt)
                    arrows('delnearest', [initpt(1, 1), initpt(1, 2)]);
                end

            case 2 %Ctrl - B
                checkdlg(hfig);
                backw;
            case 11 %Ctrl - K
                checkdlg(hfig);
                ke;
            case 18 %Ctrl - R
                checkdlg(hfig);
                ru;
            case 20 %Ctrl - T
                checkdlg(hfig);
                time;
            case 21 %Ctrl - U
                conteq('-r');
                i_figure(i_figno('phase2'));
            case 27 %Esc
                if hfig == i_figno('dialog')
                    close(hfig)
                else
                    ud = get(hfig, 'userdata');
                    if ~isempty(ud) && isfield(ud, 'stop')
                        ud.stop = 1;
                        set(hfig, 'UserData', ud);
                    end

                end

            case 5 %Ctrl - E
                checkdlg(hfig);
                eigen;
            case 64 %Shift - 2
                checkdlg(hfig);
                null;
            case 35 %Shift - 3
                checkdlg(hfig);
                null3;
            case 25 %Ctrl - Y
                checkdlg(hfig);
                hax = get(hfig, 'CurrentAxes');
                initpt = get(hax, 'CurrentPoint');
                if ~isempty(initpt)
                    arrows('nearest', [initpt(1, 1), initpt(1, 2)]);
                end

            case 89 %Shift - Y
                checkdlg(hfig);
                hax = get(hfig, 'CurrentAxes');
                initpt = get(hax, 'CurrentPoint');
                if ~isempty(initpt)
                    arrows('delnearest', [initpt(1, 1), initpt(1, 2)]);
                end

            case 127 %delete
                cls;
            case 13 %Enter
                if hfig == i_figno('dialog')
                    checkdlg(hfig);
                    drawnow;
                    ru;
                end

                
            case 78 %Shift - N (annotation updating text)
                %       hax = get(hfig, 'CurrentAxes');
                %       initpt = get(hax, 'CurrentPoint');
                %       h1=annotation(hfig,'textbox','Position',[initpt(1, 1), initpt(1, 2),100,100],'tag','Annotation','LineStyle','none','String','A=xx');
                %
                %     case 6
                %       i_movie('rewind');
                %       refresh;
                %       i_movie('play', 1);
            %otherwise
              %  checkdlg(hfig);
                  % disp(key);
               % dokeypress(hfig)
        end

        
end

function checkdlg(hfig)
if hfig == i_figno('dialog')
    close(hfig);
end

function setinitialvalues(ud)
N0=i_initvar;
N0=setinitialvalue(ud.meta.xname,ud.initpt(1),N0);
N0=setinitialvalue(ud.meta.yname,ud.initpt(2),N0);
if length(ud.initpt)>2
  N0=setinitialvalue(ud.meta.zname,ud.initpt(3),N0);
end
if isfield(ud, 'quasiss')&&ud.quasiss
    %if quasi steady state interpolate eq value from userdata
    mask = ud.mask;
    siz1 = size(ud.qssvalues);
    if length(siz1) == 2
        siz1(3) = 1;
    end  
    qq = zeros(siz1(3), 1);
    initpt=N0(~mask);
    for i = 1:siz1(3)
        qq(i) = interp2(linspace(ud.meta.xlim(1), ud.meta.xlim(2), siz1(1)), linspace(ud.meta.ylim(1), ud.meta.ylim(2), siz1(2)), ud.qssvalues(:, :, i), initpt(1), initpt(2));
    end
    N0(mask) = qq;
end
i_keep(N0);

function N0=setinitialvalue(avar, aval, N0)
if ~any(avar=='''')&&~strcontains(avar,'(t+')
    iX=i_getno(avar,true);
    N0=i_setparvar(iX, N0, aval);
end

function resetpoint(H, oldpoint)
if ishandle(H)
    set(H, 'Pointer', oldpoint);
    set(H, 'WindowButtonDown',@(hobj,ev)i_callb('mdown',hobj));
    set(H, 'WindowButtonMotionFcn',@(hobj,ev)i_callb('mmove',hobj));
end

return

function docancel(ix, iy, oldN0, p1, p2, old_Y)
global g_Y;
i_setparvar(ix, p1, []);
i_setparvar(iy, p2, []);
i_keep(transpose(oldN0));
g_Y = old_Y;

function s2 = clearcode(s1)
s2=strrep(s1,'_','');
s2=strrep(s2,'\','');
s2=strrep(s2,'{','(');
s2=strrep(s2,'}',')');
