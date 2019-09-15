%dialog box to view and change the Jacobians
function ret=i_enterjacdlg(flag,jac,makestruct)
global g_grind;
if nargin>0
    if nargin<3
        makestruct=false;
    end
    switch flag
        case 'str2jac'
            ret=str2jac(jac,[],makestruct);
        case 'jac2str'
            ret=jac2str(jac);
        otherwise
            error('grind:enterjacdlg','Unknown flag')
    end
    return;
end
h=findobj(0,'tag','i_enterjacdlg');
if isfunhandle(h)
    cancelpressed(h);
end
ud=getjacdata;
%ud.oldsyms=g_grind.syms;
ud.changed=false;
hassym=i_hastoolbox('symbolic');
name='Enter Jacobian';
if ~hassym
    vis='off';
    name='Enter Jacobian - no symbolic toolbox: some options disabled';
else
    vis='on';
end
if g_grind.statevars.vector
    ena='off';
    name='Enter Jacobian - Vector notation not supported for symbolic code generation';
else
    ena='on';
end
h=findobj(0,'Tag','i_enterjacdlg');
if ~isempty(h)
    delete(h);
end
hfig = figure(...
    'Units','characters',...
    'PaperUnits',get(0,'defaultfigurePaperUnits'),...
    'Color',[0.941176470588235 0.941176470588235 0.941176470588235],...
    'IntegerHandle','off',...
    'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
    'MenuBar','none',...
    'Name',name,...
    'NumberTitle','off',...
    'PaperPosition',get(0,'defaultfigurePaperPosition'),...
    'PaperSize',get(0,'defaultfigurePaperSize'),...
    'PaperType',get(0,'defaultfigurePaperType'),...
    'Position',[103.8 21.5384615384615 145.6 39.9230769230769],...
    'Resize','on',...
    'ResizeFcn',@resize,...
    'UserData',ud,...
    'Tag','i_enterjacdlg',...
    'CreateFcn',@(h,evnt)movegui(h, 'center'));

uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'BackgroundColor',[1 1 1],...
    'Callback',@Jacspopupclick,...
    'FontSize',10.6666666666667,...
    'Position',[6.4 36 68.4 1.53846153846154],...
    'TooltipString','Select matrix for analytical equations',...
    'String',ud.Jacnames,...
    'Style','popupmenu',...
    'Value',1,...
    'Tag','Jacspopup');
uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'Callback',@Generatepressed,...
    'TooltipString','Use MuPAD to generate symbolic code for Jacobians',...
    'Position',[77.6 35.4 18.8 2.46153846153846],...
    'String','Generate',...
    'Visible',vis,...
    'Enable',ena,...
    'Tag','Generatebutton' );
uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'TooltipString','Input function handles instead',...
    'Callback',@functpressed,...
    'Position',[117.6 35.4 18.8 2.46153846153846],...
    'String','Handles',...
    'Tag','Functbutton' );
uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'TooltipString','Clear all Jacobians',...
    'Callback',@Resetpressed,...
    'Position',[97.6 35.4 18.8 2.46153846153846],...
    'String','Clear all',...
    'Tag','Resetbutton' );
uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'HorizontalAlignment','left',...
    'Position',[6.4 32.1538461538462 38 1.07692307692308],...
    'String','Select Element of Jacobian:',...
    'Style','text',...
    'Tag','tabtext');

uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'BackgroundColor',[1 1 1],...
    'FontSize',10.6666666666667,...
    'Callback',@ElemEditChanged,...
    'HorizontalAlignment','left',...
    'Position',[6.4 30 40 1.7],...
    'String',{  '' },...
    'Style','edit',...
    'Tooltipstring','Enter the function handle to use for compting the Jacobian',...
    'visible','off',...
    'Tag','FunHanEdit');
%cell(length(ud.dx{1}),length(ud.dy{1}))
uitable(...
    'Parent',hfig,...
    'Units','characters',...
    'Data',{},...
    'CellSelectionCallback',@CellSelectCallback,...
    'Position',[6.2 20.5384615384615 88.8 11.6153846153846],...
    'ColumnName',ud.dx{1},...
    'RowName',ud.dy{1},...
    'Tooltipstring','Select element of Jacobian matrix',...
    'Tag','jactable' );

uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'Callback',@runmupadpressed,...
    'Position',[119.8 23.4615384615385 18.8 2.46153846153846],...
    'String','MuPAD',...
    'Tooltipstring','Enter a MUPAD command to adapt the current equation',...
    'Visible',vis,...
    'Tag','MuPADbutton' );
uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'TooltipString','Use the advanced Simplify option of MuPAD to simplify current (may take some time)',...
    'Callback',@(hObject,event)runmupadpressed(hObject,event,'Simplify(%,Steps=1000)'),...
    'Position',[119.8 19.1538461538462 18.8 2.46153846153846],...
    'Visible',vis,...
    'Enable',ena,...
    'String','Simplify',...
    'Tooltipstring','Use mupad to simplify the current equation',...
    'Tag','Simplitybutton');

uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'BackgroundColor',[1 1 1],...
    'FontSize',10.6666666666667,...
    'Callback',@ElemEditChanged,...
    'HorizontalAlignment','left',...
    'Position',[6.4 15.0769230769231 130.2 1.69230769230769],...
    'String',{  'Edit Text' },...
    'Tooltipstring','Edit the equation of the selected element of the Jacobian',...
    'Style','edit',...
    'Tag','ElemEdit');

uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'HorizontalAlignment','left',...
    'Position',[6.6 16.7692307692308 38 1.15384615384615],...
    'String','Edit element',...
    'Style','text',...
    'Tag','elemtext');

uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'Callback',@OKpressed,...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'Position',[24.4 1.38 18.8 2.46],...
    'String','OK',...
    'Tooltipstring','Close and keep changes',...
    'Tag','OKbutton' );

uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'Callback',@cancelpressed,...
    'FontSize',10.6666666666667,...
    'Position',[59.8 1.38 18.8 2.46],...
    'String','Cancel',...
    'Tooltipstring','Close abort changes',...
    'Tag','Cancelbutton' );

uicontrol(...
    'Parent',hfig,...
    'Units','characters',...
    'FontUnits','pixels',...
    'FontSize',10.6666666666667,...
    'Callback', 'commands(''enterjac'')',...
    'Position',[95.6 1.38 18.8 2.46],...
    'String','Help',...
    'Tooltipstring','Help on enterjac',...
    'Tag','Helpbutton');
% hpanel=uipanel('units','character','position', [6 7.23076923076923 126*3 5.69230769230769],'tag','hpanel');
% set(hpanel,'units','normalized');
h=uicontrol('parent',hfig,...
    'style','slider',...
    'units','normalized',...
    'position',[0.05 0.17 0.9 .03],...
    'min',0,...
    'max',1,...
    'callback',@hscroll_Callback,...
    'tag','hscroll');
% hax=axes('parent',hpanel,'outerposition',[0 0 1 1]);
if isoctave||~verLessThan('MATLAB','7.12')
    addlistener(h,'Value','PreSet',@(h1,h2)hscroll_Callback(h,h1));
end

hax=axes;
set(hax,'Units','character','position', [6 9 120 5.69230769230769]);
set(hax,'visible','off', 'ylim', [0 3],'ydir','reverse')
text('Interpreter','latex','String','', 'position',[0.0,1],'fontsize',11,'tag','latextext','VerticalAlignment','top');
Jacspopupclick(hfig);
%CellSelectCallback(hfig,[])

%%Callbacks

function hscroll_Callback(hObject,~)
if nargin==0
    hfig=gcf;
else
    hfig=getparentfig(hObject);
end
hscr=findobj(hfig,'tag','hscroll');
aval=get(hscr,'value');
h=findobj(hfig,'tag','latextext');
ext=get(h,'extent');
pos=get(h,'position');
if ext(3)>1
    set(h,'position',[(1-ext(3))*aval pos(2:3)]);
    set(hscr,'visible','on');
else
    set(h,'position',[0 pos(2:3)]);
    set(hscr,'visible','off');
    set(hscr,'value',0);
end


function ElemEditChanged(hObject,~)
hfig=getparentfig(hObject);
ud=get(hfig,'userdata');
h=findobj(hfig,'tag','ElemEdit');
eq=get(h,'string');
h=findobj(hfig,'Tag','Jacspopup');
aval=get(h,'Value');
if ~isfunhandle(ud.Jacs{aval})
    if isempty(ud.Jacs{aval})
        ud.Jacs{aval}=cell(length(ud.dx{aval}),length(ud.dy{aval}));
    end
    oldeq=ud.Jacs{aval}{ud.selectedcell(1),ud.selectedcell(2)};
    ud.Jacs{aval}{ud.selectedcell(1),ud.selectedcell(2)}=eq;
    ud.changed=~strcmp(oldeq,eq);
    set(hfig,'Userdata',ud);
    h=findobj(hfig,'tag','latextext');
    set(h,'string',str2latex(sprintf('%s_%d%d=%s',ud.Jacnames{aval}(1),ud.selectedcell(1),ud.selectedcell(2),eq)));
end


function Generatepressed(hObject,~)
hfig=getparentfig(hObject);
ud=get(hfig,'userdata');
set(hfig,'Pointer','watch')
drawnow();
ud.syms.errormsg={};
ud.syms=i_model2mupad(true,ud.syms);
if ~isempty(ud.syms.errormsg)
    errordlg(sprintf('Error in symbolic toolbox: %s',ud.syms.errormsg{1}))
else
    ud.changed=true;
    ud=update_Jacs(ud);
end
set(hfig,'userdata',ud);
Jacspopupclick(hObject)
set(hfig,'Pointer','arrow')

function resize(hobj,~)
hfig=getparentfig(hobj);
data={...
    'Jacspopup',[1 1 1 0],[24 351 256.5 15];...
    'Generatebutton',[ 0 1 1 0],[291 345.15 70.5 24];...
    'Functbutton',[0 1  1 0],[441 345.15 70.5 24];...
    'Resetbutton',[0 1 1 0],[366 345.15 70.5 24];...
    'tabtext',[1 0 1 0],[24 313.5 142.5 10.5];...
    'FunHanEdit',[1 0 1 0],[24 292.5 150 16.575];...
    'jactable',[1 1 1 1],[23.25 200.25 333 113.25];...
    'MuPADbutton',[0 1 0 1],[449.25 228.75 70.5 24];...
    'Simplitybutton',[0 1 0 1],[449.25 186.75 70.5 24];...
    'ElemEdit',[1 1  0 1],[24 147 488.25 16.5];...
    'elemtext',[1 0  0 1],[24.75 163.5 142.5 11.25];...
    'OKbutton',[1 0 0 1],[91.5 13.455 70.5 23.985];...
    'Cancelbutton',[1 0 0 1],[224.25 13.455 70.5 23.985];...
    'Helpbutton',[1 0 0 1],[358.5 13.455 70.5 23.985];...
    'hscroll',[1 1 0 1],[27.3 66.045 491.4 11.655];...
    };
i_resizedlg(hfig,data(:,1),data(:,2),data(:,3),[206.25 84.75 546 388.5]);
i_resizedlg(hfig,'horizcenter',{'OKbutton','Cancelbutton','Helpbutton'});


function Resetpressed(hObject,~)
%global g_grind;
hfig=getparentfig(hObject);
ud=get(hfig,'userdata');
btn=questdlg('Are you sure to clear all Jacobians?', ...
    'reset', 'Yes', 'No', 'No');
if strcmp('Yes',btn)
    f=fieldnames(ud.syms);
    for i=1:length(f)
        ud.syms.(f{i})={};
    end
    ud=update_Jacs(ud);
    ud.changed=true;
    set(hfig,'userdata',ud);
    Jacspopupclick(hObject)
end

function cancelpressed(hObject,~)
%global g_grind;
hfig=getparentfig(hObject);
% ud=get(hfig,'userdata');
% if isfield(ud,'oldsyms')
%     g_grind.syms=ud.oldsyms;
% end
delete(hfig);


function OKpressed(hObject,~)
global g_grind;
hfig=getparentfig(hObject);
ud=get(hfig,'userdata');
if ud.changed
    ud=update_ud(ud);
    
    g_grind.syms=ud.syms;
    %     g_grind.syms.Jacobian=ud.Jacs{1};
    %     g_grind.syms.Jacobianp=ud.Jacs{2};
    %     g_grind.syms.Sensitivp=ud.Jacs{3};
    %     if g_grind.solver.isimplicit
    %         g_grind.syms.Jacobian_y=ud.Jacs{4};
    %         g_grind.syms.Jacobian_yp=ud.Jacs{5};
    %     else
    %         g_grind.syms.Hessian=cell(g_grind.statevars.dim,g_grind.statevars.dim,g_grind.statevars.dim);
    %         if isfunhandle(ud.Jacs{3+1})
    %             g_grind.syms.Hessian=ud.Jacs{3+1};
    %         else
    %             allempty=true;
    %             for i=1:g_grind.statevars.dim
    %                 if ~isempty(ud.Jacs{3+i})
    %                     allempty=false;
    %                     g_grind.syms.Hessian(:,:,i)=ud.Jacs{3+i};
    %                 end
    %             end
    %
    %             if allempty
    %                 g_grind.syms.Hessian=[];
    %             end
    %         end
    %     end
end
delete(hfig);
enterjac('-?');

function ud=update_ud(ud)
hessndx=strncmp(ud.fieldnames,'Hessian,',8);
Hess1=cat(3,ud.Jacs{hessndx});
ud.fieldnames=ud.fieldnames(~hessndx);
ud.Jacs=ud.Jacs(~hessndx);
ud.fieldnames{end+1}= 'Hessian';
ud.Jacs{end+1}=Hess1;
for i=1:length(ud.Jacs)
    if i==3
        lst=sens_var_list;
    else
        lst=[];
    end
    if isstruct(ud.syms.(ud.fieldnames{i}))
        ud.syms.(ud.fieldnames{i})=str2jac(ud.Jacs{i},lst,true);
    else
        ud.syms.(ud.fieldnames{i})=str2jac(ud.Jacs{i},lst,false);
    end
end


%
% for i=1:length(ud.Jacs)
%     if ~isfunhandle(ud.Jacs{i})
%         if i==3
%             ud.Jacs{i}=str2jac(ud.Jacs{i},sens_var_list);
%         else
%             ud.Jacs{i}=str2jac(ud.Jacs{i});
%         end
%     end
% end
% for i=1:length(ud.Jacs)
%     f=strfind(ud.fieldnames{i},',');
%     if ~isempty(f)
%         if isempty(ud.Jacs{i})
%             ud.syms.(ud.fieldnames{i}(1:f-1))={};
%         else
%             j=str2double(ud.fieldnames{i}(f+1:end));
%             ud.syms.(ud.fieldnames{i}(1:f-1))(:,:,j)=ud.Jacs{i};
%         end
%     else
%         ud.syms.(ud.fieldnames{i})=ud.Jacs{i};
%     end
% end
% if iscell(ud.syms.Hessian)
% for i=1:size(ud.syms.Hessian,1)
%    for j=1:size(ud.syms.Hessian,2)
%        for k=1:size(ud.syms.Hessian,3)
%            if isempty(ud.syms.Hessian{i,j,k})
%                ud.syms.Hessian{i,j,k}='';
%            end
%        end
%    end
% end

function ud=update_Jacs(ud)
statevars=i_statevars_names;
for i=1:length(ud.fieldnames)
    f=strfind(ud.fieldnames{i},',');
    if ~isempty(f)
        j=str2double(ud.fieldnames{i}(f+1:end));
        fld=ud.fieldnames{i}(1:f-1);
        Hess=jac2str(ud.syms.(fld),statevars);
        if ~iscell(Hess)||isempty(Hess)
            ud.Jacs{i}=Hess;
        else
            ud.Jacs{i}=Hess(:,:,j);
        end
    else
        if i==3&&~isempty(ud.syms.(ud.fieldnames{i}))
            flds=sens_var_list;
        else
            flds=statevars;
        end
        ud.Jacs{i}=jac2str(ud.syms.(ud.fieldnames{i}),flds);
    end
end

function CellSelectCallback(hObject,eventdata)
hfig=getparentfig(hObject);
ud=get(hfig,'userdata');
h=findobj(hfig,'tag','Jacspopup');
aval=get(h,'value');
if ~(isa(eventdata,'matlab.ui.eventdata.CellSelectionChangeData')||isfield(eventdata,'Indices'))||isempty(eventdata.Indices)
    row=1;
    col=1;
else
    row=eventdata.Indices(1);
    col=eventdata.Indices(2);
end
if isfunhandle(ud.Jacs{aval})||isempty(ud.Jacs{aval})||row>size(ud.Jacs{aval},1)||col>size(ud.Jacs{aval},2)
    s='';
else
    s=ud.Jacs{aval}{row,col};
end
ud.selectedcell=[row col];
set(hfig,'userdata',ud);
h=findobj(hfig,'tag','ElemEdit');
set(h,'string',s);
h=findobj(hfig,'tag','latextext');
set(h,'string',str2latex(sprintf('%s_%d%d=%s',ud.Jacnames{aval}(1),row,col,s)));
%h=findobj(hfig,'tag','hscroll');
%set(h,'value',0);
hscroll_Callback(h);

% pos=get(h,'position');
% set(h,'position',[1-ext(3) pos(2:3)]);

function runmupadpressed(hObject,~,comm)
if nargin<3
    comm=inputdlg('Enter MuPAD command (%=current equation)','enterjac',1,{'Simplify(%,Steps=1000)'});
    if isempty(comm)
        return;
    end
    comm=comm{1};
end
hfig=getparentfig(hObject);
%ud=get(hfig,'userdata');
h=findobj(hfig,'tag','ElemEdit');
eq=get(h,'string');
if ~isempty(eq)
    set(hfig,'Pointer','watch')
    drawnow();
    peq=parsed_equation(eq);
    peq2=char(runmupad(peq,comm));
    set(h,'string',peq2);
    ElemEditChanged(h);
    set(hfig,'Pointer','arrow')
end

function functpressed(hObject,~)
hfig=getparentfig(hObject);
ud=get(hfig,'userdata');
prompt={'Function handle for the Jacobian','Function handle for the Jacobian of parameters','Function handle for the Sensitivities','Function handle for full Hessian'};
name='Enter Jacobian (leave empty if no handle should be used)';
numlines=1;
defaultanswer={'','','',''};
for i=1:length(defaultanswer)
    if isfunhandle(ud.Jacs{i})
        defaultanswer{i}=func2str(ud.Jacs{i});
    end
end
answer=inputdlg(prompt,name,numlines,defaultanswer);
if ~isempty(answer)
    for i=1:length(answer)
        if ~isempty(answer{i})
            ud.changed=true;
            ud.Jacs{i}=str2func(answer{i});
        elseif isfunhandle(ud.Jacs{i})
            ud.Jacs{i}='';
        end
    end
    set(hfig,'userdata',ud)
    Jacspopupclick(hfig)
end

function Jacspopupclick(hObject,~)
hfig=getparentfig(hObject);
h=findobj(hfig,'tag','Jacspopup');
aval=get(h,'value');
ud=get(hfig,'userdata');
if isfunhandle(ud.Jacs{aval})
    set(findobj(hfig,'Tag','jactable'),'Visible','off')
    set(findobj(hfig,'Tag','ElemEdit'),'Visible','off')
    set(findobj(hfig,'Tag','elemtext'),'Visible','off')
    set(findobj(hfig,'tag','hscroll'),'Visible','off')
    h1=findobj(hfig,'Tag','tabtext');
    set(h1,'string','Enter function handle')
    h1=findobj(hfig,'Tag','FunHanEdit');
    set(h1,'Visible','on')
    set(h1,'string',func2str(ud.Jacs{aval}));
    h1=findobj(hfig,'tag','latextext');
    set(h1,'string','');
else
    h=findobj(hfig,'tag','jactable');
    set(h,'Visible','on')
    set(findobj(hfig,'Tag','ElemEdit'),'Visible','on')
    set(findobj(hfig,'Tag','elemtext'),'Visible','on')
    set(findobj(hfig,'Tag','FunHanEdit'),'Visible','off')
    dx=ud.dx{aval};
    dy=ud.dy{aval};
    dat=cell(length(dy),length(dx));
    h1=findobj(hfig,'Tag','tabtext');
    if ~strcontains(ud.Jacnames{aval},'Hessian')
        set(h1,'String','Select Element of Jacobian:')
        for i=1:length(dx)
            s=sprintf(['d(%s)/d(' dx{i} ');'],dy{:});
            c= textscan(s,'%s','whitespace',';');
            dat(:,i)=c{1};
           % dat(:,i)=regexp(s(1:end-1),';','split');
        end
    else
        var3=regexp(ud.Jacnames{aval},'(?<=")[^"]*(?=")','match');
        set(h1,'String','Select Element of Hessian:')
        for i=1:length(dx)
            s=sprintf(['d2(%s)/d(' dx{i} ')/d(' var3{1} ');'],dy{:});
            dat(:,i)=regexp(s(1:end-1),';','split');
        end
    end
    set(h, 'Data',dat,'ColumnName',dx,'RowName',dy);
    CellSelectCallback(hObject,[]);
end
%%other functions

function var=sens_var_list
global g_grind;
statevars=i_statevars_names;
pars=[g_grind.pars statevars];
var=cell(length(statevars),length(pars)+1);
var(:,1)=transpose(statevars);
for i=1:length(pars)
    %efficient vectorized filling of matrix
    s=sprintf(['Sens_' pars{i} '_%s;'],statevars{:});
    var(:,i+1)=regexp(s(1:end-1),';','split');
end
var=transpose(var(:));

% function g_grind_changed(hfig)
% ud=get(hfig,'userdata');
% ud.changed=true;
% ud=getjacdata(ud);
% set(hfig,'userdata',ud);
% Jacspopupclick(hfig);
% CellSelectCallback(hfig,[])

function ud=getjacdata(ud)
global g_grind;
ud.syms=g_grind.syms;
ud.fieldnames={'Jacobian','Jacobianp','Sensitivp'};
%ud.Jacs={g_grind.syms.Jacobian,g_grind.syms.Jacobianp,g_grind.syms.Sensitivp};
funs=sprintf('f_%d,',1:g_grind.statevars.dim);
funs=regexp(funs(1:end-1),',','split');
statevars=i_statevars_names;
ud.Jacnames={'Jacobian of state variables (eigen/conteq)','Jacobian of parameters (conteq)','Sensitivity of parameters (timesens)'};
ud.dy=[{funs},{funs},{funs}];
ud.dx=[{statevars},{g_grind.pars},{[g_grind.pars statevars]}];
ud.selectedcell=[1 1];
if g_grind.solver.isimplicit
    if ~isfield(ud.syms,'Jacobian_yp')
        ud.syms.Jacobian_yp={};
    end
    if ~isfield(ud.syms,'Jacobian_y')
        ud.syms.Jacobian_y={};
    end
    ud.Jacnames=[ud.Jacnames,{'Jacobian of y','Jacobian of y'''}];
    %    ud.Jacs=[ud.Jacs,{g_grind.syms.Jacobian_y,g_grind.syms.Jacobian_yp}];
    ud.fieldnames=[ud.fieldnames,{'Jacobian_y','Jacobian_yp'}];
    funs=  regexp(sprintf('eq_%d,',1:g_grind.statevars.dim),',','split');
    funs=funs(1:end-1);
    dstates=regexp(sprintf('%s'',',statevars{:}),',','split');
    dstates=dstates(1:end-1);
    ud.dx=[ud.dx {statevars},{dstates}];
    ud.dy=[ud.dy {funs},{funs} ];
else
    if isfield(ud.syms,'Hessian')&&isfunhandle(ud.syms.Hessian)
        ud.fieldnames=[ud.fieldnames,{'Hessian'}];
        ud.Jacnames=[ud.Jacnames,{'Hessian matrix (full)'}];
    else
        ud.dy=[ud.dy repmat({funs},1,g_grind.statevars.dim)];
        ud.dx=[ud.dx repmat({statevars},1,g_grind.statevars.dim)];
        Hess=cell(1,numel(statevars));
        for i=1:length(Hess)
            Hess{i}=sprintf('Hessian(:,:,%d) derivatives of J with respect to "%s" (conteq)',i,statevars{i});
        end
        ud.Jacnames=[ud.Jacnames,Hess];
        for i=1:length(Hess)
            %             if isfield(g_grind.syms,'Hessian')&&~isempty(g_grind.syms.Hessian)
            %                 ud.Jacs=[ud.Jacs {g_grind.syms.Hessian(:,:,i)}];
            %             else
            %                 ud.Jacs=[ud.Jacs {[]}];
            %             end
            ud.fieldnames=[ud.fieldnames,{sprintf('Hessian,%d',i)}];
        end
    end
end
ud=update_Jacs(ud);
% for i=1:length(ud.Jacs)
%     if i==3
%         ud.Jacs{i}=jac2str(ud.Jacs{i},sens_var_list);
%     else
%         ud.Jacs{i}=jac2str(ud.Jacs{i});
%     end
% end


function strjacs=jac2str(jac,statevars)
global g_grind
if nargin<2||isempty(statevars)
    statevars=i_statevars_names;
end
if isfunhandle(jac)
    strjacs=jac;
    return
end
strct=isstruct(jac);
if strct
    jac1=jac;
    jac={jac1.unique(:).equation};
end
strjacs=cell(size(jac));
for j=1:numel(jac)
    strjac=jac{j};
    for i = 1:length(statevars)
        %replace statevar symbols with g_X1(i)
        strjac = regexprep(strjac, sprintf('(?<![a-zA-Z_0-9])g_X1[(]%d,[:][)]', i), statevars{i});
    end
    strjac=regexprep(strjac,'[\.]([\^*/]{1})','$1'); %remove vectorization for readability;
    %    strjac=regexprep(strjac,'([0-9]{1})[.][0]\>','$1'); %remove 1.0 ->1 or 2.0 ->2
    strjacs{j}=strjac;
end
if g_grind.solver.isimplicit
    for j=1:numel(jac)
        strjac=strjacs{j};
        for i = 1:length(statevars)
            %replace statevar symbols with g_X2(i)
            strjac = regexprep(strjac, sprintf('(?<![a-zA-Z_0-9])g_X2[(]%d,[:][)]', i), [statevars{i} '''']);
        end
        %   strjac=regexprep(strjac,'[\.]([\^*/]{1})','$1'); %remove vectorization for readability;
        %    strjac=regexprep(strjac,'([0-9]{1})[.][0]\>','$1'); %remove 1.0 ->1 or 2.0 ->2
        strjacs{j}=strjac;
    end
end
if strct
    jac=strjacs;
    strjacs=cell(jac1.size);
    strjacs(:)={'0'};
    if isfield(jac1,'unique')
        for i=1:length(jac1.unique)
            strjacs(jac1.unique(i).indices)=jac(i);
        end
    end
end



function jac=str2jac(strjac,statevars,dostruct)
jac=strjac;
if isfunhandle(jac)||isempty(jac)
    return
end
if nargin<2||isempty(statevars)
    statevars=i_statevars_names;
end
if nargin<3
    dostruct=false;
end
for i = 1:length(statevars)
    %replace statevar symbols with g_X1(i)
    jac=regexprep(jac,sprintf('(?<![a-zA-Z_0-9])%s(?![a-zA-Z_0-9])',statevars{i}),sprintf('g_X1(%d,:)', i));
end
allempty=true;
for i=1:size(jac,1)
    for j=1:size(jac,2)
        for k=1:size(jac,3)
            if ~isempty(jac{i,j,k})
                allempty=false;
                jac{i,j,k}=vectorize(jac{i,j,k});
            end
        end
    end
end
if allempty
    jac={};
end
if dostruct
    jac=makestruct(jac);
end

function res=makestruct(c)
res.size=size(c);
uni=unique(c(:));
k=1;
for i=1:length(uni)
    if ~any(strcmp({'','0'},uni{i}))
        res.unique(k)=struct('equation',uni{i},'indices',find(strcmp(uni{i},c(:))));
        k=k+1;
    end
end

function res=isfunhandle(f)
res=isa(f,'function_handle');