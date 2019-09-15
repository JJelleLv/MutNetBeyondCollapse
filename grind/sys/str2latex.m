%STR2LATEX - converts a MATLAB equation to a Latex equation
%   This function converts a string to Latex (also the simple MS Word variant)
%   All precedence rules are followed according to the MATLAB standard.
%   The resulting equation is showed in a dialog box
%
%   [latex,wordlatex] = STR2LATEX(STR) - makes a latex equation and prepares
%       a string for Word 2010 (paste in Word 2010, convert to equation and press
%      'Professional')
%   STR2LATEX(STR) - shows the latex equation in MATLAB
%   STR2LATEX({STR1,STR2,STR3}) - you can also enter a cell string
function [res, wordlatex] = str2latex(s)
debugcell =0;
%timesymb='\times ';
timessymb='\; ';
if strncmp(s, '$${', 3) %is already latex
   if nargout == 0
      viewlatex({s});
   end
   if nargout > 1
      res = s;
   end
   return;
end
if iscell(s)||isa(s,'parsed_equation')&&(length(s)>1)
   s1 = cell(length(s), 1);
   w1 = s1;
   tabs = '';
   for i = 1:length(s)
      if isa(s,'parsed_equation')
          s2=char(s(i));
      else
         s2=s{i};
      end
      f=strfind(s2,'%');
      if ~isempty(f)
          s2=s2(1:f-1);
      end
      haselse= ~isempty(regexp(s2,'\<else\>', 'once'))...
              ||~isempty(regexp(s2,'\<elseif\>', 'once'))||~isempty(regexp(s2,'\<case\>', 'once'));
      if (~isempty(regexp(s2,'\<end\>', 'once'))||haselse)&&~isempty(tabs)
      %    ~isempty(intersect(fs1,{'end','else','elseif','case'}))&&~isempty(tabs)
         tabs = tabs(1:end - 7);
      end
      s1{i} = '';
      %ok = true;
      if isa(s,'parsed_equation')
         [s1{i},w1{i}] = str2latex(s(i));
      else
         [s1{i},w1{i}] = str2latex(s{i});
      end
      if ~isempty(tabs)&&length(s1{i}) > 4
         s1{i} = [s1{i}(1:3) tabs s1{i}(4:end)];
      end
      if (~isempty(regexp(s2,'\<if\>', 'once'))||haselse)
    %  if ~isempty(intersect(fs1,{'if','else','elseif','case'}))
         tabs = [tabs '\;\;\; '];
      end
      if debugcell
         if ~checklatex(s1{i})
            fprintf('Error in line %d\n', i);
            error('Error in latex');
            s1{i} = '$${\textrm{Cannot display error in latex}}$$';
         end
      end
   end
   if nargout == 0
      viewlatex(s1, s);
   end
   if nargout >= 1
      res = s1;
   end
   if nargout >= 2
      wordlatex = w1;
   end

   return;
end
if isa(s,'parsed_equation')
    [fs,struc]=s.addbrackets;
else
    obj=parsed_equation(s);
    [fs,struc]=addbrackets(obj,true);
    %[fs,struc] = parseeq(s, true);
end
orgfs=fs;
%l_variable = 40;
isvars=find(struc>=parsed_equation.vr.funcname);
for i=1:length(isvars)
    fs{isvars(i)}=beautvar(fs{isvars(i)});
end
nooper = strcmp(fs,' ');
%if there is any two vars/numbers next to eachother it is assumed not an
%equation
if any(nooper)
   f=strcmp(fs,'«')|strcmp(fs,'»');
   comm = fs{end};
   if strncmp(comm, '%', 1)
      f(end) = true;
   else
      comm = '';
   end
   s = cleanuptext(sprintf('%s', fs{~f}));
   if ~isempty(comm)
      fs = {sprintf('\\textrm{%s}', s) comm};
   else
      fs = {sprintf('\\textrm{%s}', s) comm};
   end
end

f=find(strcmp('=', fs));
if length(f) == 1
   fstart = fs(1:f);
   fs = fs(f + 1:end);
   f3=strcmp(fstart,'«')|strcmp(fstart,'»');
   fstart = fstart(~f3);
   f=find(strcmp('''',fstart), 1);
   if ~isempty(f)
      fstart={sprintf('\\frac{d%s}{dt}=', fstart{1})};
   end
else
   fstart = {};
end
%fs = addbrackets(fs);
%make scientific notation prettier
for i = 1:length(fs)
   if isempty(fs{i})|| strcontains('0123456789', fs{i}(1))
      s1 = lower(fs{i});
      if ~strcontains(s1, '\times')
         f = strfind(s1, 'e');
         if ~isempty(f)
            fs{i}=[s1(1:f(1)-1) '\times 10^{' s1(f(1)+1:end) '}'];
         end
      end
   end
end
%handle comments
f = find(strncmp('%', fs, 1));
if ~isempty(f)
   s = cleanuptext(fs{f});
   %\text is not supported by MATLAB :-(
   fs{f} = sprintf('\\textrm{\\%s} ', s);
end
% if any(strcmp({'defextern ','defpermanent ','definepars '},fs{1}))
%    fs={sprintf('\\textrm{%s}',sprintf('%s',fs{:}))};
% end


f = find(strcmp('exp', fs));
for j = length(f):-1:1
   fs=[fs(1:f(j)-1) {'e','^'} fs(f(j)+1:end)];
end
f = find(strcmp('transpose', fs));
while ~isempty(f)
   brackets=(strcmp(fs,'(')|strcmp(fs,'«'))-(strcmp(fs,')')|strcmp(fs,'»'));
   j2 = popbrack(brackets, f(1) + 1, 1);
   if j2-f(1)<=3
   fs = [fs(1:f(1) - 1) {'«'} fs(f(1) + 2:j2 - 1) { '^\top','»'} fs(j2+1 :end) ];   
   else
   fs = [fs(1:f(1) - 1) {'«','('} fs(f(1) + 2:j2 - 1) { ')', '^\top','»'} fs(j2+1 :end) ];   
   end       
   f = find(strcmp('transpose', fs));
end

f = find(strcmp('.^', fs));
while ~isempty(f)
   brackets=(strcmp(fs,'(')|strcmp(fs,'«'))-(strcmp(fs,')')|strcmp(fs,'»'));
   fs{f(1)} = '\power\';
   j2 = popbrack(brackets, f(1) + 1, 1);
   h = 0;
   if strcmp(fs{j2},')')&&strcmp(fs{f(1)+1},'(')
      h = 1;
   end
   fs=[fs(1:f(1)) {'«','{'} fs(f(1)+h+1:j2-h) {'}','»'} fs(j2+1:end) ];
   f = find(strcmp('.^', fs));
end

fs=cellrep(fs,'\power\','.^');
f = find(strcmp('^', fs));
while ~isempty(f)
   brackets=(strcmp(fs,'(')|strcmp(fs,'«'))-(strcmp(fs,')')|strcmp(fs,'»'));
   fs{f(1)} = '\power\';
   j2 = popbrack(brackets, f(1) + 1, 1);
   h = 0;
   if strcmp(fs{j2},')')&&strcmp(fs{f(1)+1},'(')
      h = 1;
   end
   fs=[fs(1:f(1)) {'«','{'} fs(f(1)+h+1:j2-h) {'}','»'} fs(j2+1:end) ];
   f = find(strcmp('^', fs));
end

fs=cellrep(fs,'\power\','^');
f = find(strcmp('abs', fs));
while ~isempty(f)
   brackets=(strcmp(fs,'(')|strcmp(fs,'«'))-(strcmp(fs,')')|strcmp(fs,'»'));
   j2 = popbrack(brackets, f(1) + 1, 1);
   fs=[fs(1:f(1)-1) {'«','\left|'} fs(f(1)+2:j2-1) {'\right|','»'} fs(j2+1:end) ];
   f = find(strcmp('abs', fs));
end
f = find(strcmp('iif', fs));
if length(f) == 1 %if more than one
   brackets=(strcmp(fs,'(')|strcmp(fs,'«'))-(strcmp(fs,')')|strcmp(fs,'»'));
   j2 = popbrack(brackets, f(1) + 2, 1);
   cond = fs(f(1) + 2:j2);
   j3 = popbrack(brackets, j2 + 2, 1);
   eq1 = fs(j2 + 2:j3);
   j4 = popbrack(brackets, j3 + 2, 1);
   eq2 = fs(j3 + 2:j4);
   fs=[fs(1:f(1)-1) {'«' ,'\left\{ \begin{array}{ll}'} eq1 {'& \mbox{if $'} cond {'$} \\'} eq2 ...
      {'& \mbox{otherwise} \\\end{array}\right. \\ ', '»'} fs(j4+2:end)];
   %  fs=[fs(1:f(1)-1) {'«','\left|'} fs(f(1)+2:j2-1) {'\right|','»'} fs(j2+1:end) ];
   %  f = find(strcmp('iif', fs));
end
f = find(strcmp('sqrt', fs));
while ~isempty(f)
   brackets=(strcmp(fs,'(')|strcmp(fs,'«'))-(strcmp(fs,')')|strcmp(fs,'»'));
   j2 = popbrack(brackets, f(1) + 1, 1);
   fs=[fs(1:f(1)-1) {'«','\sqrt','{'} fs(f(1)+2:j2-1) {'}','»'} fs(j2+1:end) ];
   f = find(strcmp('sqrt', fs));
end

f = find(strcmp('/', fs));
while ~isempty(f)
   brackets=(strcmp(fs,'(')|strcmp(fs,'«'))-(strcmp(fs,')')|strcmp(fs,'»'));
   j1 = popbrack(brackets,f(1) - 1, - 1);
   j2 = popbrack(brackets, f(1) + 1, 1);
   if j1 == 0
      j1 = 1;
   end
   h = 0;
   if strcmp(fs{j1},'(')&&strcmp(fs{f(1)-1},')')
      h = 1;
   end
   h1 = 0;
   if strcmp(fs{j2},')')&&strcmp(fs{f(1)+1},'(')
      h1 = 1;
   end
   fs=[fs(1:j1-1) {'«','\frac{'} fs(j1+h:f(1)-1-h) {'}{'} fs(f(1)+h1+1:j2-h1) {'}','»'} fs(j2+1:end) ];
   f = find(strcmp('/', fs));
end
f=strcmp(fs,'«')|strcmp(fs,'»');
fs = [fstart fs(~f)];


fromf={'.*'    '(',     ')',       '''',    '<=',    '>=',     '~=',    '&',       '|',   '~',     '&&',   '||',     '*',  ...
   'asin', 'acos',    'atan' ,'return','function','if',   'elseif',          'else', 'end', '==','log'};
tof={' \cdot ',' \left(',' \right)','^\top ',' \leq ',' \geq ',' \neq ',' \land ',' \lor ',' \neg ',' \land ',' \lor ',timessymb,...
   '\arcsin','\arccos', '\arctan', '\textrm{return}', '\textrm{function}\;','\textrm{if}\;', '\textrm{elseif}\;','\textrm{else}','\textrm{end}','=','\ln'};
for i = 1:length(fs)
   f = find(strcmp(fs{i}, fromf), 1);
   if ~isempty(f)
      fs{i} = tof{f};
   end
end

res1=sprintf('$${%s}$$',sprintf('%s',fs{:}));
if nargout == 0
   viewlatex({res1}, {s});
end
if nargout > 0
   res = res1;
   if nargout > 1
      wordlatex = str2wordlatex(orgfs);
   end
end

function s = cleanuptext(s)
s = strrep(s,'#',' ');
s = strrep(s,'\',' ');
s = strrep(s,'\',' ');
s = strrep(s,'^',' ');
s = strrep(s,'ÿ','-');
s=strrep(s,'ö','o');
s = strrep(s,'_','\_');
s(2:end) = strrep(s(2:end),'%','-');
s = strrep(s,'&','and');
s =strrep(s,sprintf('\t'),'  ');

function s = beautvar(s)
f = strfind(s, '_');
if length(f) > 1||(length(f)==1&&f(1)==length(s))
   s=strrep(s,'_','\_');
end
if length(f) == 1
   f = strfind(s, '_');
   if f(1) < length(s) - 1
      s=[s(1:f(1)),'{' s(f(1)+1:end) '}'];
   end
end
greek={'Delta','Gamma','Lambda','Omega','Phi','Pi','Psi','Sigma',...
   'Theta','Upsilon','Xi','alpha','beta','chi',...
   'delta','epsilon','eta','gamma','iota','kappa','lambda','mu','ni','nu',...
   'omega',  'phi','pi','psi','rho','sigma','tau','theta','upsilon',...
   'varsigma','vartheta','xi','zeta'};
specials = [greek {'sin','cos','tan','ln'}];
if any(strcmp(s,specials))
    s=sprintf('\\%s',s);
end
if strcontains('1234567890', s(end))&&isletter(s(end - 1))
   h = '';
   s1 = s(1:end - 1);
   if any(strcmp(s1, greek))
      h = '\';
   end
   s = [h s(1:end - 1) '_' s(end)];
end
s=strrep(s,'#','{t+1}'); %trick to get this right
function res = str2wordlatex(fs)
if ischar(fs)
    obj=parsed_equation(fs);
    [fs]=addbrackets(obj,true);
%    fs = parseeq(fs, true);
elseif isa(fs,'parsed_equation')
    fs=fs.fs;
end
f=find(strcmp('=', fs));
if length(f) == 1
   fstart = fs(1:f);
   fs = fs(f + 1:end);
   f3=strcmp(fstart,'«')|strcmp(fstart,'»');
   fstart = fstart(~f3);
   f=find(strcmp('''',fstart), 1);
   if ~isempty(f)
      fstart={sprintf('d%s/dt =', fstart{1})};
   end
else
   fstart = {};
end

%fs = addbrackets(fs);
%make scientific notation prettier
for i = 1:length(fs)
   if strcontains('0123456789', fs{i}(1))
      s1 = lower(fs{i});
      if ~strcontains(s1, '\times')
         f = strfind(s1, 'e');
         if ~isempty(f)
            fs{i}=[s1(1:f(1)-1) '\times 10^(' s1(f(1)+1:end) ')'];
         end
      end
   end
end
f = find(strcmp('exp', fs));
for j = length(f):-1:1
   fs=[fs(1:f(j)-1) {'e','^'} fs(f(j)+1:end)];
end
f = find(strcmp('abs', fs));
while ~isempty(f)
   brackets=(strcmp(fs,'(')|strcmp(fs,'«'))-(strcmp(fs,')')|strcmp(fs,'»'));
   j2 = popbrack(brackets, f(1) + 1, 1);
   if strcmp(fs{j2},')')&&strcmp(fs{f(1)+1},'(')
      h = 1;
   end
   fs=[fs(1:f(1)-1) {'«','\abs\'} fs(f(1)+h+1:j2-h) {'\abs\','»'} fs(j2+h:end) ]; %temporary, else it is interpreted as or
   f = find(strcmp('abs', fs));
end
f=strcmp(fs,'«')|strcmp(fs,'»');
fs = [fstart fs(~f)];
specials={'Delta','Gamma','Lambda','Omega','Phi','Pi','Psi','Sigma',...
   'Theta','Upsilon','Xi','alpha','beta','chi',...
   'delta','epsilon','eta','gamma','iota','kappa','lambda','mu','ni','nu',...
   'omega',  'phi','pi','psi','rho','sigma','tau','theta','upsilon',...
   'varsigma','vartheta','xi','zeta','sqrt'};
fromf={'.*',     '''',    '<=',    '>=',     '~=',    '&',       '|',   '~',     '&&',   '||',     '*',  ...
   'asin', 'acos',    'atan'};
tof={' \cdot ','^\top ',' \leq ',' \geq ',' \neq ',' \land ',' \lor ',' \neg ',' \land ',' \lor ',' ',...
   'arcsin','arccos', 'arctan'};
for i = 1:length(fs)
   f = find(strcmp(fs{i}, fromf), 1);
   if ~isempty(f)
      fs{i} = tof{f};
   else
      f = find(strcmp(fs{i}, specials), 1);
      if ~isempty(f)
         fs{i} = sprintf('\\%s', specials{f});
      end
   end
end
fs=cellrep(fs,'\abs\','|');
res = sprintf('%s', fs{:});



function isok = checklatex(latexstring)
lastwarn('');
f = figure;
set(f, 'position', [10000, 100, 10, 10]) %somewhere out of the screen
text('Interpreter','latex','String',latexstring, 'position',[0.01,0.95]);
drawnow;
isok = isempty(lastwarn);
delete(f);

function viewlatex(latexstring, origstring)
if nargin == 1
   origstring = {};
end
ud.latex = latexstring;
ud.extents=nan(length(ud.latex),2);
ud.orig = origstring;
ch = get(0, 'children');
latexview=strcmp(get(ch,'tag'),'equationsdlg');
if any(latexview)
   h0 = ch(latexview);
   figure(h0);
   set(h0, 'userdata', ud);
   settoprow(h0, 1);
   return;
end

%figure(1)
hfig = figure('Units','pixels', ...
   'Color',[0.914 0.914 0.914], ...
   'MenuBar','none', ...
   'IntegerHandle', 'off',...
   'Name','Equations', ...
   'PaperPosition',[18 180 576 432], ...
   'PaperType','A4', ...
   'PaperUnits','points', ...
   'NumberTitle','off', ...
   'Position',[198  214  700 500], ...
   'Tag','equationsdlg', ...
   'ToolBar','none');

set(hfig, 'userdata', ud);



% panel1 = uipanel('Parent',hfig);
% panel2 = uipanel('Parent',panel1);
% set(panel1,'Position',[0 0.15 0.95 0.85]);
% set(panel2,'Position',[0 -3  3 1]);
%h = image;
axes('Parent', hfig,...
    'units','normalized',...
    'position',[0.02 0.10 0.90 0.85],...
    'units','normalized',...
    'ylim', [1 14],...
    'ydir','reverse',...
    'visible','off');
for i = 1:14
   text('Interpreter','latex','String','', 'position',[0.01,i],'fontsize',14,'tag',sprintf('latextext%d',i),'VerticalAlignment','top');
end
%pos=get(h,'Extent');
%Extent = [0.0075815 0.739479 0.284306 0.416834]
%set(h,'position',[0.01,0.95-pos(4)/2]);
h=uicontrol('Style','Slider','Parent',hfig,...
   'Units','normalized','Position',[1-0.022 0.09+0.03 0.022 1-(0.09+0.03)],...
   'Min', 1, ...
   'Max', length(ud.latex), ...
   'tag','vscroll',...
   'Value',length(ud.latex),'Callback',{@slider_Vcallback});
addlistener(h,'Value','PreSet',@(h1,h2)slider_Vcallback(h));
uicontrol('parent',hfig,...
    'style','text',...
    'units','normalized',...
    'position',[0 0 1 .09+0.03]);
h=uicontrol('parent',hfig,...
    'style','slider',...
    'units','normalized',...
    'position',[0 0.087 1-0.022 .03],...
    'min',0,...
    'max',1,...
    'callback',@slider_HCallback,...
    'tag','hscroll');
addlistener(h,'Value','PreSet',@(h1,h2)slider_HCallback(h));

uicontrol('Parent',hfig, ...
   'Units','normalized', ...
   'Callback',{@OKBtn}, ...
   'ListboxTop',0, ...
   'Position',[0.45 0.01 0.1 0.06], ...
   'String','OK', ...
   'Tag','OKBtn', ...
   'TooltipString','Close and save');
uicontrol('Parent',hfig, ...
   'Units','normalized', ...
   'Callback',{@CopyLatexBtn}, ...
   'ListboxTop',0, ...
   'Position',[0.65 0.01 0.1 0.06], ...
   'String','Copy Latex', ...
   'Tag','CopyBtn', ...
   'TooltipString','Copy latex to clipboard');
uicontrol('Parent',hfig, ...
   'Units','normalized', ...
   'Callback',{@CopyWordBtn}, ...
   'ListboxTop',0, ...
   'Position',[0.85 0.01 0.1 0.06], ...
   'String','Copy to Word', ...
   'Tag','CopyToWord', ...
   'TooltipString','Copy to Word 2010');

settoprow(hfig, 1);
% s = uicontrol('Style','Slider','Parent',hfig,...
%       'Units','normalized','Position',[0 0.1 0.95 0.05],...
%       'Value',0,'Callback',{@slider_Hcallback1});
function settoprow(hfig, no)
ud = get(hfig, 'userdata');
i = round(no);
h=findobj(hfig,'tag','vscroll');
if length(ud.latex) < 14
   set(h,'visible','off');
elseif (no == 1)
   set(h,'visible','on', 'SliderStep', [1, 5] / (length(ud.latex) - 1), ...
      'Max', length(ud.latex),'Value',length(ud.latex));
end
pos = 1;
maxext=0;
for j = 1:14
   h=findobj(hfig,'tag',sprintf('latextext%d',j));
   if (i > 0)&&(i<=length(ud.latex))
      set(h, 'String', ud.latex{i},'position',[0.01,pos]);
      if isnan(ud.extents(i,1))
          ex = get(h, 'extent'); %relatively slow
          ud.extents(i,1)=ex(3);
          ud.extents(i,2)=ex(4);
      end
      pos = pos + ud.extents(i,2) * 0.9;
      if ud.extents(i,1)>maxext
          maxext=ud.extents(i,1);
      end
   else
      set(h,'String','','position',[0.01,pos]);
      pos = pos + 1;
   end
   i = i + 1;
end
ud.maxext=maxext;
set(hfig,'userdata',ud);
slider_HCallback(h);


%%Callbacks
function slider_HCallback(src,~)
if nargin==0
    hfig=gcf;
else
    hfig=getparentfig(src);
end
ud=get(hfig,'Userdata');
hscr=findobj(hfig,'tag','hscroll');
aval=get(hscr,'value');
ch=get(get(hfig,'currentaxes'),'children');
%h=findobj(hfig,'tag','latextext');
pos=get(ch,'position');
if ud.maxext>1.05
    for i=1:length(ch)
       set(ch(i),'position',[(1-ud.maxext)*aval pos{i}(2:3)]);
    end
    set(hscr,'visible','on');
else
    for i=1:length(ch)
       set(ch(i),'position',[0 pos{i}(2:3)]);
    end
    set(hscr,'visible','off');
    set(hscr,'value',0);
end



function slider_Vcallback(src, ~)
aval = get(src, 'Value');
max = get(src, 'Max');
settoprow(get(src,'parent'), (max + 1) - aval);

function OKBtn(src, ~, ~)
close(get(src,'parent'));
function CopyLatexBtn(src, ~, ~)
ud = get(get(src,'parent'), 'userdata');
clipboard('copy',sprintf('%s\n',ud.latex{:}));
msgbox('Copied latex to clipboard','str2latex');

function CopyWordBtn(src, ~, ~)
ud = get(get(src,'parent'), 'userdata');
s = cell(length(ud.latex), 1);
for i = 1:length(ud.latex)
   s{i} = str2wordlatex(ud.orig{i});
end
clipboard('copy',sprintf('%s\n',s{:}));
msgbox('Copied string to clipboard','str2latex');

function fs =  cellrep(fs, sfrom, sto)
ks = find(strcmp(sfrom, fs));
for i = 1:length(ks)
   fs{ks(i)} = sto;
end
function j = popbrack(brack, i, ddir)
if ddir < 0
   brack = -brack;
end
if brack(i) == 0
   j = i;
else
   b = brack(i);
   j = i;
   while b > 0&&j > 0&&j<=length(brack)
      j = j + ddir;
      if j > 0
         b = b + brack(j);
      end
   end
end

