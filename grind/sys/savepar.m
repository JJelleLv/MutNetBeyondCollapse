%SAVEPAR   Save current settings
%    Save the model, the current parameter settings and the initial
%    conditions. If data is used in the model (<a href="matlab:help setdata">setdata</a> or <a href="matlab:help loaddata">loaddata</a>) 
%    it is attached to the inifile.
%    Note that the current values of the parameters and the state variables
%    are the new defaults!
%  
%    Usage:
%    SAVEPAR - saves the model and parameter settings to the current ini file
%    SAVEPAR FILE - saves the model and settings to FILE
%    SAVEPAR FILE 1 - saves the model as FILE and overwrites
%    existing files without confirmation
%    SAVEPAR('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'file' [file name] - Filename for saving the model
%     'keeprand' [logical] - Keep stochastic functions like rand in the inifile
%     'overwrite' [logical] - If true, the file is overwritten without notice
%
%  
%    See also savemodel 
%
%   Reference page in Help browser:
%      <a href="matlab:commands('savepar')">commands savepar</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function savepar(varargin)
%(afile, overwrite, keeprand)
global g_grind g_data;
fieldnams={'file', 'F', 'Filename for saving the model','';...
   'overwrite', 'l', 'If true, the file is overwritten without notice',false;...
   'keeprand', 'l', 'Keep stochastic functions like rand in the inifile',false}';
args=i_parseargs(fieldnams,'file,overwrite,keeprand', '',varargin);

if ~isfield(args,'keeprand')
   args.keeprand = false;
end
if ~isfield(args,'overwrite')
   args.overwrite = false;
end
if ~isfield(g_grind, 'inifile')
   i_errordlg('No model selected');
   error('GRIND:savepar:NoModel','No model to save');
end
if nargin == 0
   if isempty(g_grind.inifile)
      args.file=uiputfile('*.ini;*.all', 'Save model and parameters as');
   else
      if isempty(fileparts(g_grind.inifile))
         args.file=uiputfile([fullfile(cd,g_grind.inifile) ';*.all'], 'Save model and parameters as');
      else
         args.file=uiputfile([g_grind.inifile ';*.all'], 'Save model and parameters as');
      end
      if ~ischar(args.file)
          disp('Cancelled');
          return;
      end
   end
   args.overwrite = 1;
end
if isempty(strfind(args.file, '.all'))
   p = [par('full'),val()];
   for j = 1:length(p)
      s=strtrim(p{j}(1:strfind(char(p{j}), '=') - 1));
      p{j} = changecommand(s, p{j},args.keeprand);
      %      s=strtrim(p{j}(1:strfind(char(p{j}), '=') - 1));
      s2=strtrim(p{j}(strfind(char(p{j}), '=')  + 1:end));
      s2 = updateunits(s, s2);
      f = strfind(s2, ';');
      if ~isempty(f)&&(f(1)==length(s2))
         s2 = strtrim(s2(1:end - 1));
      end
      hasexp=changescheme(s,'%exp=', s2);
      if ~hasexp
         addscheme(s,'%exp=', s2);
      end
      if g_grind.statevars.vector
         ppar = evalin('base', s);
         hasr=changescheme(s, '%row=', sprintf('%d', size(ppar,1)));
         if ~hasr
            addscheme(s, '%row=', sprintf('%d', size(ppar,1)));
         end
         hasr=changescheme(s, '%col=', sprintf('%d', size(ppar,2)));
         if ~hasr
            addscheme(s, '%col=', sprintf('%d', size(ppar,2)));
         end
      end
   end
   [~, comm] = par('-groups');
   removecommands('par group ');
   if ~isempty(comm)
      g_grind.commands = [g_grind.commands , comm];
      for i = 1:length(g_grind.pars)
         changescheme(g_grind.pars{i}, '%tex=', g_grind.pargroups{i});
      end
   end
   N0 = i_initvar;
   if g_grind.statevars.vector
      for j = 1:length(g_grind.statevars.vectnames)
         kk = i_findparcom(g_grind.statevars.vectnames{j});
         if isempty(kk) || ~isempty(strfind(g_grind.commands{kk}, '['))
            ppar = evalin('base', g_grind.statevars.vectnames{j});
            s=sprintf('%s = %s;', g_grind.statevars.vectnames{j},mat2str(ppar,15));
%             s=sprintf('%s =[', g_grind.statevars.vectnames{j});
%             for k = 1:size(ppar, 1)
%                 s1=sprintf('%g, ',ppar(k,:));
%                 s=[s s1];
% %                for l = 1:size(ppar, 2)
% %                   s = sprintf('%s%g, ',s,ppar(k,l));
% %                end
%                if k < size(ppar, 1)
%                   s(length(s) - 1) = ';';
%                   s = sprintf('%s...\n    ', s);
%                end
%             end
%             s(length(s) - 1) = ']';
%             s(length(s)) = ';';
            changecommand(g_grind.statevars.vectnames{j}, s, 0);
            hasexp=changescheme(g_grind.statevars.vectnames{j}, '%exp=', mat2str(ppar,15));
            if ~hasexp
               addscheme(g_grind.statevars.vectnames{j},'%exp=', s);
            end
         end
      end
   elseif ~isempty(g_grind.statevars.names)
      for j = 1:g_grind.statevars.dim
         s = i_statevars_names(j);
         [des,un] =  par('-d', s);
         if ~isempty(un)||~isempty(des)
            s2=sprintf('%s = %0.10g; %%%s,%s', s, N0(j),des,un);
         else
            s2=sprintf('%s = %0.10g;', s, N0(j));
         end
         changecommand(s, s2, 0);
         if ~isempty(un)||~isempty(des)
            updateunits(s, s2);
         end
         hasexp=changescheme(i_statevars_names(j), '%exp=', sprintf('%0.10g', N0(j)));
         if ~hasexp
            addscheme(i_statevars_names(j), '%exp=', sprintf('%0.10g', N0(j)));
         end
      end
   end
   if ~isempty(g_data)
      if ~any(strncmpi(g_grind.commands, 'loaddata -ini', 13))
         g_grind.commands{end + 1} = 'loaddata -ini';
         %added data to inifile
      end
      if any(strncmpi(g_grind.commands, 'loaddata -ini', 13))
         dat = find(strcmp(g_grind.scheme, '%[data]'), 1);
         if isempty(dat)
            dat = length(g_grind.scheme) + 1;
            g_grind.scheme{dat} = '%[data]';
         end
         thedata = setdata('-current');
         if ~isempty(thedata.varlist)
            s = sprintf('%s\t', thedata.varlist{:});
            dat = dat + 1;
            g_grind.scheme{dat} = s(1:end - 1);
            for i = 1:size(thedata.matrix, 1)
               s = sprintf('%.10g\t', thedata.matrix(i, :));
               dat = dat + 1;
               g_grind.scheme{dat} = s(1:end - 1);
            end
         end
      end
   end
end
savemodel(args.file, args.overwrite);
function s2 =  updateunits(s, s2)
if strcontains(s2, '%')
   s2 = strtrim(s2(1:strfind(s2, '%') - 1));
   [des,uni] = par('-d', s);
   if ~isempty(uni)
      hasuni=changescheme(s,'%uni=', uni);
      if ~hasuni
         addscheme(s,'%uni=', uni);
      end
   end
   % des = par('-d', s);
   if ~isempty(des)
      hasdes=changescheme(s,'%des=', des);
      if ~hasdes
         addscheme(s,'%des=', des);
      end
   end
end
% 
% function icomm = findcommand(pname)
% global g_grind;
% icomm = [];
% f = regexp(g_grind.commands, sprintf('\\<%s\\>', pname));
% fndx = find(~cellfun(@isempty, f));
% for k = length(fndx):-1:1
%    f=strfind(char(g_grind.commands{fndx(k)}), '=');
%    if ~isempty(f)
%       s2 = strtrim(g_grind.commands{fndx(k)}(1:f(1) - 1));
%       if strcmp(pname, s2)
%          icomm = fndx(k);
%          return;
%       end
%    end
% end

function comm = changecommand(pname, pnew,keeprand)
global g_grind;
k = i_findparcom(pname);
comm = pnew;
if ~isempty(k)
   if keeprand
      keeprand =  ~(isempty(strfind(g_grind.commands{k}, 'rand('))&& isempty(strfind(g_grind.commands{k}, 'randn(')) ...
         && isempty(strfind(g_grind.commands{k}, 'randi('))&& isempty(strfind(g_grind.commands{k}, 'nestedmatrix(')));
   end
   if ~keeprand
      g_grind.commands{k} = pnew;
   else
      comm = g_grind.commands{k};
   end
elseif ~strcontains(pnew,'[]') %never save an uninitialized parameter
   %if not found then add
   g_grind.commands = [g_grind.commands, {pnew}];
end


function removecommands(pname)
global g_grind;
if length(g_grind.commands) > 1
   g_grind.commands = g_grind.commands(~strncmp(pname, g_grind.commands, length(pname)));
else
   f = strncmp(pname, g_grind.commands, length(pname));
   if f
      g_grind.commands = '';
   end
end

function addscheme(pname, plabel, pnew)
global g_grind;
if isfield(g_grind, 'scheme')&&~strcmp(pnew,'[]')
   k=find(strcmp(g_grind.scheme, sprintf('$sym=%s', pname)));
   % k=1;
   while k < length(g_grind.scheme)
      s = g_grind.scheme{k};
      if strncmp(s, '%sym=', 5) && strcmp(s(6:end), pname)
         if (k < length(g_grind.scheme)) &&~(strncmp(g_grind.scheme{k}, '%[', 2))
            g_grind.scheme = {g_grind.scheme{1:k}, sprintf('%s%s',plabel,pnew),g_grind.scheme{k + 1:end}};
         end
         return;
      end
      k = k + 1;
   end
end

function found = changescheme(pname, plabel, pnew)
global g_grind;
found = 0;
if isfield(g_grind, 'scheme')&&~strcmp(pnew,'[]')
   k=find(strcmp(g_grind.scheme, sprintf('%%sym=%s', pname)));
   if ~isempty(k)
      while (k < length(g_grind.scheme)) &&~(strncmp(g_grind.scheme{k}, '%[', 2) || strncmp(g_grind.scheme{k}, plabel, 5))
         k = k + 1;
      end
      if (k < length(g_grind.scheme))  && strncmp(g_grind.scheme{k}, plabel, 5)
         s = g_grind.scheme{k};
         found = 1;
         g_grind.scheme{k} = [s(1:5), pnew];
      end
   end
end
