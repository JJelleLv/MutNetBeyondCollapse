%implements USE/MODEL. if used by program optionally no clear of global variables
function i_use(amodelname, doclear, doclose, docheck, doerase)
global g_grind;
if nargin < 3
   doclose = 1;
end

if nargin  < 4
   docheck = 1;
end

if nargin  < 5
   doerase = 1;
end

doerase = doerase && ~isempty(g_grind);
%if strcmp(pwd, fullfile(matlabroot, 'work'))
%   grindp = grindpath(0);
%   cd(grindp{1});
%end

%addpath(grindpath);
if ~isempty(amodelname)
   if ~strcontains(amodelname, '.')
      amodelname = [amodelname '.ini'];
   end

   if strcontains(amodelname,'*')
     [pathname, amodelname, iscurdir] = findgrindfiles(amodelname);
     if isempty(amodelname)
         error('GRIND:model:Cancelled','No files found');
     end

   else
     [pathname,  iscurdir] = findgrindfile(amodelname);
   end

   if isempty(pathname)
      i_errordlg(['Inifile: ' amodelname ' does''t exist']);
      error('GRIND:model:NoFile','File "%s" doesn''t exist',amodelname);
   elseif ~iscurdir
      disp(['Could not find the model "' amodelname '" in "' cd '"']);
      disp(['Changing to:"' pathname '"']);
      cd(pathname);
      [~,amodelname,ext]=fileparts(amodelname);
      amodelname=[amodelname ext];
      %i_use(amodel, doclear, doclose);
     % return;
   end

   replayall('-removeinfo')
   if doclose
      disp(['Using ' amodelname]);
      if doclear
         clear global;
      end

      [modl, comman, schem] = i_loadinifile(amodelname);
      i_makemodel(modl, comman, amodelname, schem);
      i_parcheck(1);
   else
      i_moddlg(doclear, docheck, ~doclose, doerase);
      if ~isempty(amodelname)
         h = findobj(gcf,'Tag','FileName');
         if ~strcontains(amodelname, '.')
            amodelname = [amodelname '.ini'];
         end

         set(h, 'String', amodelname);
         feval(get(h,'callback'),h,[])
      end

   end

else
   i_moddlg(doclear, docheck, ~doclose, doerase);
end


%FINDGRINDFILE - searches grind subdirectories for file/directory
%
%   Usage:
%   [pathname, iscurdir] = findgrindfile(afile)
function [pathname, filename, iscurdir] = findgrindfiles(afilename)
paths1 = {};
[pathname,filename] = fileparts(afilename);
if ~isempty(pathname)
   iscurdir=strcmp(pathname,pwd);
   return;
end

%Current directory
adir = pwd;
ss = dir([adir filesep afilename]);
ss = {ss.name};
for i=1:length(ss)
   if isempty(strfind(ss{i},filesep))
      ss{i}=[adir filesep ss{i}];
   end

end

paths1 = [paths1 ss];

grindp = grindpath(0);
for i = 1:length(grindp)
   paths1 = findinsubdirs(grindp{i}, afilename, paths1);
end


paths1=sort(paths1);
ndx=ones(size(paths1));
for i=length(paths1)-1:-1:1
   ndx(i)= ~strcmp(paths1{i},paths1{i+1});
end

paths1=paths1(logical(ndx)); 
if length(paths1) > 1
   curp=find(strcmp(paths1,[pwd filesep afilename]));
   if isempty(curp)
      curp=1;
    end

   [pathname,v] = listdlg('Name','More than one file found with that name','Promptstring','Select a file:','SelectionMode','single',...
      'ListString',paths1,'InitialValue',curp,'ListSize',[400,300]);
   if ~v
      pathname = '';
   else
      pathname = paths1{pathname};
   end

elseif ~isempty(paths1)
   pathname = paths1{1};
end

[pathname, filename, ext] = fileparts(pathname);
filename=[filename ext];
iscurdir = strcmp(pathname, pwd);

function paths1 = findinsubdirs(adir, afile,paths)
dddir = dir(adir);
ss = dir([adir filesep afile]);
ss = {ss.name};
for i=1:length(ss)
   if isempty(strfind(ss{i},filesep))
      ss{i}=[adir filesep ss{i}];
   end

end

paths1 = [paths ss];
for i = 1:length(dddir)
   if dddir(i).isdir && ~strncmp(dddir(i).name,'...',length(dddir(i).name))
      paths1 = findinsubdirs([adir filesep dddir(i).name], afile, paths1);
   end

end


