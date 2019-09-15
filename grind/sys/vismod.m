%VISMOD   Create a model as a Forrester diagram
%
%
%   Select different elements of a diagram to create a Forrester diagram.
%   Elements are:
%  *<a href="matlab:commands concept_statevar">State variables</a>, name and initial values of state variables.
%  *<a href="matlab:commands concept_auxvar">Auxiliary variables</a>, name and equation of auxiliary (or help variables.
%   Such variable contains an equation.
%  *<a href="matlab:commands concept_parameter">Parameters</a>, name and values of parameters.
%  *<a href="matlab:commands concept_externvar">External variables</a>, name and values of external variables.
%  *<a href="matlab:commands concept_externfunc">External or used-defined functions</a>, name of external function.
%  *<a href="matlab:commands concept_flow">Continuous flow of substance</a>, a part of a differential equation.
%  *<a href="matlab:commands concept_train">Discrete flow of substance</a>, a part of a difference equation.
%  *<a href="matlab:commands concept_connector">Connectors</a>, information flow between components.
%  *<a href="matlab:commands concept_cloud">Cloud</a>, a source or sink
%
%   Usage:
%   VISMOD - create new model or edit the current model.
%   VISMOD FILE - edit inifile FILE.
%   VISMOD('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'file' [string] - name of the file
%   VISMOD('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' - Clears the current model from memory (and deletes currently temporary files).
%
%
%   Note:
%   VISMOD works only on Microsoft Windows systems.
%  
%  
%   See also model, modelpanel
%
%   Reference page in Help browser:
%      <a href="matlab:commands('vismod')">commands vismod</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function vismod(varargin)
global g_grind;
if ~exist('i_use','file')
  addpath([grindpath filesep 'sys2']);
end
if ~strncmpi(computer,'PCWIN',5)
    error('GRIND:vismod:NoPCWIN','"vismod" works only on PC systems, please use <a href="matlab:model">model</a> instead to create a model');    
end
fieldnams={'file', 's', 'name of the file',''}';
args=i_parseargs(fieldnams,'file(+)',{'-c'},varargin);

if any(strcmp(args.opts,'-c'))
    model('-c'); % clear the current model
    return;
end
g=grindpath(2);
if isempty(g_grind)&&~strcmp(pwd, g)
   cd(g);
   fprintf('Changed directory to %s\n',g);
end
tmpf = 'visualgrind%d.tmp';
i = 0;
tmpfile = sprintf(tmpf, i);
while exist(fullfile(grindpath, tmpfile), 'file')  && (i < 10)
   i = i + 1;
   tmpfile = sprintf(tmpf, i);
end
if nargin == 0
   if ~isempty(g_grind)&&isfield(g_grind,'inifile')&&~isempty(g_grind.inifile)
      s = g_grind.inifile;
      args.file = sprintf('curfile#%s', g_grind.inifile);
 %     err=1;
 %     try        
 %         i_parcheck;
 %         err=0;
 %     catch 
 %     end
 %     if ~err
       if ~isempty(g_grind.model)
        savepar(fullfile(grindpath, tmpfile), 1, 1);
       elseif ~isempty(g_grind.inifile)
           copyfile(g_grind.inifile,fullfile(grindpath, tmpfile));
       end
 %     end
      g_grind.inifile = s;
   else
      args.file = 'new';
   end
else
   if isempty(strfind(args.file, '.'))
      args.file = [args.file '.ini'];
   end
   [pathname, iscurdir] = findgrindfile(args.file);
   if isempty(pathname)
      i_errordlg(['Inifile: ' args.file ' does''t exist']);
      error('GRIND:vismod:NoInifile','Inifile "%s" doesn''t exist',model);
   elseif ~iscurdir
      disp(['Could not find the model "' args.file '" in "' cd '"']);
      disp(['Changing to:"' pathname '"']);
      cd(pathname);
   end
end
oldpath = cd;
cd(grindpath);
if (length(args.file) > 8) && strcmp(args.file(1:8), 'curfile#')
   name = args.file;
   ext = '';
   p = oldpath;
else
   [p, name, ext] = fileparts(args.file);
end
if isempty(p)
   p = oldpath;
end
s = sprintf('visualgrind vismod "%s" "%s" "%s"', p, [name, ext], tmpfile);
%disp(s);
disp('Close the vismod window to continue working in MATLAB');
dos(s);
cd(grindpath);
fid = fopen(tmpfile, 'r');
if fid > 0
   line = myfgetl(fid);
   if strcmp(line, 'OK')
      oldpath = myfgetl(fid);
      inifile = myfgetl(fid);
      fclose(fid);
      delete(tmpfile);
      cd(oldpath);
      use(inifile);
      newpath=fileparts(inifile);
      if ~isempty(newpath)
         cd(newpath);
      end
      try
         analunits('-c')
      catch err
          analunits('-l')
          fprintf(2,'WARNING: The units of parameters/variables are inconsistent.\nYou can use <a href="matlab:analunits">analunits</a> to change the units, or use <a href="matlab:vismod">vismod</a> to edit the equations.\nIt is however possible to run the model anyway\n\n');
         %disp(err.message)
           rethrow(err);
      end
   else
      disp('Cancelled')
      fclose(fid);
      delete(tmpfile);
      cd(oldpath);
   end
else
   error('GRIND:vismod:commerror','Failed to communicate with visualgrind.exe');
end
function s = myfgetl(fid)
s = fgetl(fid);
if ~isempty(s) && (int16(s(end)) == 13)
   s = s(1:end - 1);
end

