%LOADDATA   Load data for external variables or state variables
%   Load a data file with data of external variables or state variables (for 
%   data fitting, see <a href="matlab:help optimpars">optimpars</a>). The first row of the file should contain 
%   the names of the variables (t=time). If there is no time entered, it is 
%   assumed that the data are equally spaced with time step 1. The format of the
%   ASCI file should be TAB delimited, comma delimited or space delimited. You can
%   also save the data to a MAT file which should contain an array with the data and
%   a cell array with the variable names (or use table or dataset instead).
%   It is also possible to read data from an ini file (see <a href="matlab:help savepar">savepar</a>)
%
%   Usage:
%   LOADDATA - The user can select a file.
%   LOADDATA FILE- Loads the data file named FILE.
%   LOADDATA('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'file' [file name] - name of the data file.
%     'varlist' [cell] - names of the variables in the columns of the file
%   LOADDATA('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-ini' - read data from the current ini file.
%     '-s' - Skips columns with names that do not exist in the model.
%
%
%   See also setdata, defextern, optimpars, savepar
%
%   Reference page in Help browser:
%      <a href="matlab:commands('loaddata')">commands loaddata</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function loaddata(varargin)
 global g_grind;
if nargin == 0
   if isfield(g_grind, 'loaddata')
      filter=[g_grind.loaddata.path filesep '*.csv' pathsep '*.dat' pathsep '*.txt' pathsep '*.ini' pathsep '*.mat'];
   else
      filter=['*.csv' pathsep '*.dat' pathsep '*.txt' pathsep '*.ini' pathsep '*.mat'];
   end
   [filename,path] = uigetfile(filter, 'Get data file (delimited or CSV)');
   if filename ~= 0
      args.file = [path filename];
   else
      return;
   end
else
    fieldnams={'file', 'F', 'name of the data file.','';...
   'varlist', 'c', 'names of the variables in the columns of the file',''}';
    args=i_parseargs(fieldnams,'file','-ini,-s',varargin);
end

if any(strcmp(args.opts, '-ini'))
   args.file = g_grind.inifile;
end
if ~isempty(strfind(args.file, ';'))
   args.file  = args.file(1:strfind(args.file, ';') - 1);
end
[~, ~, ext] = fileparts(args.file);
if strcmp(ext, '.mat')
   [amatrix, varlist] = loadmat(args.file);
else
   [varlist, amatrix] = i_loaddata(args.file);
end
if ischar(varlist)
   varlist = str2cell(varlist);
end
if isfield(args,'varlist')
    varlist=args.varlist;
end
if any(strcmp(args.opts, '-s'))  %skip variables that are unknown
   for i = 1:length(varlist)
      if ~strcmp(varlist{i}, 't')
         p = i_getno(varlist{i});
         if isempty(p.no)||p.ispar
            varlist{i} = '###skip###';
         end
      end
   end
end
if ~(isempty(varlist)||isempty(amatrix))
   setdata(amatrix, varlist);
else
   error('grind:loaddata','Error reading file');
end
function  [g_amatrix, g_varlist] = loadmat(g_filename)
%g_ before names to avoid conflicts
g_amatrix = [];
g_varlist = {};
d = whos('-file', g_filename);
for i = 1:length(d)
   d1 = d(i);
   if strcmp(d1.class, 'double')
      load(g_filename,d1.name) %as late as possible
      g_amatrix = eval(d1.name);
   elseif strcmp(d1.class, 'cell')
      load(g_filename,d1.name)
      g_varlist = eval(d1.name);
   elseif strcmp(d1.class, 'dataset') %statistics toolbox
      if ~i_hastoolbox('stats')
            error('grind:loaddata','Statistics toolbox needed for reading this file: "%s"',g_filename);
      end
      load(g_filename,d1.name)
      g_i_dat = eval(d1.name);
      g_amatrix = double(g_i_dat);
      g_varlist = g_i_dat.Properties.VarNames;
   elseif strcmp(d1.class, 'table')
      load(g_filename,d1.name)
      g_i_dat = eval(d1.name);
      g_amatrix = g_i_dat{:, :};
      g_varlist = g_i_dat.Properties.VariableNames;
   end
end
