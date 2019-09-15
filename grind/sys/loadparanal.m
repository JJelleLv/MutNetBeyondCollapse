function loadparanal(varargin)
global g_grind; 
fieldnams={'file', 's', 'file to load',''}';
args=i_parseargs(fieldnams,'file','',varargin);
if ~isfield(args,'file')
     [args.file, pathname] = uigetfile('*.mat', 'Load paranal file');
     args.file=[pathname args.file];
end
inif=g_grind.inifile;
load(args.file);
if ~strcmp(inif,g_grind.inifile)
  use(g_grind.inifile);
%  global g_paranal g_grind;
  load(args.file);
end
paranal 1
