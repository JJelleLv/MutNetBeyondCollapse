%USE   Use model
%   Opens an inifile and use the model. If an wildcard (asterisk) is used,
%   a list of possibilities is displayed
%  
%   Usage:
%   USE opens the model dialog box (similar as MODEL).
%   USE FILE - Opens FILE if in current directory, else it searches for the file.
%   USE AMOD* - searches for files using wildcards * and ? AMOD*.
%   USE('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'file' [string] - name of the inifile
%   USE('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' - close the current model and clear memory
%
%   See also model
%
%   Reference page in Help browser:
%      <a href="matlab:commands('use')">commands use</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function use(varargin)
%amodel = [];
%if the model name has spaces it can be read as separate arguments
fieldnams={'file', 's', 'name of the inifile',''}';
args=i_parseargs(fieldnams,'file(+)',{'-c'},varargin);
if isfield(args,'file')&&iscell(args.file)
    args.file=strtrim(sprintf('%s ', args.file{:}));
end
if ~isfield(args,'file')
    args.file=[];
end
if isempty(varargin)
   model
   return;
elseif  any(strcmp(args.opts,'-c'))
   model -clear;
   return;
end

finishgrind;
if ~exist('i_use','file')
  addpath([grindpath filesep 'sys2']);
end
i_use(args.file,1,1,1,1);
par('modelonly');
