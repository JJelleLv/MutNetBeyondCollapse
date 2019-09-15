% VARPASTE- Paste any matrix from clipboard
%   Paste any matrix from the Windows clipboard, which was stored 
%   as text tab delimited (for instance with Microsoft Excel). At default 
%   if there are text strings in the clipboard a cell array (values/strings) 
%   is created, else a numerical matrix.
%
%   Usage:
%   A=VARPASTE - Paste the contents of the clipboard to the variable A (cell array or matrix).
%   A=VARPASTE(CLASS) - Paste the contents of the clipboard to the variable A, force A 
%   to be of a class CLASS.
%   A=VARPASTE('-file=filename', CLASS) - Load the contents of a tab-delimited text file to the variable A, force A 
%   to be of a class CLASS.
%   The following classes are supported:
%     'numeric' - matrix of numbers (default if there are no missing values)
%     'cell' - cell matrix with numbers and strings (default if there are strings)
%     'string' - cell list with one string per line
%     'logical' - matrix with logical values
%     'dataset' - dataset, see Statistics Toolbox.
%     'table' - table (matlab versions > 2014).
%     'char' - character matrix
%     'int8','uint8','int16','uint16','int32','uint32','int64','uint64'-  integer matrix
%     'single' - single precision matrix
%     'double' - double precision matrix (the same as 'numeric')
%    VARPASTE('argname',argvalue,...) - Valid argument name-value pairs [with type]:
%     'class' [numeric | cell | char | string | dataset | table | logical | int8 | uint8 | int16 | uint16 | int32 | uint32 | int64 | uint64 | single | double or empty] - the class of the resulting variable (default='numeric')
%     'delimiter' [string] - the delimiter for columns (default \t or whitespace)
%     'whitespace' [string] - the whitespace character. Repeating whitespaces are ignored (default '')
%   VARPASTE('-opt1','-opt2',...) - Valid command line options:
%     '-f=filename' - use the file filename instead of the clipboard
%
%
%See also varcopy, <a href="matlab:help class">class</a>
%
%   Reference page in Help browser:
%      <a href="matlab:commands('varpaste')">commands varpaste</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function A = varpaste(varargin)%forcetype, filename)

fieldnams={'class', 'e[n+umeric|ce+ll|ch+ar|string|da+taset|t+able|l+ogical|int8|uint8|int16|uint16|int32|uint32|int64|uint64|si+ngle|do+uble]#E', 'the class of the resulting variable (default=''numeric'')';...
   'delimiter', 's', 'the delimiter for columns (default \t or whitespace)';...
   'whitespace', 's', 'the whitespace character. Repeating whitespaces are ignored (default '''')'}';
args=i_parseargs(fieldnams,'class', '-f',varargin);
if ~isfield(args,'delimiter')
   if isfield(args,'whitespace')
       args.delimiter=args.whitespace;
   else
       args.delimiter='\t';
   end
end
if ~isfield(args,'whitespace')
   args.whitespace='';
end
f1=strncmp(args.opts,'-f',2);
if any(f1)
   f=strfind(args.opts{f1}, '=');
   if ~isempty(f)
      if isfield(args,'filename')
         h = args.filename;
      end
      args.filename = args.opts{f1}(f(1) + 1:end);
%       if nargin == 2
%          forcetype = h;
%       else
%          forcetype = '';
%       end
   else
%       if nargin == 2
%          h = filename;
%       end
      [filename,pathname]=uigetfile('*.*','What is the name of the file');
      args.filename = [pathname filename];
%       if nargin == 2
%          forcetype = h;
%       end
   end
else
   args.filename = 'clipboard';
end
if ~isfield(args,'class')
   args.class = '';
end
% if ~isempty(forcetype) && (forcetype(1)=='-')
%    forcetype = forcetype(2:end);
% end


forcetype2 = args.class;
if any(strcmp(args.class,{'int8','uint8','int16','uint16','int32','uint32','int64',...
      'uint64','single','double'}))
   args.class = 'numeric';
end
% if ~any(strcmp(args.class,{'','numeric','cell','char','string','dataset','table','logical','int8','uint8','int16',...
%       'uint16','int32','uint32','int64',...
%       'uint64','single','double'}))
%    error('GRIND:varpaste','Unsupported class "%s" for conversion',args.class);
% end

if strcmpi(args.filename, 'clipboard')
   lines = str2cell(clipboard('paste'));
else
   f=strfind(args.filename, '=');
   if ~isempty(f)
      args.filename = args.filename(f(1) + 1:end);
   end
   fid = fopen(args.filename, 'r');
   try
      %       i = 1;
      %       while ~feof(fid)
      %          lines{i} = fgetl(fid);
      %          i = i + 1;
      %       end
      try
        lines=textscan(fid,'%s', 'delimiter',args.delimiter, 'whitespace',args.whitespace);
      catch err
        if ~isempty(strfind(err.identifier,'BufferOverflow'))
            lines=textscan(fid,'%s', 'delimiter',args.delimiter, 'whitespace',args.delimiter, 'bufsize',1000000); %#ok<BUFSIZE>
        else
            rethrow(err);
        end
      end
%       if  verLessThan('matlab','8')
%           lines=textscan(fid,'%s','delimiter','\n','bufsize',100000); %much faster
%       else
%           lines=textscan(fid,'%s','delimiter','\n'); %much faster
%       end
      lines = lines{1};
      fclose(fid);
   catch err
      fclose(fid);
      rethrow(err)
   end
end
if strncmpi(args.class, 's', 1) %string
   A = lines;
   return;
end
if strncmpi(args.class, 'ch', 2) %char
   A = char(lines);
   return;
end
if ~isempty(lines)
   line = lines{1};
else
   line='';
end
if ~isempty(line)&&((line(1)=='[') || (line(1)=='{'))
   line = sprintf('%s', lines{:});
   try
      A = eval(line);
   catch
      A = stabread(lines,args.delimiter,args.whitespace);
   end
else
   A = stabread(lines,args.delimiter,args.whitespace);
end
if strncmpi(args.class, 'n', 1) && iscell(A) %numeric
   A = cell2num(A);
end
if strncmpi(args.class, 'c', 1) && isnumeric(A) %cell
   A = num2cell(A);
end
if strncmpi(args.class, 'l', 1) %logical
   M = false(size(A));
   hasnan = false;
   if iscell(A)
      for i = 1:size(A, 1)
         for j = 1:size(A, 2)
            if (isnumeric(A{i, j}) || islogical(A{i, j}))
               if isnan(A{i, j})
                  hasnan = true;
               else
                  M(i, j) = logical(A{i, j});
               end
            elseif any(strcmpi(A{i,j},{'on','yes'}))
               M(i, j) = true;
            end
         end
      end
   end
   if isnumeric(A)
      nans = isnan(A);
      if any(nans)
         M(nans) = 0;
         M(~nans) = logical(A(~nans));
         hasnan = 1;
      else
         M = logical(A);
      end
   end
   if hasnan
      warning('GRIND:varpaste:logicalNaN','Logical NaN are not defined: NaNs converted to false');
   end
   A = M;
end
if strncmpi(args.class, 'd', 1) ||strncmpi(args.class, 't', 1)%dataset
   if strncmpi(args.class, 't', 1)
      ds = table;
   else
      ds = dataset;
   end
   h = 1;
   fields = cell(size(A, 2), 1);
%   rownames1 = 0;
   for i = 1:size(A, 2)
      if iscell(A)&&~isempty(A{1, i}) && ischar(A{1, i})
         h = 2;
         afield = strtrim(A{1, i});
         J = str2num(afield); %#ok<ST2NM>
         if ~isnumeric(J) %matlab bug
             J=[];
         end
         if (~isempty(J) && ~isnan(J))
            afield = sprintf('field%d', i);
         end
         spattern={'"','';...
             '+','p';...
             '-','m';...
             '^[^A-Za-z]','V';...
             '[^A-Za-z0-9_]','_'};
         afield=regexprep(afield,spattern(:,1),spattern(:,2));
%          if i == 1
%             rownames1 = 1;
%          end
         %afield = sprintf('field%d', i);
      end
      fields{i} = afield;
   end
   for i = 1:size(A, 2)
      if iscell(A)
         AA = cell2num(A(h:end, i));
      else
         AA = A(h:end, i);
      end
      if all(isnan(AA))
         AA = A(2:end, i);
      end
      ds.(fields{i}) = AA;
   end
%    if rownames1&&iscell(ds.field1)
%       if isa(ds, 'dataset')
%          ds.Properties.ObsNames = ds.field1;
%       else
%          ds.Properties.RowNames = ds.field1;
%       end
%       ds.field1 = [];
%    end
   A = ds;
end
if ~strcmp(args.class, forcetype2)
   A =  cast(A, forcetype2);
end


function A = stabread(lines,delim,whitespace)
hasstr = 0;
delim=sprintf(delim);
if isempty(lines)
    A=[];
else
    if ~isempty(whitespace)
       %remove double spaces
       if strcmp(delim,whitespace)
          lines= regexprep(lines,'[ ]*',' ');
       else
          lines= regexprep(lines,'[ ]*','');
       end
    end
    f = strfind(lines{1}, delim);
    A = cell(length(lines), length(f) + 1);
    for j = 1:length(lines)
        line = lines{j};
        if j == 1
            if ~strcontains(line, delim)
                delim = ',';
                if ~strcontains(line, delim)
                    delim = ';';
                end
            end
        end
        if ~isempty(line)
            [aa, ~, err] = sscanf(line, '%f');
            if isempty(err)&&length(aa)==size(A, 2)
                if isempty(aa)
                    A = lines;
                    hasstr = 1;
                else
                    for i = 1:length(aa)
                        A{j, i} = aa(i);
                    end
                end
            else
                f = [0,strfind(line, delim),length(line) + 1];
                for i = 1:length(f) - 1
                    s = line(f(i) + 1:f(i + 1) - 1);
                    J = str2double(s);
                    if length(strfind(s, '-')) > 1
                        J = [];
                    end
                    if  (~isempty(J) && ~isnan(J))
                        res = J;
                    else
                        res = s;
                        hasstr = 1;
                    end
                    A{j, i} = res;
                end
            end
        end
    end
    if ~hasstr
        A = cell2num(A);
    end
end


