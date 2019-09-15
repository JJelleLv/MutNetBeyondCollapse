%VARCOPY   Copy variable to clipboard
%   Copy a numeric or a cell array (containing only values or character strings) 
%   or a dataset(<a href="matlab:commands toolboxes">statistics toolbox</a>) or a table(version >2014b) to the Windows clipboard, 
%   as text tab delimited.
%
%   Usage:
%   VARCOPY(g_t) - Copy the variable g_t to the clipboard.
%   VARCOPY g_Y - Copy variable g_Y to clipboard.
%   VARCOPY(ds) - Copy dataset ds to clipboard.
%   VARCOPY(tab) - Copy table tab to clipboard.
%
%   See also varpaste
%
%   Reference page in Help browser:
%      <a href="matlab:commands('varcopy')">commands varcopy</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function varcopy(A, precision)
if nargin == 0
   prompt = {'Enter variable:', 'Number of digits to copy:'};
   defaultanswer = {'', '15'};
   answer = inputdlg(prompt, 'varcopy - copy variable to clipboard', 1, defaultanswer);
   A = evalin('base', answer{1});
   precision = str2double(answer{2});
elseif nargin < 2
   precision = 15;
end
if ~ischar(A) && ~isnumeric(A) && ~islogical(A) &&~isa(A, 'dataset')&&~isa(A, 'table') &&~iscell(A)
   error('GRIND:varcopy:typeNotSupp','Variables of type "%s" are not supported', class(A));
end
if ischar(A)
   try
      M = evalin('base', A);
   catch
      M = [];
   end
   if ~isempty(M)
      A = M;
   end
end
if ischar(A)
   clipboard('copy', A);
   return;
end
s = [];
pattrn = sprintf('%%.%dg\\t', precision);
for i = 1:size(A, 1)
   if isnumeric(A) || islogical(A)
      s1 = sprintf(pattrn, A(i, :));
   elseif iscell(A) || isa(A, 'dataset')|| isa(A, 'table')
      if isa(A,'table')&&~isempty(A.Properties.RowNames)
          s1=sprintf('%s\t',A.Properties.RowNames{i});
      elseif isa(A,'dataset')&&~isempty(A.Properties.ObsNames) %statistics toolbox
          s1=sprintf('%s\t',A.Properties.ObsNames{i});
      else
          s1 = [];
      end
      for j = 1:size(A, 2)
         Aij=A{i,j};
         if ischar(Aij)
            s1 = sprintf('%s%s\t', s1, Aij);
         elseif isa(Aij,'function_handle')
            s1 = sprintf('%s@%s\t', s1, func2str(Aij));
         else
            s1 = sprintf(['%s' pattrn], s1, Aij);
         end
      end
   end
   s = sprintf('%s%s\n', s, s1(1:end - 1));
end
if isa(A, 'dataset')
   vars = get(A, 'VarNames');
   if ~isempty(A.Properties.ObsNames)
      s1 = sprintf('%s\t','',vars{:});
   else
      s1 = sprintf('%s\t', vars{:});
   end
   s = sprintf('%s\n%s', s1(1:end - 1), s);
elseif  isa(A, 'table')
   vars = A.Properties.VariableNames;
   if ~isempty(A.Properties.RowNames)
      s1 = sprintf('%s\t','',vars{:});
   else
      s1 = sprintf('%s\t', vars{:});
   end
   s = sprintf('%s\n%s', s1(1:end - 1), s);
end
clipboard('copy', s);



