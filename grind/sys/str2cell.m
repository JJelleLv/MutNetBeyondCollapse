%STR2CELL - Convert a string with end-of-line marks to a cell structure
%
%    Example:
%    c=str2cell(sprintf('line 1\n line 2\ntest'))
%    c =
%         'line 1'
%         ' line 2'
%         'test'
%    See also SPRINTF
function c = str2cell(s)

if isempty(s)
    c={};
elseif size(s, 1) > 1
   c=cellstr(s);
%    c=cell(size(s,1),1);
%    for i=1:size(s,1)
%       c{i}=char(s(i,:));
%    end
%    c=strtrim(c);
else
%     try
%         c=textscan(s,'%s', 'delimiter','\n', 'whitespace','');
%     catch err
%         if ~isempty(strfind(err.identifier,'BufferOverflow'))
%             c=textscan(s,'%s', 'delimiter','\n', 'whitespace','', 'bufsize',10000000);
%         else
%             rethrow(err);
%         end
%     end
%     c=c{1};
    % Maybe slightly slower but more robust than textscan
    c=transpose(regexp(s,'\n','split'));
    if isempty(c{end})  %trailing \n ?
        c=c(1:end-1);
    end
end
%    f = strfind(s, sprintf('\n'));
%    if isempty(f)
%       c{1} = s;
%    else
%       c = cell(length(f), 1);
%       p = 1;
%       for i = 1:length(f)
%          c{i} = s(p:f(i) - 1);
%          p = f(i) + 1;
%       end
%       if p <= length(s)
%          c{length(f) + 1} = s(p:length(s));
%       end
%    end
% end
