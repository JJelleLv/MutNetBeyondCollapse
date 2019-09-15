%f = findstrings(strng, s, wholewords)
%like strfind but can find whole words only
%pattern can be a cell of patterns!
function f = findstrings(pattern, s, wholewords)
if iscell(pattern)
    if wholewords
        spattern=sprintf('(\\<%s\\>)|',pattern{:});
    else
        spattern=sprintf('(%s)|',pattern{:});
    end
    f=regexp(s,spattern,'start');
else
    if wholewords
        f=regexp(s,['\<' pattern '\>'],'start')';
    else
        f=strfind(s,pattern);
    end
end
    
% if wholewords
%     f=strfind(s,pattern);
% end
% if ~iscell(pattern)
%    f = strfind(s, pattern);
%    if ~isempty(f) && wholewords
%       f = checkwhole(pattern, s, f);
%    end
% 
% else
%    f = [];
%    for i = 1:length(pattern)
%       f1 = strfind(s, pattern{i});
%       if ~isempty(f1) && wholewords
%          f1 = checkwhole(pattern{i}, s, f1);
%       end
% 
%       f = [f f1];
%    end
% 
% end
% 
% 
% function f = checkwhole(substr, s, f1)
% len = length(substr);
% f = [];
% for i = 1:length(f1)
%    if ((f1(i)==1) || ~isnumletter(s(f1(i) - 1))) &&  ((f1(i)+len - 1==length(s)) || ~isnumletter(s(f1(i)+len)))
%       f = [f f1(i)];
%    end
% 
% end
% 
% function res = isnumletter(ch)
% res = isletter(ch);
% if ~res
%    res = strcontains('1234567890_', ch);
% end
% 
