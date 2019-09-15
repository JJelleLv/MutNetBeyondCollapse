function ndx=i_findparcom(pname)
global g_grind;
%remove comments:
comm=regexp(g_grind.commands, '^[^%]*','match','once');
%check lines where the parameter is directly followed by =
ndx= find(~cellfun('isempty',regexp(comm,sprintf('\\<%s\\>(?=[ ]*[=])',pname))));
if isempty(ndx)
    %if not found allow parameter to have something between the = sign
   ndx= find(~cellfun('isempty',regexp(comm,sprintf('\\<%s\\>(?=.*[=])',pname))));
end

ndx2=strcontains(g_grind.commands(ndx),'%');%~cellfun('isempty',strfind(g_grind.commands(ndx),'%')); 
if any(ndx2)
    ndx=ndx(ndx2);
    ndx=ndx(1);
elseif ~isempty(ndx)
     ndx=ndx(end);
end

% global g_grind;
% ndx = [];
% f = regexp(g_grind.commands, sprintf('\\<%s\\>', pname));
% fndx = find(~cellfun(@isempty, f));
% for k = length(fndx):-1:1
%    f=strfind(char(g_grind.commands{fndx(k)}), '=');
%    if ~isempty(f)
%       s2 = strtrim(g_grind.commands{fndx(k)}(1:f(1) - 1));
%       if strcmp(pname, s2)
%          ndx = fndx(k);
%          return;
%       end

%    end

% end

