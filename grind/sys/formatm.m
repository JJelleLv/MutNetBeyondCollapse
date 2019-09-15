%FORMATM - set the inendation of a m file
%   This function can mess up the inendation, therefore
%    always check the results before approving
function formatm(afile,incfunct)
if nargin == 0
   afile = uigetfile;
end
if nargin <2
    incfunct=0;
else
    if char(incfunct)
        incfunct=str2num(incfunct);  %#ok<ST2NM>
    end
end
if strcmp(afile, '1')
   d = what(pwd);
   d.m = sort(lower(d.m));
   for i = 1:length(d.m)
      formatm(char(d.m(i)));
   end
   return;
end
afile = which(afile);
if isempty(afile)
    error('formatm:file','File not found');
end
nincrem = 3;
fid = fopen(afile, 'r');
incr = 0;
line = '';
nmax = 100;
k = 1;
lines = cell(nmax, 1);
iscomment = 0;
while ~feof(fid)
   f=strfind(line,'...');
   if ~iscomment && ~isempty(f) && (f(end)==length(line) - 2)
      tmpinc = nincrem;
   else
      tmpinc = 0;
   end
   line = strtrim(fgetl(fid));
   %    if ~ischar(line), break, end
   if ~isempty(line)
      iscomment = (line(1) == '%');
   else
      iscomment = 0;
   end
   if ~iscomment
      if strcmp(line, 'end') || strcmp(line, 'end;')
         incr = incr  - nincrem;
      end
      if strcmp(line, 'catch') || strncmp(line, 'catch ', 6) || strcmp(line, 'else') || strncmp(line, 'elseif ', 7) || strncmp(line, 'case ', 5)
         tmpinc = -nincrem;
      end
      if strncmp(line, 'case ', 5) || strcmp(line, 'otherwise')
         tmpinc = -round(nincrem / 2);
      end
      if incr < 0
         incr = 0;
      end
      %repair bug
      line=strrep(line,'./','./');
      line=strrep(line,'.*','.*');
      
      %spacing (not done in lines with >1 strings)
      line=checkspace(line, ' ',{'-'},' ',{'= -','=-','E-','e-','* -','+ -','(-','( -',': -',':-'});
      line=checkspace(line, ' ', {'*'},' ',{'.*'});
      line=checkspace(line, ' ', {'/'},' ',{'./'});
      line=checkspace(line, ' ', {'\'},' ',{'.\'});
      line=checkspace(line, ' ', {'~=', '>=', '<=', '+', '.*', '==', '.\','|', '&', './'}, ' ',{'|','&'}); %No  convert
      line=checkspace(line, ' ', {'&&', '||'}, ' ', {'&','|'}); %No convert
      line=checkspace(line, ' ', {'&', '|'}, ' ', {'&&', '||'}); %No convert
      line=checkspace(line, ' ', {'='},' ',{'~=','>=','<=','=='});
      line=checkspace(line, ' ', {'>'},' ',{'>='});
      line=checkspace(line, ' ', {'<'},' ',{'<='});
      line = checkspace(line, '', {',', ';'}, ' ', {});
      line = checkspace(line, ' ', {'...'}, '', {});
   end
   lines{k} = [char(ones(1, incr + tmpinc) * ' ') line];
   k = k + 1;
   if k > nmax
      nmax = nmax + 100;
      lines = [lines; cell(100, 1)]; 
   end
   if incfunct&&~iscomment && ~strcontains(line,'; end') && ~strcontains( line,', end') && ... 
       strncmp(line, 'function ', 9) 
      incr = incr + nincrem;
   end
   if ~iscomment && ~strcontains(line,'; end') && ~strcontains( line,', end') && ... 
      (strncmp(line, 'if ', 3) || strncmp(line, 'while ', 6) || ...
      strncmp(line, 'for ', 4)) || strncmp(line, 'switch ', 7)...
      || strncmp(line, 'switch(', 7) || strcmp(line, 'try')||...
     strncmp(line, 'classdef ', 9)||strcmp(line, 'properties')||...
     strncmp(line, 'properties ',11)||strcmp(line, 'methods')||...
     strncmp(line, 'methods ',8)   
      incr = incr + nincrem;
   end
end
fclose(fid);
for i = 1:k - 1
   disp(lines{i});
end
disp(' ');
disp('Please check the above file');
overwrite = input('Press 1 if you want to overwrite the file > ');
if overwrite
   disp(['overwrite ' afile]);
   LF = sprintf('\n');  %#ok<SPRINTFN>
   fid = fopen(afile, 'w');
   for i = 1:k - 1
      fwrite(fid,[lines{i}, LF]);
   end
   fclose(fid);
end

function line2 = checkspace(line, spbefor, keywords, spafter, exclkw)
line2 = line;
fparen = strfind(line2,'''');
nparen = length(fparen);
if nparen < 3
   for i = 1:length(keywords)
      kw2 = [char(spbefor) keywords{i} char(spafter)];
      fkw = strfind(line2, keywords{i});
      if ~isempty(fkw)
         fkw2 = strfind(line2, kw2);
         k = 1;
         while isempty(fkw2) && (k <= length(exclkw))
            fkw2 = strfind(line2, exclkw{k});
            k = k + 1;
         end
         if isempty(fkw2)
            if (nparen <= 1)
               line2 = strtrim(strrep(line2, keywords{i}, kw2));
            else
               ok = 1;
               for j = 1:length(fkw)
                  if (fkw(j) > fparen(1)) && (fkw(j) < fparen(2))
                     ok = 0;
                  end
               end
               if ok
                  line2 = strtrim(strrep(line2, keywords{i}, kw2));
               end
            end
         end
      end
   end
end

