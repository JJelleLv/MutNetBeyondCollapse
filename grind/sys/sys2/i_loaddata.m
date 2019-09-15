function [varlist, amatrix] = i_loaddata(filename)
global g_grind;
fid = fopen(filename, 'r');
if (fid == -1)
   error('GRIND:loaddata:NoFile','File not found or permission denied.');
end
try
   [g_grind.loaddata.path, g_grind.loaddata.name, g_grind.loaddata.ext] = fileparts(filename);
   TAB = sprintf('\t');
   firstline = fgetl(fid);
   if strcmp(firstline, '%model') %inifile
      while ~feof(fid) && ~strcmp(firstline, '%[data]')
         firstline = fgetl(fid);
      end

      if feof(fid)
         varlist = {};
         amatrix = [];
         warning('GRIND:loaddata:EmptyIniFile','The inifile "%s" has no data', filename);
         return;
      end

      firstline = fgetl(fid);
   end

   isTAB = 0;
   if strcontains(firstline, TAB)
      dlm = TAB;
      isTAB = 1;
   elseif strcontains(firstline, ',')
      dlm = ',';
   elseif strcontains(firstline, ' ')
      dlm = ' ';
   else
      error('GRIND:loaddata:FileError','File format incorrect: no delimiter found in the file');
   end

   if ~isTAB
      firstline = strrep(firstline, dlm, TAB);
   end

   amatrix = cell(2000, 1);
   k = 1;
   firstline= strrep(firstline,' ','_');
   f = strfind(firstline, dlm);
   if max(isletter(firstline)) %if firstline in A ...Z or a..z
      varlist = strrep(firstline,sprintf('\t'),sprintf('\n'));  %#ok<SPRINTFN>
      amatrix = {};
   else
      varlist = sprintf('Column_%d\n', 1:length(f) + 1);
      amatrix{k} = firstline;
      k = k + 1;
   end

   varlist = transpose(str2cell(varlist));
   changeddec = 0;
   while 1
      line = fgetl(fid);
      if ~ischar(line)
         break
      else
         if ~isTAB
            line = strrep(line, dlm, TAB);
         elseif strcontains(line, ',') %the user made a TAB delimited file with decimal commas
            line=strrep(line,',','.');
            changeddec = 1;
         end

         f1 = strfind(line, TAB); %fill missing tabs at the end of lines (Excell problem)
         if length(f1) < length(f)
            line = [line ones(1, length(f) - length(f1)) .* TAB]; 
         end

         line = addnans(TAB, line);
         amatrix{k} = line;
         k = k + 1;
         % amatrix = [amatrix, {line}]; 
      end

   end
   amatrix = amatrix(1:k - 1);
   if changeddec
      warning('GRIND:loaddata:commas','Decimal commas are changed into dots');
   end

   fclose(fid);
catch err
   fclose(fid);
   rethrow(err);
end

function line = addnans(char, line)
%char=sprintf('\t');
f = strfind(line, char);
lenline = length(line);
df = [ diff(f) 0];
dubf=f(df == 1);
for i = length(dubf):-1:1
   line = sprintf('%sNaN%s', line(1:dubf(i)), line(dubf(i) + 1:end));
end

if f(1) == 1
   line = ['NaN' line];
end

if f(end) == lenline
   line = [ line 'NaN'];
end

