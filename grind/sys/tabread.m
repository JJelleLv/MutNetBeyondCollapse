% TABREAD - function to read a tab delimited file to a cell or a matrix
%
% Usage: A = tabread(filename) (if there are not numerical parts A becomes a cell otherwise a matrix)
%        A = tabread(filename,1) force that A becomes a matrix (Non-numerical fields become NAN)
%        A = tabread(filename,-1) force that A becomes a cell with strings only)
%        A = tabread(filename,1,[300000,2]) optionally add a size to (greatly) improve the speed (big files!)
function A = tabread(filename, forcematrix ,asize)
if nargin < 2
   forcematrix = 0;
end
if nargin > 2
   if forcematrix == 1
      A = zeros(asize);
   else
      A = cell(asize);
   end
end
j = 0;
hasstr = 0;
fid = fopen(filename, 'r');
if (fid == -1)
   error('GRIND:tabread:NoFile','File not found or permission denied');
end
while ~feof(fid)
   line = fgetl(fid);
   j = j + 1;
   if ~isempty(line)
      if forcematrix == -1
         err = 1;
      else
         [aa, ~, err] = sscanf(line, '%f');
      end
      if isempty(err)
         if forcematrix == 1
            A(j, :) = transpose(aa);
         else
            for i = 1:length(aa)
               A{j, i} = aa(i);
            end
         end
      else
         f = [0,strfind(line, sprintf('\t')),length(line) + 1];
         for i = 1:length(f) - 1
            s = line(f(i) + 1:f(i + 1) - 1);
            if strncmpi(line, '[initial', 8)
               J = [];
            else
               J = str2num(s);  %#ok<ST2NM>
               if ~isnumeric(J)
                   J=[];
               end                   
               %            J = 1; %??????
            end
            if (forcematrix==1) || (~isempty(J) && ~isnan(J) && (forcematrix ~= -1))
               res = J;
            else
               res = s;
               hasstr = 1;
            end
            if forcematrix == 1
               if isempty(res)
                  A(j, i) = NaN;
               else
                  A(j, i) = res;
               end
            else
               A{j, i} = res;
            end
         end
      end
   end
end
if ~hasstr && ~(forcematrix==1)
   A = cell2num(A);
end
fclose(fid);

