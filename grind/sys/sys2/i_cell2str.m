function s = i_cell2str(acell)
if isempty(acell)
   s = '{}';
else
   isstrcell = 0;
   isnumcell = 0;
   for i = 1:length(acell)
      if isempty(acell{i})
         isstrcell = 1;
         isnumcell = 1;
         break;
      elseif ischar(acell{i})
         isstrcell = 1;
      else
         isnumcell = 1;
      end

   end

   if isstrcell && ~isnumcell
      s = sprintf('''%s'' ',acell{:});
      s = sprintf('{%s}', s(1:end - 1));
   elseif isnumcell && ~isstrcell
      s = sprintf('%g ', acell{:});
      s = sprintf('{%s}', s(1:end - 1));
   else
      s = '{';
      for i = 1:length(acell)
         if ischar(acell{i})
            s=sprintf('%s''%s'' ',s,acell{i});
         else
            if isempty(acell{i})
               s =  sprintf('%s[] ', s);
            else
               if ~isempty(acell{i})
                  s=sprintf('%s[%s] ',s,sprintf('%g ',acell{i}));
               else
                  s = sprintf('%s%g ', s, acell{i});
               end

            end

         end

      end

      s = sprintf('%s}', s(1:end - 1));
   end

end

