function A1=cell2num(A)
   A1 = zeros(size(A));
   for i = 1:size(A, 2)
      for j = 1:size(A, 1)
         aa=A{j, i};
         if islogical(aa) || isnumeric(aa)
            if isempty(aa)
               A1(j,i)=NaN;
            else
               A1(j,i)= aa;
            end
         else
            A1(j, i) = str2double(aa);
         end
      end
   end
