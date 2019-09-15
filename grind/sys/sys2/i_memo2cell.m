function res = i_memo2cell(avar)
res = avar;
if ~isempty(res)
   if ~iscell(res)
      res = num2cell(res, 2);
   end

   while ~isempty(res) && iscell(res{1})
      res = res{1};
   end

   if size(res, 1) > size(res, 2)
      res = transpose(res);
   end

end


