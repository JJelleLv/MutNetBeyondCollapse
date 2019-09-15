function res=i_checklog(s)
res=zeros(1,2);
if ischar(s)
   if strcmpi(s,'logx')
      res(1)=1;
   elseif strcmpi(s,'logy')
      res(2)=1;
   elseif strcmpi(s,'logxy')||strcmpi(s,'logyx')
      res=[1 1];
   else
      res=[];
   end

else
   res=[];
end

