function anum=i_checkstr(astr)
if isempty(astr)
   anum=[];
   return;
end
if ischar(astr)
  anum=str2num(astr);  %#ok<ST2NM>
  if isempty(anum)||~(isnumeric(anum)||islogical(anum))
     try
       anum=evalin('base',astr);
     catch %#ok<CTCH>
       anum=[];
     end
  end
else
   anum=astr;
end




   
