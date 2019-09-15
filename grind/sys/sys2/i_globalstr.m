function apar=i_globalstr(AList,AVar,astart)
if nargin>1
   AList=AList(~strcmp(AList,AVar));
end

for i = 1:length(AList)
   f=strfind(AList{i},'(');
   if ~isempty(f) 
       if (~isempty(strfind(AList{i}, '(1)'))||~isempty(strfind(AList{i},'(1,1)')))
          AList{i}= AList{i}(1:f(1)-1);
       else
          AList{i}='';
       end

   end
end
AList=AList(~strcmp(AList,''));
if nargin<3
   apar = sprintf('%s ','global',AList{:});
else
   apar = sprintf('%s ',astart,AList{:});
end

if ~isempty(apar)
    apar(end) = ';';
end

