%function to get the cellnumber of a matrix
function res=cellno(A)
global g_t;
if size(A,1)==size(g_t,1)
   res=repmat((1:size(A,2)),size(A,1),1);
else
   res=1:length(A(:));
end
