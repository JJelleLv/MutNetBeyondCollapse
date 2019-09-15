%my implemntation of the contains function
function res=strcontains(fullstr,pattern)
%contains in Matlab 2017 is slower
if ischar(fullstr)
    res=~isempty(strfind(fullstr,pattern)); %#ok<STREMP>
elseif iscell(fullstr)
    res=~cellfun('isempty',strfind(fullstr,pattern));  %#ok<STRCL1> is faster for chars/cells 
elseif isempty(fullstr)
    res=false;
else
    res=[];
end