%true if s1=s2
function [comp, diffields] = struccmp(s1, s2)
comp = true;
diffields={};
if (isstruct(s1)&&~isstruct(s2))||(~isstruct(s1)&&isstruct(s2))
    comp=false;
    return;
elseif isempty(s1)&&isempty(s2)
    return;
elseif ~isstruct(s1)||~isstruct(s2)
    comp=comparefields(s1,s2);
    return;
end

f1 = fieldnames(s1);
f2 = fieldnames(s2);
f3 = intersect(f1,f2);
if length(f1) ~= length(f2) || length(f3)~=length(f1)
    comp = false;
    diffields=setxor(f1,f2)';
end
for i = 1:length(f3)
    v1 =  s1.(f3{i});
    v2 =  s2.(f3{i});
    [comp1,diffsubfields]=comparefields(v1,v2);
    if comp1==0
        comp=false;
        if isempty(diffsubfields)
            diffields=[diffields,f3(i)];
        else
            for j=1:length(diffsubfields)
                diffields=[diffields,{sprintf('%s.%s',f3{i},diffsubfields{j})}];
            end
            
        end
        
    end
    
end




function [comp1,diffsubfields]=comparefields(v1,v2) %compare fields of any class
classv1=class(v1);
sv1=size(v1);
sv2=size(v2);
diffsubfields='';
comp1=true;
if length(sv1)~=length(sv2)||(any(sv1 ~= sv2))
    comp1 = false;
elseif ~strcmp(classv1,class(v2))
    comp1=false;
else
    switch classv1
        case 'char'
            comp1=strcmp(v1,v2);
        case 'cell'
            for j=1:length(v1)
                if ~comparefields(v1{j},v2{j}) %recursion
                    comp1=false;
                    break;
                end

            end

        case 'logical'
            comp1=all(v1==v2);
        case 'struct'
            [comp1,diffsubfields]=struccmp(v1,v2); %recursion on other level
        case 'function_handle'
            comp1=strcmp(char(v1),char(v2));
        otherwise
            if isnumeric(v1)
                ff=v1 ~= v2&~(isnan(v1)&isnan(v2));
                if any(ff(:))
                    comp1 = 0;
                end

            end

    end

end

