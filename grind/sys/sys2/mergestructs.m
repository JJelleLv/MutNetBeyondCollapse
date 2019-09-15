%merge two structs
%astruct=mergestructs(astruct,newstruct)
%
%fields of a astruct will be replaced by newstruct, new fields will be
%added if necessary
%nested for substructures
function fullstruct=mergestructs(fullstruct,substruct)
if isempty(fullstruct)
    fullstruct=substruct;
elseif ~isempty(substruct)&&numel(substruct)==numel(fullstruct)
    fnew=fieldnames(substruct);
    fold=fieldnames(fullstruct);
    for i=1:length(fnew)
        for j=1:numel(substruct)
            if isstruct(substruct(j).(fnew{i}))&&any(strcmp(fold,fnew{i}))&&isstruct(fullstruct(j).(fnew{i}))
                fullstruct(j).(fnew{i})=mergestructs(fullstruct(j).(fnew{i}),substruct(j).(fnew{i}));
            else
                fullstruct(j).(fnew{i})=substruct(j).(fnew{i});
            end

        end

    end

else
    error('grind:mergestructs','Cannot merge strucs of different sizes')
end

