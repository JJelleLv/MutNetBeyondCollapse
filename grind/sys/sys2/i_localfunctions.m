function [funcs,types]=i_localfunctions(g_model)
%remove comments
g_model=regexp(g_model, '^[^%]*','match','once');
idents=regexp(g_model,'[a-zA-Z_][a-zA-Z0-9_]*[\(]?','match');
if ischar(idents)
    idents={idents};
end

if iscell(idents{1})
  idents=[idents{:}];
end

idents=unique(idents);
ndx=strcontains(idents,'(')|strcmp(idents,'if')|strcmp(idents,'for')|strcmp(idents,'end'); 
funcs=regexp(idents(ndx),'[a-zA-Z0-9_]*','match', 'once');
%funcs=regexp(g_model,'([a-zA-Z_][a-zA-Z0-9_]*(?=[\(]))|(\<if\>)|(\<for\>)|(\<end\>)|(\<while\>)','match');%slower
if isempty(funcs)
    funcs={};
    types={};
    return
end

% if ischar(funcs)
%     funcs={funcs};
% end

% if iscell(funcs{1})
%   funcs=[funcs{:}];
% end


%ndx=cellfun('isempty',strfind(funcs,'%'));
%funcs=unique(funcs);
types=cell(size(funcs));
types(:)={'local'};
for i=1:length(types)
    ext=exist(funcs{i},'builtin');
    if ext==5
        types{i}='built-in';
    else
        pat=which(funcs{i});
        if strncmp(pat,matlabroot,length(matlabroot))
            types{i}='matlab';
        elseif strncmp(pat,grindpath,length(grindpath))
            types{i}='grind';
        end
    end

end

