function i_resetlabels(hfig,~,htag)
if nargin==0
    hfig=gcf;
end

if nargin<3
    h=findobj(hfig,'type','text');
    par=get(h,'parent');
    for j=length(h):-1:1
        if ~strcmpi(get(par{j},'type'),'axes')||strcmp(get(par{j},'tag'),'legend')
            h(j)=[];
        end

    end

else %more safe to have a special tag
    h=findobj(hfig,'tag',htag);
end

set(h,'visible','on');
% for j=1:length(h)
%     if ~strcmp(get(h(j),'units'),'data')
%         set(h(j),'units','data')
%         ud=get(h(j),'userdata');
%         if isfield(ud,'textcoordinates')
%             set(h(j),'position',ud.textcoordinates);
%         end

%     end
% end



