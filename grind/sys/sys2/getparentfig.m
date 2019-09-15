function hfig=getparentfig(hobj)
if ~ishandle(hobj)
    error('grind:getparentfig','Argument should be a valid handle');
end
for i=1:10 %max 10 levels
    if strcmp(get(hobj,'type'),'figure')
        hfig=hobj;
        return
    elseif isempty(hobj)||~ishandle(hobj)
        hfig=gcf;
        if isempty(hfig)
            hfig=gcbf;
        end           
        return;
    else
        hobj=get(hobj,'parent');
    end
end
error('grind:getparentfig','Cannot find parent figure');

