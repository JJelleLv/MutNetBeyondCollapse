function i_adjustlabels(hfig,~,htag)
if nargin==0
    hfig=gcf;
end

if nargin<3||isempty(htag)
    h=findobj(hfig,'type','text');
    par=get(h,'parent');
    if ~iscell(par)
        par={par};
    end

    for j=length(h):-1:1
        if ~strcmpi(get(par{j},'type'),'axes')||strcmp(get(par{j},'tag'),'legend')
            h(j)=[];
        end

    end

else
    h=findobj(hfig,'tag',htag);
end

allowedoverlap=0.15;
ud=get(hfig,'userdata');
if isstruct(ud)&&isfield(ud,'allowedoverlap')
    allowedoverlap=ud.allowedoverlap;
end
for j=1:length(h)
    if strcmp(get(h(j),'units'),'data')
        ud=get(h(j),'userdata');
        if isfield(ud,'textcoordinates')
            set(h(j),'position',ud.textcoordinates);
        else
            data=get(h(j),'position');
            ud.textcoordinates=data;
            set(h(j),'userdata',ud);
        end

    else
        set(h(j),'units','data')
        ud=get(h(j),'userdata');
        if isfield(ud,'textcoordinates')
            set(h(j),'position',ud.textcoordinates);
        end

    end
end


set(h,'units','pixel');
set(h,'verticalAlign','bottom','horizontalAlign','left');
%default: text at NE location
set(h,'visible','off')

%remove labels with the same name
strngs=get(h,'string');
if iscell(strngs)
   [~,ndx]=unique(strngs);
   huni=h(ndx);
else
   huni=h;
end
set(huni,'visible','on');
extents=get(huni,'extent');
if ~iscell(extents)
    extents={extents};
end
for i=2:length(huni)
    for j=1:i-1
        if overlap(extents{i},extents{j})>allowedoverlap
            %try SW location
            set(huni(i),'verticalAlign','top','horizontalAlign','right');
            extents{i}=get(huni(i),'extent');
            for k=1:j
                if overlap(extents{i},extents{k})>allowedoverlap
                    set(huni(i),'visible','off');
                    break;
                end
            end
        end
    end
end

for j=1:length(h)
    if ~strcmp(get(h(j),'units'),'data')
        set(h(j),'units','data')
        ud=get(h(j),'userdata');
        if isfield(ud,'textcoordinates')
            set(h(j),'position',ud.textcoordinates);
        end

    end
end


function res=overlap(r1,r2)
%e1,e2=Position and size of text. A four-element vector that defines the size and position of the text string:
%[left,bottom,width,height]
%left and bottom are the distance from the lower left corner 
XA1=r1(1);
XA2=r1(1)+r1(3);
YA1=r1(2);
YA2=r1(2)+r1(4);
XB1=r2(1);
XB2=r2(1)+r2(3);
YB1=r2(2);
YB2=r2(2)+r2(4);
S=(XA2-XA1)*(YA2-YA1)+(XB2-XB1)*(YB2-YB1);
SI= 2*max(0, min(XA2, XB2) - max(XA1, XB1)) * max(0, min(YA2, YB2) - max(YA1, YB1));
res=SI/S;

