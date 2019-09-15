function pos=subplotpos(rows,cols,x1,y1,spacing)
%tight position for subplots, you can only set the spacing between the
%figures (simplified from subaxis)
%   s=subplot('position',subplotpos(rows,cols,x1,y1),'tag',sprintf('subplot(%d,%d)',x1,y1));
%
if nargin<5
    spacing=[0.02 0.02];
end
if length(spacing)==1
    spacing=spacing+zeros(2,1);
end

Args=struct('Holdaxis',0, ...
        'SpacingVertical',spacing(1),'SpacingHorizontal',spacing(2), ...
        'MarginLeft',.2,'MarginRight',.1,'MarginTop',0.05,'MarginBottom',.1, ...
        'rows',rows,'cols',cols); 


cellwidth=((1-Args.MarginLeft-Args.MarginRight)-(Args.cols-1)*Args.SpacingHorizontal)/Args.cols;
cellheight=((1-Args.MarginTop-Args.MarginBottom)-(Args.rows-1)*Args.SpacingVertical)/Args.rows;
xpos1=Args.MarginLeft+cellwidth*(x1-1)+Args.SpacingHorizontal*(x1-1);
xpos2=Args.MarginLeft+cellwidth*x1+Args.SpacingHorizontal*(x1-1);
ypos1=Args.MarginTop+cellheight*(y1-1)+Args.SpacingVertical*(y1-1);
ypos2=Args.MarginTop+cellheight*y1+Args.SpacingVertical*(y1-1);

pos=[xpos1 1-ypos2 xpos2-xpos1 ypos2-ypos1];



