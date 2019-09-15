%COMBFIG   Combine several figures
%   Combine several existing figures in a NxN matrix of figures. Note: 
%   <a href="matlab:help era">era</a> (erase) does not clear this figure.
%
%   Usage:
%   COMBFIG - Combine all currently opened figures in one figure.
%   COMBFIG [FIGNO1 FIGNO2.. FIGNOn] - combine the listed figure numbers.
%   COMBFIG [FIGNO1 FIGNO2.. FIGNOn] [ROW COL] - combine the listed 
%   figures, but override the default number of rows/columns (max 4).
%   COMBFIG [FIGNO1 FIGNO2.. FIGNOn] [ROW COL] START - Define a 
%   starting position (default 1). This is useful for combining several
%   versions of the same figure, for instance to combine 4 versions of
%   Figure 1: COMBFIG 1 [2 2] 1; COMBFIG 1 [2 2] 2;
%   COMBFIG 1 [2 2] 3; COMBFIG 1 [2 2] 4 [..];
%  
%   COMBFIG('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'axis' [x | y | b] - name of axis (x=x axis, y=yaxis, b=both), see '-r'.
%     'cells' [integer>0 and length(integer)<=2] - [ROW COL] number of rows and columns.
%     'figs' [handle] - list of figure handles or numbers
%     'figno' [integer>0] - number of the figure that is created.
%     'spacing' [number] - spacing between panels (see '-p')
%     'start' [integer>0] - starting position of the first figure.
%   COMBFIG('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-p' - repositions the panels made with combfig or subplot such that 
%         there is less whitespace.
%     '-p' SPACING  - keep SPACING between the axes in relative units,
%        SPACING = 0.02 is default, SPACING=0 is axis without spacing.
%     '-p' [0 0.05]  - spacing is 0 vertically and 0.05 horizontally.
%     '-r' - Remove overlapping axes and tick texts (assuming that all 
%        axes are the same)
%     '-r' x - Remove only x-axis (y=y-axis; b = both (default)).
%     '-l' - adjust tick labels for zero space ('-p' 0)
%      
%See also copyfig, era, <a href="matlab:help subplot">subplot</a>, replayall
%
%   Reference page in Help browser:
%      <a href="matlab:commands('combfig')">commands combfig</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
function ndx = combfig(varargin)
%(fignrs, Cells, starti)
if ~exist('i_figno','file')
    addpath([grindpath filesep 'sys2'])
end
fieldnams={'figs','h','list of figure handles or numbers',double(get(0,'children'));...
   'cells','i>0&length(i)<=2','[ROW COL] number of rows and columns.',[];...
   'figno','i>0','number of the figure that is created.',i_figno('combfig');...
   'start','i>0','starting position of the first figure.',1;...
   'axis','e[x|y|b]','name of axis (x=x axis, y=yaxis, b=both), see ''-r''.','b';...
   'spacing','n','spacing between panels (see ''-p'')',[]}';
args=i_parseargs(fieldnams,'if hasoption(''-p''),deffields=''spacing'';elseif hasoption(''-r'')||hasoption(''-l''),deffields=''axis'';else,deffields=''figs,cells,start'';end;','-p,-r,-l',varargin);
if nargout == 1&&isfield(args,'figs')&&isfield(args,'cell')&&isfield(args,'start')
    ndx = sub2ind(flipud(args.figs(:)), args.start, args.cells); %irritating that we need flipud as the cells are different defined
    return;
end
if ~isfield(args,'figs')
    args.figs=get(0, 'children');
end
if any(strcmp(args.opts, '-l')) %adjust tick labels for zero space
    hs = getsubplotgrid(gcf);
    if ~isfield(args,'axis')
        orientation = 'b';
    else
        orientation = args.axis;
    end
    [rows, cols] = size(hs);
    if strncmpi(orientation,'b',1) || strncmpi(orientation,'x',1)
        for i = 1:rows - 1
            for j = 1:cols
                if ishandle(hs(i, j))&&(hs(i, j)~=0)
                    adjustticklabels(hs(i, j), 'X')
                end
            end
        end
    end
    if strncmpi(orientation,'b',1) || strncmpi(orientation,'y',1)
        for i = 1:rows
            for j = 1:cols
                if ishandle(hs(i, j))&&(hs(i, j)~=0)
                    adjustticklabels(hs(i, j), 'Y')
                end
            end
        end
    end
    return;
end
if  any(strcmp(args.opts, '-r')) %remove axis
    hs = getsubplotgrid(gcf);
    if numel(hs)<2
        disp('Too few axes');
        return;
    end
    if ~isfield(args,'axis')
        orientation = 'b'; %both
    else
        orientation = args.axis;
    end
    [rows, cols] = size(hs);
    if strncmpi(orientation,'b',1) || strncmpi(orientation,'x',1)
        for i = 1:rows - 1
            for j = 1:cols
                if ishandle(hs(i, j))&&(hs(i, j)~=0)
                    set(get(hs(i,j),'xlabel'),'string','')
                    set(hs(i, j), 'XTickLabel', []);
                end
            end
        end
        %       for i = 2:rows
        %          for j = 1:cols
        %             if ishandle(hs(i, j))&&(hs(i, j)~=0)
        %                set(get(hs(i,j),'title'),'string','')
        %             end
        %          end
        %       end
    end
    if strncmpi(orientation,'b',1) || strncmpi(orientation,'y',1)
        for i = 1:rows
            for j = 2:cols
                if ishandle(hs(i, j))&&(hs(i, j)~=0)
                    set(get(hs(i,j),'ylabel'),'string','')
                    set(hs(i, j), 'YTickLabel', []);
                end
            end
        end
    end
    set(get(gcf,'children'),'fontsize',12)
    aa=get(get(gcf,'children'),'xlabel');
    for i = 1:length(aa)
        set(aa{i}, 'fontsize', 12);
    end
    aa=get(get(gcf,'children'),'ylabel');
    for i = 1:length(aa)
        set(aa{i}, 'fontsize', 12);
    end
    return;
end
if any(strcmp(args.opts, '-p')) %adjust positions
    set(gcf,'PaperPositionMode','auto')
    set(gcf,'units','normalized');
    figpos = get(gcf, 'position');
    vert = 0.5469; %normal size of figure
    hor = 0.41;
    hs = getsubplotgrid(gcf);
    if isempty(hs)
        disp('No axes to adjust');
        return;
    elseif length(hs)==1
        return;
    end
    if ~isfield(args,'spacing')
        spacing = 0.02;
    else
        spacing = args.spacing;
    end
    if length(spacing) == 1
        spacing = spacing + zeros(2, 1);
    end
    [rows, cols] = size(hs);
    if cols<=4&&rows<=4
        newfigpos =  [figpos(1:2) hor * cols / 2 vert * rows / 2];
    else
        rat = rows / cols;
        if rat < 1
            newfigpos =  [figpos(1:2) hor vert * rat];
        else
            newfigpos =  [figpos(1:2) hor / rat vert];
        end
    end
    if newfigpos(4)>1
        newfigpos(2)=-0.2;
    end
    set(gcf, 'position', newfigpos);
    for i = 1:rows
        for j = 1:cols
            h=findobj(gcf,'tag',sprintf('subplot(%d,%d)',i,j));
            if any(ishandle(h))&& any(h~=0)
                set(h, 'position', subplotpos(rows, cols,j, i, spacing));
            end
        end
    end
    return;
end
if ~isfield(args,'start')
    args.start = 1;
end
nfig = length(args.figs);
if ~isfield(args,'cells')
    if nfig  == 0
        i_errordlg('Not enough figures to combine');
        error('GRIND:combfig:TooFewFigs','Not enough figures to combine');
    elseif nfig == 1
        args.cells = [1 1];
    elseif nfig == 2
        args.cells = [2 1];
    elseif nfig == 3
        args.cells = [3 1];
    elseif nfig == 4
        args.cells = [2 2];
    elseif nfig <= 6
        args.cells = [3 2];
    elseif nfig <= 8
        args.cells = [4 2];
    else
        args.cells = [3 3];
    end
end
scale = zeros(1, 2);
for i = 1:2
    if i == 1
        j = 2;
    else
        j = 1;
    end
    if args.cells(i) == 1
        scale(j) = 1;
    elseif args.cells(i) == 2
        scale(j) = 0.45;
    elseif args.cells(i) == 3
        scale(j) = 0.27;
    else
        scale(j) = 0.2;
    end
end
if isfield(args,'figno')&&~isempty(args.figno)
    h = figure(args.figno);
elseif exist('i_figno','file')
    h = i_figure(i_figno('combfig'));
    cm=get(args.figs(1),'colormap');
    set(h, 'colormap', cm);
else
    h = figure;
end
set(h, 'Name', 'Combined figure');
ud.rows = args.cells(1);
ud.cols = args.cells(2);
set(h, 'userdata', ud);

for i = 1:nfig
    %   s = subplot(args.cells(1), args.cells(2), i + args.start - 1);
    [i1, i2] = ind2sub(args.cells, i + args.start - 1);
    s=subplot(ud.rows,ud.cols,i + args.start - 1,'tag',sprintf('subplot(%d,%d)',i1,i2));
    %  s=subplot('position',subplotpos(args.cells(1),args.cells(2),i1,i2),'tag',sprintf('subplot(%d,%d)',i1,i2));
    %hs = zeros(ud.rows, ud.cols) - 1;
    
    pos = get(s, 'Position');
    tag = get(s, 'tag');
    delete(s);
    if ~ishandle(args.figs(i))
        error('GRIND:combfig:cannotfindfig','Error combfig: cannot find figure %g', args.figs(i));
    else
        hax=findall(args.figs(i),'type','axes');
        if isempty(hax)
            error('grind:combfig:noaxes','Cannot combine figures without axes')
        end
        %    arrs=findall(args.figs(i),'tag','Arrow');
    % arrows are not added correctly. Too complex to repair it seems.
        tag1 = get(hax, 'tag');

        hax=hax(~ (strcmp(tag1,'legend')|strcmp(tag1,'Colorbar')));
        %pos2 = get(ax, 'position');
        %       l=findobj(ax,'tag','legend');
        %       if ~isempty(l)
        
        hs = copyobj(hax, h);
        set(hs, 'Position', pos);
        set(hs, 'tag', tag);
%         if ~isempty(arrs)
%           hfig=get(hs(2),'parent');
%           arrowsX=get(arrs,'X');
%             arrowsY=get(arrs,'Y');
%             newarrows=findall(hfig,'tag','Arrow');
%             delete(newarrows);
%             for k=1:length(arrowsX)
%                 pnt=[mean(arrowsX{k}),mean(arrowsY{k})];
%                 arrows(hs(1),'nearest_normalized',pnt);
%             end
%         end

        %        for j = 1:length(hs)
        %           pos2 = get(hs(j), 'Position');
        %           pos2(1:2) = pos(1:2) + scale(1:2) .* (pos2(1:2) - 0.13);
        %           pos2(3:4) = scale(1:2) .* pos2(3:4);
        %           set(hs(j), 'Position', pos2);
        %        end
    end
end
i_grindlegend(13);
if isfield(args,'spacing')
    combfig('-p',args.spacing);
end
function hs = getsubplotgrid(h)
ch = get(h, 'children');
hs=[];
if ~isempty(ch)
    tags = get(ch, 'tag');
    ch = ch(~strcmp(tags, 'legend'));
    types = get(ch, 'type');
    if ~iscell(types)
        types={types};
    end
    poss = get(ch, 'position');
    if ~iscell(poss)
        poss={poss};
    end
    ipos = zeros(length(ch), 4);
    for i = length(ch):-1:1
        if  strcmp(types{i}, 'axes')
            ipos(i, :) = poss{i};
        else
            ipos(i, :) = [];
        end
    end
    colpos = sort(unique(sort(ipos(:, 1))), 'ascend');
    rowpos = sort(unique(sort(ipos(:, 2))), 'descend');
    hs = zeros(length(rowpos), length(colpos));
    for i = 1:length(ch)
        if strcmp(types{i}, 'axes')
            arow=find(rowpos == ipos(i, 2));
            acol=find(colpos == ipos(i, 1));
            hs(arow, acol) = ch(i);
            set(ch(i),'tag',sprintf('subplot(%d,%d)',arow,acol));
        end
    end
end
function adjustticklabels(hax, orient)
axis(hax, 'tight')
newtick = get(hax, [orient 'tick']);
tickdiff = (newtick(2) - newtick(1));
newlim = [newtick(1) - tickdiff newtick(end) + tickdiff];
axis(hax, 'manual')
set(hax, [orient 'lim'], newlim);
set(hax, [orient 'tick'], [newtick(1) - tickdiff newtick]);

