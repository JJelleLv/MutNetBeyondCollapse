function res = i_outfundlg( def)

   if nargin < 1
      def.fun = '';
      def.opt = 1;
   end

   if iscell(def.fun)
       ndx=strfind(def.fun,' ');
       ndx= find(~cellfun('isempty',ndx));
       if ~isempty(ndx)
           for i=1:length(ndx)
               def.fun{ndx(i)}=['''' def.fun{ndx(i)} ''''];
           end
       end
       def.fun=strtrim(sprintf('%s ',def.fun{:}));
   end

   hfig = figure('Units','points', ...
      'Color',[0.914 0.914 0.914], ...
      'MenuBar','none', ...
      'Name','Enter function and source for outfun', ...
      'NumberTitle','off', ...
      'PaperPosition',[18 180 576 432], ...
      'PaperUnits','points', ...
      'Position',[357 384.75 309 142.5], ...
      'Tag','outfun', ...
      'ToolBar','none',...
    'CreateFcn',@(h,evnt)movegui(h, 'center'));

   uicontrol('Parent',hfig, ...
      'Units','points', ...
      'BackgroundColor',[0.914 0.914 0.914], ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[21 117 264 9], ...
      'String','Enter function(s), state/auxiliary variable(s) etc. (space delimited):', ...
      'Style','text', ...
      'Tag','StaticText1');
   uicontrol('Parent',hfig, ...
      'Units','points', ...
      'BackgroundColor',[0.914 0.914 0.914], ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'Position',[21 65.25 122.25 15.75], ...
      'String','Source of data:', ...
      'Style','text', ...
      'Tag','StaticText2');
   hedit = uicontrol('Parent',hfig, ...
      'Units','points', ...
      'BackgroundColor',[1 1 1], ...
      'HorizontalAlignment','left', ...
      'ListboxTop',0, ...
      'String',def.fun, ...
      'Position',[21 99 261 15], ...
      'Style','edit', ...
      'Tooltipstring','Enter equation(s) (including any variable from the model), delimited with space',...
      'Tag','EditText1');
   hopt = uicontrol('Parent',hfig, ...
      'Units','points', ...
      'BackgroundColor',[1 1 1], ...
      'ListboxTop',0, ...
      'Position',[21 52.5 261 15], ...
      'String',{'normal run','paranal','conteq','mcarlo'},... %,'conteq2d'}, ...
      'Style','popupmenu', ...
      'Tooltipstring','The origin of the data for output',...
      'Tag','PopupMenu1', ...
      'Value', def.opt);
   uicontrol('Parent',hfig, ...
      'Units','points', ...
      'Callback',@c_ok,...
      'ListboxTop',0, ...
      'Position',[216 15.75 63 17.25], ...
      'Tooltipstring','Execute outfun',...
      'String','OK', ...
      'Tag','Pushbutton1');
   uicontrol('Parent',hfig, ...
      'Units','points', ...
      'Callback',@c_cancel,...
      'ListboxTop',0, ...
      'Position',[22.5 14.25 63 17.25], ...
      'String','Cancel', ...
      'Tooltipstring','Abort outfun',...
      'Tag','Pushbutton1');
   uiwait(hfig);
   res = [];
   if ishandle(hfig)
      ud = get(hfig, 'userdata');
      if ud
         res.fun = get(hedit, 'string');
         res.opt = get(hopt, 'value');
         res.optstr = get(hopt, 'string');
         res.optstr = res.optstr{res.opt};
      end

      close(hfig);
   end

    function c_cancel(hobject,~)
   set(getparentfig(hobject), 'userdata', 0);
   uiresume;
        function c_ok(hobject,~)
   set(getparentfig(hobject), 'userdata', 1);
   uiresume;


