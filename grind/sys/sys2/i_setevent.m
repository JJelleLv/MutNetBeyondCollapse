function nextt = i_setevent(flag, event, dt, varargin)
global g_grind;
switch flag
%  case 'update'
%    g_grind.event.queue = g_grind.event.events;
%    i_setevent sortqueue;
%    %   g_grind.event.events=g_grind.event.queue;
%    i_setevent filldlg;
%  case 'EditButton'
%    v=get(findobj(gcf,'Tag','EventsList'),'value');
%    ev=i_setevent('check',get(findobj(gcf,'Tag','EventEdit'),'string'),...
%       get(findobj(gcf,'Tag','TimeEdit'),'string'),...
%       evalin('base',get(findobj(gcf,'Tag','ArgEdit'),'string')));
%    g_grind.event.events(v) = ev;
%    i_setevent update;
%  case 'AddButton'
%    i_setevent('addnew',get(findobj(gcf,'Tag','EventEdit'),'string'),...
%       get(findobj(gcf,'Tag','TimeEdit'),'string'),...
%       evalin('base',get(findobj(gcf,'Tag','ArgEdit'),'string')));
%  case 'DeleteButton'
%    f = gcf;
%    if strcmp(questdlg('Are you sure to delete the current event?','Setevent','Yes','No','No'),'Yes')
%       v=get(findobj(f,'Tag','EventsList'),'value');
%       g_grind.event.events = g_grind.event.events([1:v - 1 v + 1:end]);
%       i_setevent update;
%    end

%  case 'ListButtonDwn'
%    v=get(findobj(gcf,'Tag','EventsList'),'value');
%    if isfield(g_grind, 'event')
%       if ~isempty(g_grind.event.queue)
%          ev = g_grind.event.queue(v);
%          set(findobj(gcf,'Tag','EventEdit'),'string',ev.event);
%          set(findobj(gcf,'Tag','TimeEdit'),'string',ev.strt);
%          set(findobj(gcf,'Tag','ArgEdit'),'string',i_cell2str(ev.arg));
%       end

%    end

%  case 'filldlg'
% %   set(findobj(gcf,'Tag','EventsList'),'string',s);
%    if isfield(g_grind, 'event')
%       s = cell(1,length(g_grind.event.queue));
%       for i = 1:length(g_grind.event.queue)
%          ev = g_grind.event.queue(i);
%          if length(ev.arg) == 1
%             s{i}=sprintf('Event %d: %s(%s) t=%s', i, ev.event, tostr(ev.arg{1}), ev.strt);
%          elseif length(ev.arg) > 1
%             s{i}=sprintf('Event %d: %s(%s,%s) t=%s',i,ev.event,tostr(ev.arg{1}),tostr(ev.arg{2}),ev.strt);
%          else
%             s{i}=sprintf('Event %d: %s t=%s', i, ev.event, ev.strt);
%          end

%       end

%    else
%        s={};
%    end

%    set(findobj(gcf,'Tag','EventsList'),'string',s);
%  case 'ClearallButton'
%    if strcmp(questdlg('Are you sure to delete all events?','Setevent','Yes','No','No'),'Yes')
%       i_setevent clear;
%       i_setevent filldlg;
%    end

 case 'check'
   if ~isempty(event)
      if isempty(which(event))
         s = sprintf('Event handler for "%s" not found, must be a function',event);
         i_errordlg(s);
         error('GRIND:setevent:FuncNotFound',s);  %#ok<SPERR>
      end

      ev.event = event;
      if ~isfield(ev, 'enabled')
         ev.enabled = 1;
      end

      if ischar(dt)
         try
            ev.t = i_checkstr(dt);
         catch 
            ev.t = 0;
         end

          ev.strt = dt;
      else
         ev.t = dt;
         ev.strt = num2str(dt);
      end

      if nargin > 3
         ev.arg = varargin{:};
      else
         ev.arg = [];
      end

   else
      ev = [];
   end

   nextt = ev;
 case 'addnew'
   button = [];
   for i = length(varargin{1}):-1:1
      if strcmpi(varargin{1}{i}, '-replace')
         button = 'Replace';
         varargin{1} = {varargin{1}{1:i - 1} varargin{1}{i + 1:end}};
      elseif strcmpi(varargin{1}{i}, '-add')
         button = 'Add';
         varargin{1} = {varargin{1}{1:i - 1} varargin{1}{i + 1:end}};
      end

   end

   ev =  i_setevent('check', event, dt, varargin{:});
   if ~isempty(ev)
      if isfield(g_grind, 'event') && ~isempty(g_grind.event.events)
         btn = 'Add';
         for i = 1:length(g_grind.event.events)
            if g_grind.event.events(i).t == ev.t
               if ~isempty(button)
                  btn = button;
               else
                  btn=questdlg('There is already an event at that time','setevent','Add','Replace','Cancel','Replace');
               end

               if strcmp(btn, 'Replace')
                  ii = i;
                  break
               end

            end

         end

         if strcmp(btn, 'Replace')
            g_grind.event.events(ii) = ev;
         elseif strcmp(btn, 'Add')
            g_grind.event.events(end + 1) = ev;
         end

      else
         g_grind.event.events = ev;
      end

   end

   g_grind.event.t = 0;
   if ~isfield(g_grind.event, 'events')
      g_grind.event.events = [];
   end

   g_grind.event.queue = g_grind.event.events;
   i_setevent('sortqueue');
   % g_grind.event.events = g_grind.event.queue;
   g_grind.solver.hasevents = ~isempty(g_grind.event.queue);
   g_grind.checks.lastsettings = [];
 case 'list'
   %disp(sprintf('Last time: %d',g_grind.event.t));
   for i = 1:length(g_grind.event.events)
      ev = g_grind.event.events(i);
      args='';
      if length(ev.arg) >= 1
         args='(';
         for j=1:length(ev.arg)-1
            args=sprintf('%s%s,',args,tostr(ev.arg{j}));
         end

         args=sprintf('%s%s)',args,tostr(ev.arg{end}));
      end

      if ev.enabled
         fprintf('Event %d: %s%s t=%s\n',i,ev.event,args,ev.strt);
      else
         fprintf('Disabled event %d: %s%s t=%s\n',i,ev.event,args,ev.strt);
      end
   end

 case 'runnext'
   if nargin == 1
      if ~isempty(g_grind.event.queue)
         dt = g_grind.event.queue(1).t;
      else
         nextt = NaN;
         return;
      end

   end

   if dt < g_grind.event.t
      g_grind.event.queue = g_grind.event.events;
      for i = 1:length(g_grind.event.queue)
         g_grind.event.queue(i).t = i_checkstr(g_grind.event.queue(i).strt);
      end

      i_setevent('sortqueue');
      g_grind.solver.hasevents = ~isempty(g_grind.event.queue);
   end

   g_grind.event.t = dt;
   i = 1;
   while (i <= length(g_grind.event.queue)) && (g_grind.event.queue(i).t <= g_grind.event.t+0.0001)
      ev = g_grind.event.queue(i);
      if ~isempty(ev.arg)
         next = feval(ev.event, ev.t, ev.arg{:});
      else
         next = feval(ev.event, ev.t);
      end

      if next < ev.t
         next = NaN;
      end

      g_grind.event.queue(i).t = next;
      i = i + 1;
   end

   i_setevent('sortqueue');
   if ~isempty(g_grind.event.queue)
      nextt = g_grind.event.queue(1).t;
   else
      nextt = NaN;
   end

 case 'sortqueue'
   if ~isempty(g_grind.event.queue)
       g_grind.event.queue=g_grind.event.queue([g_grind.event.queue(:).enabled]==1);
   end

   if ~isempty(g_grind.event.queue)
%       for i = 1:size(g_grind.event.queue)
%          if ~g_grind.event.queue(i).enabled
%             g_grind.event.queue(i).t = NaN;
%          end

%       end

      [~, ix] = sort([g_grind.event.queue(:).t]);
      e = length(g_grind.event.queue);
      g_grind.event.queue = g_grind.event.queue(ix);
      while (e > 0) &&  isnan(g_grind.event.queue(e).t)
         e = e-1;
      end

      g_grind.event.queue = g_grind.event.queue(1:e);
   end

 case 'enable'
   for i = 1:length(g_grind.event.events)
      g_grind.event.events(i).enabled = event;
   end

   g_grind.event.queue = g_grind.event.events;
   i_setevent sortqueue;
   g_grind.solver.hasevents = ~isempty(g_grind.event.queue);
 case 'clear'
   g_grind.event.queue = [];
   g_grind.event.events = [];
   g_grind.event.t = 0;
   g_grind.solver.hasevents = 0;
   g_grind.checks.lastsettings = [];
end

function s = tostr(a)
if ischar(a)
   s = a;
else
   s = num2str(a(1));
end

