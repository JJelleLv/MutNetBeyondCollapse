function h = i_waitbar(step, nsteps, caption, title, updatetime, cancelbtn, varargin)
%The cancelbutton sets the userdata of the waitbar to 1
%make sure the program reacts to that
persistent g_waitbar;
global g_grind;
if isfield(g_grind,'figopts')
    f=find(strcmpi(g_grind.figopts,'visible'));
    if ~isempty(f)&&strcmpi(g_grind.figopts{f+1},'off')
        if nargout>0
            h=[];
        end
        return;
    end
end
if isempty(step)
   if ~isempty(g_waitbar)&&~isempty(g_waitbar.h)&&ishandle(g_waitbar.h)
      close(g_waitbar.h);
   end
   g_waitbar = [];
elseif ~isempty(g_waitbar)&&(nargin==1||(nargin==2)&&ischar(nsteps))
   g_waitbar.n = g_waitbar.n+double(step);
   fractiondone=g_waitbar.n/g_waitbar.nsteps;
   atoc = toc;
   if isempty(g_waitbar.h)
      timeleft = (atoc / fractiondone-atoc);
      if timeleft < 1 %too short calculation time to show the bar
         return;
      else
         if ishandle(g_waitbar.h)
            g_waitbar.h = waitbar(fractiondone, g_waitbar.title,'name',g_waitbar.caption,g_waitbar.varargin{:});
         end
         if ~isempty(g_grind.figopts)
             set(g_waitbar.h,g_grind.figopts{:});
         end

         g_waitbar.title = sprintf('%s, time left: %%s%%s',g_waitbar.title);
     end

   elseif ~ishandle(g_waitbar.h)
       g_waitbar=[];
       if nargout>0
           h=[];
       end
       %if the user closes the waitbar no new waitbar is generated.
       return;
   end
   if nargin==2
      g_waitbar.title = sprintf('%s, time left: %%s%%s',nsteps);
   end
   if (atoc > g_waitbar.last + g_waitbar.updatetime)
 %     if ~ishandle(g_waitbar.h)
 %        g_waitbar.h = waitbar(fractiondone);
 %     end

      g_waitbar.last = atoc;
      if ishandle(g_waitbar.h)
         waitbar(fractiondone, g_waitbar.h,gettimeleft(g_waitbar,fractiondone,atoc));
      end
   end
elseif nargin>2&&isempty(g_waitbar)||(step==0)
   tic;
   g_waitbar.nsteps = double(nsteps);
   g_waitbar.n = 0;
   g_waitbar.last = toc;
   if nargin < 3
       caption='';
   end

   if nargin < 4
     title = 'Running';
   end

   if nargin < 5
     updatetime = 1;
   end
   if nargin <6
       cancelbtn=false;
   end
   if nargin < 7
     g_waitbar.varargin = {};
   else
     g_waitbar.varargin=varargin;
   end

   g_waitbar.title =title;
   g_waitbar.updatetime = updatetime;
   g_waitbar.h = [];
   g_waitbar.caption = caption; 
  % if  g_waitbar.nsteps<10
       if cancelbtn
          g_waitbar.h = waitbar(0, g_waitbar.title,'name',g_waitbar.caption,'CreateCancelBtn',@CancelBtnPressed);
          set(g_waitbar.h,'CloseRequestFcn','closereq',g_waitbar.varargin{:});
       else
           g_waitbar.h = waitbar(0, g_waitbar.title,'name',g_waitbar.caption,g_waitbar.varargin{:});
       end
       set(g_waitbar.h,'userdata',0);
       g_waitbar.title = sprintf('%s, time left: %%s%%s',g_waitbar.title);
  % end


end

if nargout > 0
    if isempty(g_waitbar)||~isfield(g_waitbar, 'h')
        h=[];
    else
        h = g_waitbar.h;
    end
end
function s = gettimeleft(g_waitbar, fractiondone, atoc)
if fractiondone>1E-4
    timeleft = (atoc / fractiondone-atoc); %sec
    if timeleft==0
        s=sprintf(g_waitbar.title,0);
    elseif timeleft < 60
        s=sprintf(g_waitbar.title,datestr(timeleft/86400,'SS'),' s');
    elseif timeleft < 3600
        s=sprintf(g_waitbar.title,datestr(timeleft/86400,'MM:SS'),' m');
    elseif timeleft < 86400
        s=sprintf(g_waitbar.title,datestr(timeleft/86400,'HH:MM:SS'),' hr');
    else
        s=sprintf('%dd %s',floor(timeleft/86400),datestr(timeleft/86400,'HH:MM:SS'),'');
    end
else
    s=sprintf(g_waitbar.title,NaN);
end
function CancelBtnPressed(hobject,~)
hfig=getparentfig(hobject);
set(hfig,'Userdata',1);
set(hfig,'Visible','off');
drawnow

