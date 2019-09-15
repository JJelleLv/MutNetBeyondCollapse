%SETDEFAULTS   Defaults settings
%   Save/load default settings for plots
%  
%   Usage:
%   SETDEFAULTS saves the default settings
%   SETDEFAULTS 1 reads the default settings
%   SETDEFAULTS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'load' [logical] - if true load the saved default settings else save
%   SETDEFAULTS('-opt1','-opt2',...) - Valid command line options:
%     '-1' - resets the default settings to the system defaults
%
%   Reference page in Help browser:
%      <a href="matlab:commands('setdefaults')">commands setdefaults</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function def = setdefaults(varargin)
%(doload)
global g_grind;
fieldnams={'load', 'l', 'if true load the saved default settings else save',false}';
args=i_parseargs(fieldnams,'load','-1',varargin);
if ~isfield(args,'load')
   args.load = 0;
end
filenam = fullfile(grindpath , 'defaults.mat');
if any(strcmp(args.opts,'-1'))
   delete(filenam);
   setsystemdefaults;
   disp('default settings reset');
elseif args.load
   g_grind.odefile=get_odefile('curr_ode');
   g_grind.xaxis.lim = [0, 10];
   g_grind.yaxis.lim = [0, 10];
   g_grind.zaxis.lim = [0, 10];
   id = fopen(filenam, 'r');
   if id > 0
      fclose(id);
      load(filenam,'defaults');
      if exist('defaults','var')
         g_grind.drawnow = defaults.drawnow;
         g_grind.slowdown = defaults.slowdown;
         g_grind.tstep = defaults.g_tstep;
         g_grind.ndays = defaults.g_ndays;
         if isfield(defaults.g_pen,'nextpen')
             g_grind.pen.cycle = defaults.g_pen.nextpen;
         else
             g_grind.pen.cycle = defaults.g_pen.cycle;
         end
         g_grind.pen.colormap = defaults.g_pen.colormap;
         g_grind.pen.linewidth = defaults.g_pen.linewidth;
         g_grind.pen.fontsize = defaults.g_pen.fontsize;
         g_grind.pen.tickdir = defaults.g_pen.tickdir;
      end
   else
      setsystemdefaults;
   end
else
   defaults.version = 1;
   defaults.drawnow = g_grind.drawnow;
   defaults.slowdown = g_grind.slowdown;
   defaults.g_tstep = g_grind.tstep;
   defaults.g_ndays = g_grind.ndays;
   defaults.g_pen.cycle = g_grind.pen.cycle;
   defaults.g_pen.colormap = g_grind.pen.colormap;
   defaults.g_pen.linewidth = g_grind.pen.linewidth;
   defaults.g_pen.fontsize = g_grind.pen.fontsize;
   defaults.g_pen.tickdir = g_grind.pen.tickdir;
   disp(defaults)
   save(filenam, 'defaults');
end
if nargout == 1
   def = defaults;
end

function setsystemdefaults
global g_grind;
g_grind.slowdown = 0;
g_grind.drawnow = 1;
g_grind.tstep = NaN;
g_grind.ndays = 1000;
g_grind.pen = nextpen([]);
return

function odef = get_odefile(odefile)
j = length(odefile);
while (j > 0) && strcontains('1234567890', odefile(j))
   j = j - 1;
end
odefile1 = [odefile(1:j) '%d.m'];
i = 0;
while exist(fullfile(grindpath , sprintf(odefile1, i)), 'file') 
   i = i + 1;
end
odef = sprintf([odefile(1:j) '%d'], i);

