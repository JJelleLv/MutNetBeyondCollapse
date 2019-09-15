%PARANAL2D   2D parameter analyser
%   Change two parameters step-by-step and show the attractor by simulation. 
%   Contour plots are created with mean value of the attractor. Can also be done
%   in the advanced options of <a href="matlab:help paranal">paranal</a>.
%
%   Usage:
%   PARANAL2D - user is prompted for information
%   PARANAL2D 1 - replot the previous analysis
%   PARANAL2D -1 - reanalyse the model in the opposite direction
%   PARANAL2D -out (-o) - change the default output in a dialog box.
%   PARANAL2D -out plotno [<param1> <param2> funy2 funz1] [minx maxx] [miny maxy] [minz maxz] 
%   - sets the output in a command line: plotno  = number of plot,  is first parameter
%   funy1 = 1st variable or function for y axis [minx maxx] range for xaxis, etc.
%   PARANAL2D -default (-d) - reset the default output.
%   PARANAL2D -list (-l) - list the outputs for paranal.
%   PARANAL2D -save=filename (-s) - save the results of the last or current paranal to a mat 
%   file with name "filename.mat".
%   PARANAL2D -load=filename (-lo) - load the results of a previous paranal session 
%   (opens if needed the same inifile).
%       PARANAL2D('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'lines' [integer>=0] - 0 = scatter plot, 1= line plot, 2 = contour plot (2D only), 
%       3= surf plot (2D only) 
%     'nend' [number] - the final values of the parameters (2 values)
%     'outputtype' [integer>0 or string] - the kind of output (number or Unchanged Mean Median Maxima   
%       Minima Minima+Maxima SD CV Autoregr Skewness Sum Perc05 Perc95 Range90 Min Max Range)
%     'par' [parameter] - cell with the two parameter for changing
%     'previous' [integer] - replot the previous analysis (1 = replot, -1 = backwards)
%     'stabil' [integer>=0] - number of time units for stabilizing
%     'start' [number] - start value of the parameters (2 values)
%     'steps' [integer>0] - number of steps for changing the parameters
%     'writing' [integer>=0] - number of time units for writing the result (0= one value)
%   PARANAL2D('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-d' - default, reset the default output.
%     '-l' - list (-l) - list the outputs for paranal.
%     '-lo' - -lo=filename load the results of a previous paranal session
%     '-s' - -s=filename save the results of the last or current paranal to a mat 
%      file with name "filename.mat".
%  
%   See also paranal, conteq, forward_stabil
%
%   Reference page in Help browser:
%      <a href="matlab:commands('paranal2d')">commands paranal2d</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function paranal2d(varargin)
%evalin('base','global g_paranal2d;');
%implicit i_parseargs
global g_grind g_paranal;
if isfield(g_grind, 'pars')&&numel(g_grind.pars)<2
   error('GRIND:paranal2d:NoPars','Not enough parameters to analyse');
end
i_parcheck;
if nargin == 0
   if isfield(g_grind, 'paranal')
      answer = g_grind.paranal.dlg;
   end
   answer.paranal2d = 1;
         a1 = i_paranaldialog(answer);
      if isempty(a1)
         return
      elseif ~isempty(answer) && isfield(g_paranal.run, 'Y') && ~isempty(g_paranal.run.Y) && ~isempty(g_paranal.run.parvalues) && ...
            prod(double(a1.start==answer.start)) && strcmp(a1.par{1}, answer.par{1}) ...
            && strcmp(a1.par{2}, answer.par{2}) && prod(double(a1.steps==answer.steps)) ...
            && (a1.stabil==answer.stabil) && (a1.writing==answer.writing) && ...
            prod(double(a1.nend==answer.nend)) && struccmp(g_paranal.settings, par('-v', 0)) && ...
            strcmp(questdlg('Do you want to use data of the previous paranal run?','paranal',...
            'Yes','No','No'),'Yes')
         g_grind.paranal.dlg=a1;
         paranal(1)
         return;
      end
   paranal(a1);   
elseif ischar(varargin{1})&&strncmpi(varargin{1}, '-l', 2)
   if ~isfield(g_grind, 'paranal')
      paranal('-defaults');
   end
   for i = 1:length(g_grind.paranal.plots2d)
      fprintf('paranal2d -out %d [''%s'' ''%s'' ''%s''] %s %s %s\n',i, strtrim(sprintf('%s ',g_grind.paranal.plots2d{i}.xaxis{:})),...
         strtrim(sprintf('%s ',g_grind.paranal.plots2d{i}.yaxis{:})),strtrim(sprintf('%s ',g_grind.paranal.plots2d{i}.zaxis{:})),...
         mat2str(g_grind.paranal.plots2d{i}.xlim), mat2str(g_grind.paranal.plots2d{i}.ylim), mat2str(g_grind.paranal.plots2d{i}.zlim));
   end
elseif ischar(varargin{1})&&strncmpi(varargin{1}, '-o', 2)
   if nargin == 1
      i_paranaloutdlg('paranal2d');
      paranal2d
   else
      p.xaxis = {''};
      p.yaxis = {''};
      p.zaxis = {''};
      p.xlim = [0 10];
      p.ylim = [0 10];
      p.zlim = [0 10];
      ip = 1;
      j = 0;
      for i = 2:nargin
         v = str2num(varargin{i});  %#ok<ST2NM>
         if isempty(v)
            %analyse axes
            s = varargin{i};
            if s(1) == '['
               s = s(2:end - 1);
            end
            f=strfind(s,'''');
            if isempty(f)
               f = strfind(s, ' ');
               ff = repmat(f, 2, 1); %this is a trick to get double values
               f = [0; ff(:); length(s) + 1];
            end
            if length(f) > 1
               p.xaxis = mystr2cell(s(f(1) + 1:f(2) - 1));
            end
            if length(f) > 3
               p.yaxis = mystr2cell(s(f(3) + 1:f(4) - 1));
            end
            if length(f) > 5
               if f(5) < f(6) - 1
                  p.zaxis = mystr2cell(s(f(5) + 1:f(6) - 1));
               end
            end
         elseif length(v) == 1
            ip = v;
         elseif length(v) == 2
            switch j
             case 0
               p.xlim = v;
             case 1
               p.ylim = v;
             case 2
               p.zlim = v;
            end
            j = j + 1;
         end
      end
      if (isempty(p.zaxis)||isempty(p.zaxis{1}))&&(isempty(p.yaxis)||isempty(p.yaxis{1}))
         error('GRIND:paranal:noYaxis','Cannot plot with only the x axis defined');
      end
      g_grind.paranal.plots2d{ip} = p;
      g_grind.paranal.defaultplots = 0;
      g_grind.paranal.currno = 1;
   end
   return;
else
   paranal(varargin{:});
end
function A = mystr2cell(s)
s = outf('changeshortcut', s);
A=str2cell(strrep(s,' ',sprintf('\n'))); %#ok<SPRINTFN>

