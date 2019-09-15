%TIM - Create a time plot
%   Show the data of the last run in a time plot. On the y-axis are the
%   variables or functions that are determined by OUT. It is possible to update
%   several figures simulaneously, see OUT how to define this. If any parameters or
%   initial conditions have been changed, first the model is run for g_grind.ndays.
%
%   Usage:
%   TIM - runs the model for g_grind.ndays time steps.
%   TIM N - (or tim(N)) runs the model for N time steps.
%   TIM -OPT runs the model using the option OPT. Valid options are (short version 
%          between brackets):
%         -add (-a) continue with the previous run.
%         -hold (-h) do not overwrite previous plot.
%         -file=afile (-f)  write the data to a tab delimited file (if the filename is 
%          omitted, a saveas dialog appears.
%         -nocheck (-n)  do not check for changed option, just show last results.
%         -paranal (-p1)  show the results of the last run of PARANAL.
%         -paranal2 (-p2)  show the results of the last two runs of PARANAL
%         -out (-o) select output for current time plot.
%         -out VAR1 FUN2 select VAR1 FUN2 for current time plot.
%         -run (-r) rerun the model always.
%         -silent (-s) do not plot the results.
%
%   TIM  N -OPT1 -OPT2 - combine options and time steps 
%   RES = tim('N','-OPT') - save the results to the matrix RES 
%    
%
%   See also:
%   RU, OUT, SIMTIME, ADDMODE
function [outmat, leg] = tim(varargin)
if nargout>0
   [outmat, leg] = i_time(varargin{:});
else
   i_time(varargin{:});
end


