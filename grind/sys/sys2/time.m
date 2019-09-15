%TIME   Create a time series plot
%   Show the data of the last run in a time plot. On the y-axis are the
%   variables or functions that are determined by OUT. It is possible to update
%   several figures simultaneously, see <a href="matlab:help out">out</a> how to define this. If any parameters or
%   initial conditions have been changed, first the model is run for g_grind.ndays.
%
%   Usage:
%   TIME - runs the model for g_grind.ndays time units.
%   TIME NDAYS - (or time(N)) runs the model for NDAYS time units.
%   RES = time('N','-OPT') - save the results to the matrix RES
%   TIME('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'ndays' [number>0] - runs the model for NDAYS time units.
%     'timevars' [string] - list of equations, use in combination with '-out'
%   TIME('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-1 or -2 or -value' - use in combination with -o to define the number of the time plot (see <a href="matlab:help out">out</a>)
%     '-a' - add: continue with the previous run.
%     '-f=afile' - write the data to a tab delimited file (if the filename is omitted, a saveas dialog appears.
%     '-h' - hold: do not overwrite previous plot.
%     '-n' - nocheck: do not check for changed option, just show last results.
%     '-o' - select output for current time plot.
%     '-p' - makes it possible to replay a paranal
%     '-p1' or '-pa' - show the results of the last run of PARANAL.
%     '-p2' or '-paranal2' - show the results of the last two runs of PARANAL
%     '-r' - rerun the model always.
%     '-s' - silent (-s) do not plot the results.
%   
%
%   See also ru, out, simtime, addmode, outfun, replayall
%
%   Reference page in Help browser:
%      <a href="matlab:commands('time')">commands time</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function [outmat, leg] = time(varargin)
if nargout>0
   [outmat, leg] = i_time(varargin{:});
else
   i_time(varargin{:});
end
