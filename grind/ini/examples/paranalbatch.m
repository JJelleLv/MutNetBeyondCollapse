% example of a batch script for paranal
%
% open the model
use may76
% make a structure with the dialog results
% (don't need to enter every answer)
%
s.par= {'r'  ''};      %parameter 1 and parameter 2 (advanced options) 
                       %if par{2} is empty only one parameter is varied
s.start=[0 0];         %start value parameter 1 and 2
s.nend=[4 10];         %end value parameter 1 and 2
s.steps=[500 50];      %number of steps parameter 1 and 2
s.stabil= 60;          %time of stabilizing per step
s.writing=30;          %number of time units writing
s.outputtype=1;        %1= unchanged output; 2=Mean 3= Median 4=Maxima
                       %5= Minima; 6= Minima+Maxima 7= SD; 8=CV; 9=Sum; 10=Perc05; 11=Perc95
                       %12=Range90; 13=Min; 14=Max; 15=Range; 16=Returntime'
s.lines= 0;            %0=scatter plot (dots); 
                       %1= line plot; 2= contour (2D only); 3=surface plot(2D only)
paranal(s,'-save=may76.mat'); %perform paranal and write results to may76.mat
era;
use transcrit
y=10;
s.par= 'a';            %other way of selecting only one parameter
s.start=5;
s.nend=0;
s.steps=50;
s.stabil= 60;
s.writing=300;
s.lines=1;
paranal(s,'-save=transcrit.mat');
era;

%load previous paranal
paranal -load=may76.mat
%
% to show the other graph: paranal -lo=transcrit.mat


