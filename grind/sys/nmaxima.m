%NMAXIMA   Characterize the attractor and number of maxima
%   This function analyses whether an attractor is periodic or not and if it
%   is periodic, how many maxima there are. If the attractor is not periodic
%   it can be chaotic, stochastic or quasi-periodic. It is important to
%   stabilize the model before using nmaxima, see <a href="matlab:help stabil">stabil</a>.
%   It return a structure with a characterisation of the attractor and it
%   can evaluate if two structures are from the same attractor.
%
%   Usage:
%   N=NMAXIMA - creates the structure N with a characterisation of the 
%   attractor. Fields:
%     N.noattractor (0/1) - is 1 if there is no attractor
%     N.cyclic (0/1) - is 1 if there is a cyclic attractor
%     N.nonperiodic (0/1) - is 1 if there are no periodicity detected (chaos?)
%     N.meanperiod - length of a period in time units
%     N.nmaxima - number of maxima per cycle
%     N.indexmaxima - index of the maxima
%     N.maxspecies - the number of the state variable used for maxima
%     N.Ycycle – state variables of one cycle;
%     N.tcycle - times of one cycle;
%   N=NMAXIMA(NSTEPS) - NSTEPS is the number of time units for the run, 
%   NMAXIMA('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'mincorr' [number] - (0.95) min lagged correlation coefficient to be (possibly) periodic
%     'mincv' [number] - (0.02) min coefficient of variation in time lags to be periodic
%     'mindifference' [number>0] - (0.001) for difference equations: the minimum abs difference
%     'mindist' [number>0] - (0.2) minimum Euclidean distance to be cyclic.
%     'nsteps' [number>0] - the number of time units for the run
%     'poverlap' [number] - (0.6) a cycle is different if less than poverlap of the values are overlapping
%     'speciesno' [number] - NaN; species number for maxima (if NaN, the species with the maximum 
%    mean abundance is taken)
%     'stepinterp' [number>0] - (0.5) interpolation step 
%     'struc1' [struct] - structure with the result of a previous run
%     'struc2' [struct] - structure with the result of a previous run
%   NMAXIMA('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-c' STRUC1 STRUC2 - compare two attractors STRUC1 and STRUC2 as analysed with nmaxima.
%     '-d' STRUC1 - display the kind of attractor (STRUC1 is a result of nmaxima)
%     '-n' - suppress rerunning if the parameters were changed
%
%      
%
%   See also stabil, lyapunov, z1test
%
%   Reference page in Help browser:
%      <a href="matlab:commands('nmaxima')">commands nmaxima</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function result = nmaxima(varargin)
%(nsteps, options, par2, par3)
%
%
global g_Y g_t t g_grind;
if isfield(g_grind,'ndays')
    nsteps=g_grind.ndays;
else
    nsteps=1000;
end

fieldnams={'nsteps', 'n>0', 'the number of time units for the run',nsteps;...
   'mindist', 'n>0', '(0.2) minimum Euclidean distance to be cyclic.',0.2;...
   'mindifference', 'n>0', '(0.001) for difference equations: the minimum abs difference',0.001;...
   'mincorr', 'n', '(0.95) min lagged correlation coefficient to be (possibly) periodic',0.95;...
   'poverlap', 'n', '(0.6) a cycle is different if less than poverlap of the values are overlapping',0.6;...
   'mincv', 'n', '(0.02) min coefficient of variation in time lags to be periodic',0.02;...
   'stepinterp', 'n>0', '(0.5) interpolation step',0.5;...
   'speciesno', 'n', 'NaN; species number for maxima (if NaN, the species with the maximum',NaN;...
   'struc1', 'r', 'structure with the result of a previous run',[];...
   'struc2', 'r', 'structure with the result of a previous run',[]}';
args=i_parseargs(fieldnams,'if(hasoption(''-c'')||hasoption(''-d'')),deffields=''struc1,struc2'';else,deffields=''nsteps'';end;','-c,-d,-n',varargin);
 
if any(strcmp(args.opts,'-c'))
   n1 = args.struc1;
   n2 = args.struc2;
   opt.mindiff = 0.05;
   opt.mindistance  = 0.05;
   opt.poverlap = 0.6;
   if nargin > 3
      options = addtostruct(opt, par3);
   else
       options=opt;
   end
   result =  ~(xor(n1.noattractor, n2.noattractor) || xor(n1.cyclic, n2.cyclic) || ...
      xor(n1.nonperiodic, n2.nonperiodic));
   if result
      if ~n1.cyclic && ~n1.noattractor
         result = ((length(n1.Ycycle) == length(n2.Ycycle)) && (sqrt(sum((n1.Ycycle-n2.Ycycle).^2)) < options.mindiff));
      else
      if n1.cyclic && ~n1.nonperiodic
         %  too often it happens that the meanperiod is different
         %  result = (abs(n1.meanperiod - n2.meanperiod) < 1.5) || ...
         %           (abs(n1.meanperiod/2 - n2.meanperiod)<1.5) || ...
         %           (abs(n1.meanperiod - n2.meanperiod/2)<1.5) || ...
         %           (abs(n1.meanperiod/4 - n2.meanperiod)<1.5) || ...
         %           (abs(n1.meanperiod - n2.meanperiod/4)<1.5);
         
         if result
            mindifference = options.mindistance * size(n1.Ycycle, 2);
            result = ((size(n1.Ycycle, 2) == size(n2.Ycycle, 2)) && (sqrt(sum((mean(n1.Ycycle) - mean(n2.Ycycle)).^2)) < mindifference)) ...
               && (sqrt(sum((min(n1.Ycycle) - min(n2.Ycycle)).^2)) < mindifference) && (sqrt(sum((max(n1.Ycycle) - max(n2.Ycycle)).^2)) < mindifference);
         end
      end
      if n1.cyclic && n1.nonperiodic
         %this is a bit rough: it is considered different if of any
         %variable less than options.poverlap of the values are overlapping
         nrow1=size(n1.Ycycle,1);
         nrow2=size(n2.Ycycle,1);
         n=20;
         k1 = floor(linspace(1,nrow1,n));
         k2 = floor(linspace(1,nrow2,n));
         withindist=0;
         betweendist=99999;
         result=0;
         for i=1:n
             dist1=sqrt(sum((n1.Ycycle - repmat(n1.Ycycle(k1(i), :), nrow1, 1)).^2, 2));
             dist2=sqrt(sum((n2.Ycycle - repmat(n2.Ycycle(k2(i), :), nrow2, 1)).^2, 2));
             dist=max(min(dist1(dist1>0)),min(dist2(dist2>0)));
             if dist>withindist
                 withindist=dist;
             end
             dist1=sqrt(sum((n2.Ycycle - repmat(n1.Ycycle(k1(i), :), nrow2, 1)).^2, 2));
             dist2=sqrt(sum((n1.Ycycle - repmat(n2.Ycycle(k2(i), :), nrow1, 1)).^2, 2));
             dist=min(min(dist1),min(dist2));
             if dist<betweendist
                 betweendist=dist;
             end
             if withindist>betweendist
                 result=1;
                 break;
             end
         end
%          p = (1 - options.poverlap) / 2;
%          prc1 =  prctile(n1.Ycycle, [p, 1 - p], 1);
%          prc2 =  prctile(n2.Ycycle, [p, 1 - p], 1);
%          result = (size(n1.Ycycle, 2) == size(n2.Ycycle, 2)) && ~any(prc1(1,:) > prc2(2,:)) && ~any(prc1(2,:) < prc2(1,:)) ...
%             &&  ~any(prc2(1, :) > prc1(2, :)) && ~any(prc2(2, :) < prc1(1, :));
      end
      if n1.noattractor
         result = NaN;
      end
      end
   end
   if nargout == 0
      if isnan(result)
         disp('no attractor');
      elseif result
         disp('same/similar attractor');
      else
         disp('different attractor')
      end
   end
   return;
elseif any(strcmp(args.opts,'-d'))
   res = args.struc1;
   if res.noattractor
      disp('No attractor')
   elseif ~res.cyclic
      disp('No cycles');
   else
      if ~res.nonperiodic
         if res.nmaxima == 1
            fprintf('Limit cycle with period of %g time units\n', res.meanperiod);
         elseif res.nmaxima < 10
            fprintf('Complex cycle, %d maxima, period %g time units\n',res.nmaxima,res.meanperiod);
         end
      else
         if res.nmaxima > 20
            s = '>';
         else
            s = ' ';
         end
         fprintf('Not periodic (may be stochastic, res.nonperiodic or quasiperiodic cycle), %s%d maxima\n',s, res.nmaxima);
      end
   end
   return;
end
i_parcheck;
%default settings:
defopt=struct('nsteps',g_grind.ndays,'mindist',0.2,'mindifference',0.001,'mincorr',0.95,'mincv',0.02,'stepinterp',0.5,'speciesno',NaN);
if ~isnan(g_grind.tstep)
   defopt.stepinterp = defopt.nsteps / g_grind.tstep;
end
options=mergestructs(defopt,args);

oldndays = g_grind.ndays;
if options.nsteps < 1E-30
   error('GRIND:nmaxima:TooFewTimeSteps','Number of time units must be larger than 0');
end
try
   g_grind.ndays = options.nsteps;
   N0 = i_initvar;
   if i_settingschanged(N0)&&~any(strcmp(args.opts,'-n'))
      i_ru(t, g_grind.ndays, N0, 1);
   end
   g_grind.ndays = oldndays;
catch err
   g_grind.ndays = oldndays;
   rethrow(err);
end

if g_grind.solver.isdiffer
   tt = transpose(0:1:options.nsteps);
   Y = interp1(g_t, g_Y, tt);
else
   tt = transpose(0:options.stepinterp:options.nsteps);
   Y = interp1(g_t, g_Y, tt);
end

%prefill to get a logical order
res.noattractor = false;
res.cyclic = false;
res.nonperiodic = false;
res.meanperiod = NaN;
res.nmaxima = 0;
res.indexmaxima = [];
res.maxspecies = 1;
res.Ycycle = [];
res.tcycle = 0;


if isnan(options.speciesno)
   iX = find(max(mean(Y)) == mean(Y));
   if ~isempty(iX)
      res.maxspecies = iX(1);
   else
      res.maxspecies = 1;
   end
else
   res.maxspecies = options.speciesno;
end
%distances = [];
if isnan(sum(Y(end, :))) || isinf(sum(Y(end, :)))
   res.noattractor = true;
   res.Ycycle = g_Y(end, :);
else
   if g_grind.solver.isdiffer
      % find number of unique values
      YY = sort(Y(:, res.maxspecies));
      res.nmaxima = length(YY(diff(YY) > options.mindifference));
      res.cyclic = res.nmaxima > 0;
      res.nonperiodic = res.nmaxima > 32;
      if res.cyclic && ~res.nonperiodic
         res.meanperiod = res.nmaxima + 1;
         Y = Y(end - res.nmaxima:end, :);
         tt = tt(end - res.nmaxima:end);
      else
         res.meanperiod = NaN;
      end
   else
      %euclidean distances between the last value and all other points
      distances = sqrt(sum((Y - repmat(Y(end, :), size(Y, 1), 1)).^2, 2));
      distances = distances(1:end - 1);
      res.meanperiod = NaN;
      res.indexmaxima = [];
      %if all distances are small the run is not cyclic (or not stabilized)
      res.cyclic =  max(distances) > size(Y, 2) * options.mindist;
      res.nonperiodic = 0;
      if res.cyclic
         mindists = find(diff( sign( diff([0; distances; 0]) ) ) > 0);
         res.nonperiodic = 1;
         corrc = zeros(length(mindists), 1);
         for i = 1:length(mindists)
            alag = mindists(i) - mindists(1);
            corrc(i) = min(min(corrcoef(Y(1:end - alag,:), Y(1 + alag:end,:))));
         end
         minperiod = tt(end) + 1;
         mincc = options.mincorr;
         maxsame=0;
         for mincorr = linspace(options.mincorr, 0.999999, 20)
            cc=corrc >= mincorr;
            tims = diff(tt(mindists(cc)));
            if length(tims) > 1
               minper = median(tims);
               res.meanperiod = mean(tims(abs(tims-minper)<2));
               nsame = sum(abs(tims-minper)<2);
               cvperiod = std(tims) / res.meanperiod;
               avcc = mean(corrc(cc));
               if (cvperiod < options.mincv) && (avcc >= mincc) && (nsame >= maxsame)
                  res.nonperiodic = 0;
                  maxsame = nsame;
                  mincc = avcc;
                  minperiod = res.meanperiod;
               end
            end
         end
         if ~res.nonperiodic
            %  times_of_cycles = tt(mindists)-tt(mindists(1));
            %  tim= diff(times_of_cycles(mod(times_of_cycles,firstperiod)<2));
            res.meanperiod = minperiod;
            ndx = tt > (tt(end) - res.meanperiod);
            Y = Y(ndx, :);
            tt = tt(ndx);
         else
            if isempty(mindists)
               res.nonperiodic = 0;
               res.cyclic = 0;
               res.noattractor = 1;
            end
         end
         %to find the maxima of this species:
         res.indexmaxima = find( diff( sign( diff([0; Y(:, res.maxspecies); 0]) ) ) < 0 );
         if (length(res.indexmaxima) > 1) && (res.indexmaxima(1) == 1)
            res.indexmaxima = res.indexmaxima(2:length(res.indexmaxima));
         end
         if (length(res.indexmaxima) > 1) && (res.indexmaxima(end) == length(Y))
            res.indexmaxima = res.indexmaxima(1:length(res.indexmaxima) - 1);
         end
         res.nmaxima = length(res.indexmaxima);
      else
         res.nmaxima  = 0;
      end
   end
end
   if ~res.cyclic && ~res.noattractor
      h = min(100, size(Y,1) - 1);
      res.Ycycle  = mean(Y(end - h:end, :));
      res.tcycle  = 0;
   else
      res.Ycycle = Y;
      res.tcycle = tt - tt(1);
      res.Ycycle(end+1,:)=Y(1,:);
      res.tcycle(end+1)=res.meanperiod;
   end
   result = res;

   
% function r = addtostruct(a, b)
% r = a;
% f = fieldnames(b);
% for i = 1:length(f)
%    r.(f{i}) = b.(f{i});
% end



