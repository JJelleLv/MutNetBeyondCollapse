function i_ru(at, ndays, N0, updatesettings)
global g_t g_Y g_grind t g_data;
if ~isempty(g_grind.onstart.funs)
   onstart;
end
if g_grind.solver.isstochastic
    dwiener('-i');
    djump('-i');
end
if numel(ndays)>1
    ts=ndays;
    ndays=ts(end)-ts(1);
elseif ~isnan(g_grind.tstep)
    ts = transpose(at:ndays / g_grind.tstep:at + ndays);
    if length(ts)==2
        ts=[ts(1);ts(1)+(ts(2)-ts(1))/2;ts(2)];
    end  
else
    ts = [at; at + ndays];
end

try
   % includepars('-u')
    rednoise('-u');
    implicitdisperse;
catch
    %these are not critical errors
end

if ~isempty(g_grind.permanent)
    defpermanent('-initiate',[at,ndays]);
end

if updatesettings
 %   g_func = [];
    g_grind.checks.lastsettings = i_getsettings(N0,ndays);
end

if g_grind.solver.reset_at_data
   if isfield(g_data, 'obs')
      tt = g_data.t;
      if at < tt(1)
         tt = [at; tt];
      end

      g_tt = at;
      g_YY = transpose(N0);
      oldtstep = g_grind.tstep;
      try
         g_grind.solver.reset_at_data = 0;
         g_grind.tstep = 3;
         for i = 1:length(tt) - 1
            N0_obs = transpose(g_data.obs(i, :)); %reset to the observation
            N0_obs(isnan(N0_obs)) = N0(isnan(N0_obs));
            i_ru(tt(i),tt(i + 1) - tt(i), N0_obs,0);
            g_tt = [g_tt; g_t(1, :) + 0.01; g_t(end, :)];
            g_YY = [g_YY; g_Y(1, :); g_Y(end, :)];
            N0 = transpose(g_Y(end, :));
         end

         if ndays - at > tt(end)
            g_grind.tstep = oldtstep - length(g_tt);
            i_ru(tt(end),ndays - tt(end), N0_obs,0);
            g_tt = [g_tt; g_t];
            g_YY = [g_YY; g_Y];
         end

         g_Y = g_YY;
         g_t = g_tt;
         g_grind.tstep = oldtstep;
         g_grind.solver.reset_at_data = 1;
      catch
         g_grind.solver.reset_at_data = 1;
         g_grind.tstep = oldtstep;
      end

      return;
   end

end
if g_grind.solver.addmode
   oldg_t = g_t;
   oldg_Y = g_Y;
   if g_grind.solver.haslags&&isfield(g_grind.solver,'newhist')&&~isempty(g_grind.solver.newhist)
       if isstruct(g_grind.solver.history)&&isstruct(g_grind.solver.newhist)
           g_grind.solver.history.x=[g_grind.solver.history.x g_grind.solver.newhist.x(2:end)];
           g_grind.solver.history.y=[g_grind.solver.history.y g_grind.solver.newhist.y(:,2:end)];
           g_grind.solver.history.yp=[g_grind.solver.history.yp g_grind.solver.newhist.yp(:,2:end)];
       else
          g_grind.solver.history=g_grind.solver.newhist;
       end

        g_grind.solver.newhist=[];
   end

   if ~isempty(g_t)
      if g_grind.solver.backwards
         ts = ts-t-g_t(1);
         N0 = transpose(g_Y(1, :));
         i_keep(N0);
      else
         ts = ts-t+g_t(end);
         N0 = transpose(g_Y(size(g_Y, 1), :));
         i_keep(N0)
      end
      at=ts(1);
   end

end

if g_grind.solver.backwards
      if g_grind.solver.haslags
         g_grind.solver.name = 'ode23';
      end
end
 odefile=i_getodehandle(0,'');
try
   if isfield(g_grind, 'boxcar')
       g_grind.boxcar.gcycl=zeros(size(g_grind.boxcar.gcycl));
   end

   if g_grind.solver.hasevents
      settings = i_getsettings(N0);
      N00=N0;
      [~,NP0] = i_initvar;
     try
         nextevent=min(ndays,i_setevent('runnext','',-1));
         t0 = at;
         te = at + ndays;
         g_t = [];
         g_Y = [];
         oldt=t;
         while t0 < te
            if t0 <= nextevent
               [N0,NP] = i_initvar;
               if ~isnan(g_grind.tstep)
                  %create equally spaced output
                  ts1 = transpose((t0:ndays /g_grind.tstep:nextevent));
                  if ts1(end) < nextevent
                     ts1 = [ts1; nextevent];
                  end

                  if length(ts1) < 5 %at least 5 steps
                     ts1 = transpose((t0:(nextevent - t0) / 5:nextevent));
                  end

               else
                  ts1 = [t0; nextevent];
               end

               if (length(ts1) > 1) && (ts1(1) < ts1(end))
                  if ~isempty(g_grind.permanent) 
                     defpermanent('-activate',NP);
                  end

                  i_keep(N00) %needed if you want to access the initial value from the program (val)
                  g_grind.hfun.curr=i_getodehandle(0,'');
                  [tt, YY] = feval(str2func(g_grind.solver.name), g_grind.hfun.curr, ts1, N0, g_grind.solver.opt);
                  g_t = [g_t; tt];
                  g_Y = [g_Y; YY];
                  ke; %now the event can change the values of the state variables
                  t=g_t(end);
               end

               t0 = nextevent + 1E-10;
            end

            nextevent=i_setevent('runnext','',nextevent);
            if ~isempty(g_grind.permanent)
               defpermanent('-updatevars');
            end
            if isnan(nextevent) || (nextevent > te)
               nextevent = te;
            end

         end

         if numel(ts)>2
      %      ts = transpose((at:ndays / g_grind.tstep:te));
            g_Y(end,:) = transpose(i_initvar);
            g_Y = interp1(g_t, g_Y, ts);
            g_t = ts;
         end

         i_setsettings(settings);
         if ~isempty(g_grind.permanent)
           defpermanent('-s',NP0);
         end

         t=oldt;
     catch err
%          err=lasterror;
          i_setsettings(settings);
          if ~isempty(g_grind.permanent)
            defpermanent('-s',NP0);
          end

          t=oldt; 
         rethrow(err);
      end

   else


%       if g_grind.version.isoctave
%          g_grind.solver.opt.StepSize = 0.2; %necessary??
%          g_grind.solver.opt.OutputFcn = [];
%       end

      %#function ode45, euler, i_differ
      [g_t, g_Y] = feval(str2func(g_grind.solver.name),odefile, ts, N0, g_grind.solver.opt);
      if isempty(g_t)
          g_t=NaN;
          g_Y(end+1,:)=NaN(1,size(g_Y,2));          
      end
      if g_t(end)<ts(end)
          g_t(end+1)=ts(end);
          g_Y(end+1,:)=NaN(1,size(g_Y,2));
      end
      if size(g_t,2)>1
         g_t=g_t(:);
      end

   end

catch err
%    err =lasterror;
 %    g_grind.lastsettings.initvar = [];
    rethrow(err);
end

if g_grind.solver.backwards
   if g_grind.solver.haslags
       if g_grind.dde.isvariable
          g_grind.solver.name = 'ddesolsd';
       else
          g_grind.solver.name = 'ddesol';
       end
   end

   g_t = -flipud(g_t);
   g_Y = flipud(g_Y);
end

if g_grind.solver.addmode
   if g_grind.solver.backwards
      g_t = [g_t; oldg_t(2:end)];
      g_Y = [g_Y; oldg_Y(2:end,:)];
   else
      g_t = [oldg_t(1:end-1); g_t];
      g_Y = [oldg_Y(1:end-1,:); g_Y];
   end

end




