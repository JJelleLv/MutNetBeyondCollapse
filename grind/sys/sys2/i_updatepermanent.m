%internally used to update (run twice first with at, then without)y
function i_updatepermanent(at)
global g_permanent g_grind;
%update g_permanent
permanent=g_grind.permanent; %update in savemodel
if ~isempty(g_permanent) && ~isempty(g_permanent.t) && g_permanent.active
   if nargin == 1 
      %the problem is that ode45 and others can go back in time tha trials should be ignored then
      if g_permanent.lastt(g_permanent.lasti) > at
    %  fprintf('%g   %g\n',g_permanent.lastt(g_permanent.lasti),at)
         i = g_permanent.lasti - 1;
         g_permanent.lastt(g_permanent.lasti) = -99999;
         if i < 1
            i = 10;
         end

         while ~isnan(g_permanent.lastt(i)) &&  (g_permanent.lastt(i) > at)  && (i~=g_permanent.lasti)
            g_permanent.lastt(i) = -99999;         
            if i == 1
               i = 10;
            else
               i = i - 1;
            end

         end

         g_permanent.lasti =i;
         Y = g_permanent.lastY{g_permanent.lasti};
         if ~isempty(Y)
            for j = 1:length(permanent)
               permanent{j}.currvalue = Y(permanent{j}.from:permanent{j}.to);
               permanent{j}.currvalue = reshape(permanent{j}.currvalue, permanent{j}.dims);
            end

         end

         while (g_permanent.nextt>1) && (g_permanent.t(g_permanent.nextt - 1) >= at)
            g_permanent.nextt = g_permanent.nextt - 1;
         end

      end

      if g_permanent.lasti < 10
         g_permanent.lasti = g_permanent.lasti + 1;
      else
         g_permanent.lasti = 1;
      end
      g_permanent.lastt(g_permanent.lasti)  = at;
   else
      at = g_permanent.lastt(g_permanent.lasti);
      Y = g_permanent.lastY{g_permanent.lasti};
      for i = 1:length(permanent)
         curr=transpose(permanent{i}.currvalue(:));
         lcurr=length(curr);
         ato=permanent{i}.to;
         afrom=permanent{i}.from;
         if lcurr==0
            Y(permanent{i}.from:permanent{i}.to) = NaN;
         elseif (lcurr==1) || (lcurr == ato - afrom + 1)
            Y(afrom:ato) = curr;
         else
            size(curr)
            error('GRIND:defpermantent:NoSizeMatch','Size of permanent variable %s does not match', permanent{i}.name);
         end

      end

      g_permanent.lastY{g_permanent.lasti} = Y;
      it = g_permanent.nextt;
      if it > length(g_permanent.t)
         it = length(g_permanent.t);
      end
      if  (at >= g_permanent.t(it))
         while (it < length(g_permanent.t)) && (at>=g_permanent.t(it+1))
            for i = 1:length(permanent)
               g_permanent.Y(it, :) = Y;
            end

            it = it + 1;
         end

         g_permanent.Y(it, :) = Y;
         g_permanent.nextt = it + 1;
      end

   end

end


