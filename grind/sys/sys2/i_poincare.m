function [t_poincare] = i_poincare(aplane, increasing)
%determines the interpolated times where the plane is crossed
global g_Y g_t g_grind t;
fs=parsed_equation(aplane).fs;
ndx=~strcmp(fs,';');
fs=fs(ndx);
eqs={'==','=>','=<','=','<','>'};
for i=1:length(eqs)
    f=strcmp(fs,eqs{i});
    if any(f)
        f=find(f,1);
        if f==length(fs)-1&&strcmp(fs{end},'0')
            fs=fs(1:end-2);
        else
            fs=[fs(1:f-1) {'-','('} fs(f+1:end) {')'}];
        end

    end

end

aplane=char(minbrackets(parsed_equation(sprintf('%s',fs{:}))));
N0 = i_initvar;
if i_settingschanged(N0)
   i_ru(t, g_grind.ndays, N0, 1);
end

yy = outfun(aplane);
t_poincare = zeros(size(g_t));
ipoin = 0;
if increasing
   for i = 3:size(yy, 1) - 2
      if (yy(i - 1) < 0) && (yy(i) >= 0)
         ipoin = ipoin + 1;
         try
            t_poincare(ipoin) = interp1(yy(i - 2:i + 2), g_t(i - 2:i + 2), 0, 'pchip');
         catch
            t_poincare(ipoin) = g_t(i);
         end

      end

   end

else
   for i = 3:size(g_Y, 1) - 2
      if (yy(i - 1) > 0) && (yy(i) <= 0)
         ipoin = ipoin + 1;
         try
            t_poincare(ipoin) = interp1(yy(i - 2:i + 2), g_t(i - 2:i + 2), 0, 'pchip');
         catch
            t_poincare(ipoin) = g_t(i);
         end

      end

   end

end

t_poincare = t_poincare(1:ipoin, :);

