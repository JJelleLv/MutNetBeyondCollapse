function s = i_disptext(s)
if strcontains(s,'{\') %no double tex
    return;
end

if strcontains(s, 'timesens(')
    s=strrep(s,' ,',',');
    s=strrep(s,'''-p'',','');
    obj=parsed_equation(s);
    p1=obj.fs;
    %   p1=parseeq(s);
    s1=find(strcmp(p1,'timesens'));
    popbrack=cumsum(strcmp(p1,'(')-strcmp(p1,')'));
    for j=length(s1):-1:1
        pops=find(popbrack==popbrack(s1(j)));
        f=find(pops==s1(j));
        f1=pops(f(1));
        f2=pops(f(1)+1);
        eq=p1(f1:f2);
        s2=eq{find(strncmp(eq,'''',1),1)};
        %s2 is the text of timesens
        parvar=regexp(s2,'[#A-Za-z0-9_\(\),]*(?![\[])','match');
        if length(parvar)==2 %normal d[par]/d[var]
            var = parvar{1};
            par1 = strrep(parvar{2},'(0)','');
            par2=strrep(par1,'(0)','_0');
            if var(1)=='#'
                g_sens=evalin('base','g_sens');
                ipar= strcmp(g_sens.pars,par1);
                ivar=str2double(var(2:end));
                var=i_statevars_names(g_sens.maxndx(ivar,ipar));
            end

            latextext=sprintf('\\partial%s/\\partial%s',var,par2);
            if f1>1||(f2<length(p1)&&isempty(strfind('*+-/',p1{f2+1})))
                p1=[p1(1:f1-1) {'(',latextext,')'} p1(f2+1:end)];
            else
                p1=[p1(1:f1-1) {latextext} p1(f2+1:end)];
            end

        else  %statevar with a rank: statevar([#varno@par])
            var = parvar{2};
            par1 = parvar{3};
            par1=strrep(par1,'(0)','');
            if var(1)=='#'
                g_sens=evalin('base','g_sens');
                ipar= strcmp(g_sens.pars,par1);
                ivar=str2double(var(2:end));                
                var=i_statevars_names(g_sens.maxndx(ivar,ipar));
            end

            p1=[p1(1:f1-1) {var} p1(f2+1:end)];
        end

        s=sprintf('%s',p1{:});
    end

end

specials={'Delta','Gamma','Lambda','Omega','Phi','Pi','Psi','Sigma',...
   'Theta','Upsilon','Xi','alpha','beta','chi',...
   'delta','epsilon','eta','gamma','iota','kappa','lambda','mu','ni','nu',...
   'omega',  'phi','pi','psi','rho','sigma','tau','theta','upsilon',...
   'varsigma','vartheta','xi','zeta'};

nspec = 37; %length(specials);
syms=i_symvar(s,{'x','y','z','i','j'});

% suppress a warning if the string ends with '_'
f = strfind(s, '_');
if ~isempty(f) && (f(end)==length(s))
   s = [s(1:end - 1) '\_'];
end


for j = 1:length(syms)
   i = 1;
   char = syms{j}(1);
   while (i < nspec) && (specials{i}(1) < char)
      i = i + 1;
   end

   while (i<=nspec) && (specials{i}(1)==char)
      %   for i=1:length(specials)
      if strcmp(specials{i}, syms{j})
         s=strrep(s,specials{i},['{\' specials{i} '}']);
      end

      i = i + 1;
   end

   %  end

end


f = strfind(s, 'val(''');
if ~isempty(f)
   f2=strfind(s(f(1):end), '''');
   if length(f2)>=2
       f2=f2+f(1)-1;
       s = [s(1:f(1)-1) s(f2(1) + 1:f2(2) - 1) '_0' s(f2(2) + 2:end)];
   end

end

f = strfind(s, '_statevar(');
if ~isempty(f)
   s = [s(1:f(1)) ' of ' s(f(1) + 10:length(s) - 1)];
   s=strrep(s,'''','');
end

f = strfind(s, 'outf(');
if ~isempty(f)
   if f(1)>1
       s1=s(1:f(1)-1);
   else
       s1='';
   end

   f2 = strfind(s,',');
   if ~isempty(f2)
      f3=f(1)+6;
      while f3<length(s) && s(f3)~=''''
          f3=f3+1;
      end

      f4=f2(1)+2;
      while f4<length(s) && s(f4)~=''''
          f4=f4+1;
      end
      fun = s(f(1) + 6:f3);      
      if strncmp(fun,'gini',4)
          fun='Gini coefficient';
      end

      if strncmpi(fun, 'mov_', 4)
         fun = sprintf('moving %s', fun(5:end));
      elseif strncmpi(fun, 'perc', 4)
         fun = sprintf('%s%% percentile', fun(5:end));
      elseif strncmpi(fun, 'cover', 5)
         if (length(fun)>5)&&((fun(6)=='>')||(fun(6)=='<'))
            op = fun(6);
         else
            op = '>';
         end

         if length(f2) == 2
            aval = s(f2(2) + 1:f3 - 1);
            s = s(1:f2(2));
         else
            aval = '0.01';
         end

         fun = sprintf('cover%s%s', op, aval);
      end


      if f4==length(s)-2
         s = [s1 fun ' of ' s(f2(1) + 1:f4)];
      else
         s = sprintf('%s %s of %s %s',s1,fun,s(f2(1) + 1:f4),s(f4+2:end));
      end
   elseif strcmpi(s,'outf(''domeigen'')')
      s='dominant eigenvalue';
   elseif strcmpi(s,'outf(''realeigen'')')
      s='real eigenvalues';
   elseif strcmpi(s,'outf(''imageigen'')')
      s='imaginary eigenvalues';
   end

   s=strrep(s,'''','');
end


%this function always checks whether s2 is in s1


