function g_res1 = i_getoutfun(g_afun,g_l_symvar_afun)
global g_t g_Y g_grind g_data;
%as the parameters are declared as globals we need to use
%wierd names for local variables
if nargin<2
    %faster:
   % [g_l_eq,g_l_fn] = parseeq(g_afun);
  %  g_l_symvar_afun=g_l_eq(g_l_fn==50);
   g_l_symvar_afun=symvar(parsed_equation(g_afun));   
   %slower:
    %g_l_symvar_afun = symvar(g_afun2);
end



g_afun=strtrim(g_afun);
if isempty(g_Y)
    g_res1=[];
    return;
end

if strncmpi(g_afun, 'Observed ',9)
    if ~isempty(g_data)
      ivar2 = strtrim(g_afun(10:end));
      if strcmp(ivar2,'t')
          g_res1=g_data.t;
      else 
          indx= strcmp(g_data.varlist,ivar2);
          g_res1 = g_data.obs(:,indx);
      end

    else
      g_res1=[nan;nan];
    end

   return;
end

g_l_X = i_getno(g_afun);
if ~isempty(g_l_X.no) %shortcuts for simple output
   if g_l_X.ispar
      g_res1 = evalin('base', g_afun);
      g_res1 = repmat(transpose(g_res1(:)),size(g_Y, 1), 1);
      return;
 %  elseif g_l_X.isfun
%       i_update_g_func;
%       if isempty(g_l_X.ndx)
%          g_l_i = g_grind.funcnames.dims{g_l_X.vecno}.from:g_grind.funcnames.dims{g_l_X.vecno}.to;
%       else         
%          g_l_i = g_l_X.no;
%       end

%       g_res1 = g_func(:, g_l_i);
   %   error('Auxilary variables are no longer supported for output, use defpermanent instead');
  %    return;
   elseif g_l_X.isvar
      if ~isempty(g_l_X.vecno)&&~strcontains(g_afun, '(')
         g_res1 = g_Y(:,g_grind.statevars.dims{g_l_X.vecno}.from:g_grind.statevars.dims{g_l_X.vecno}.to);
      else 
         g_res1 = g_Y(:, g_l_X.no);
      end

      return;
   elseif g_l_X.isperm
      g_res1 = defpermanent('-g',g_l_X.no);
   end
end

if strcmp(g_afun, 't')
   g_res1 = g_t;
   return;
end

if ~isempty(g_grind.pars)
   eval(i_globalstr(g_grind.pars));
end

[g_afun, ~] = checkvecfuncs(g_afun);

if ~isempty(g_Y)
    if g_grind.statevars.vector
      for g_i=1:length(g_grind.statevars.vectnames)
         eval(sprintf('%s = g_Y(:,%d:%d);',g_grind.statevars.vectnames{g_i},...
            g_grind.statevars.dims{g_i}.from,g_grind.statevars.dims{g_i}.to));
      end

   else
      for g_l_i = 1:g_grind.statevars.dim
         eval([char(i_statevars_names(g_l_i)) '= g_Y(:,' num2str(g_l_i) ');']);
      end

   end
end

if ~isempty(g_grind.externvars)
   for g_l_m = 1:length(g_grind.externvars)
      eval(sprintf('%s=externvar(%d, %s,g_t);',g_grind.externvars{g_l_m}.name, ...
         g_l_m, g_grind.externvars{g_l_m}.default));
   end

end

perms={};
if ~isempty(g_grind.permanent)
   perms=cell(length(g_grind.permanent),1);
   for g_l_m = 1:length(g_grind.permanent)
      eval(sprintf('%s=defpermanent(''-get'',%d);', g_grind.permanent{g_l_m}.name,g_l_m));
      perms{g_l_m}=g_grind.permanent{g_l_m}.name;
   end

end



%try
t = g_t; 
%temporary variable t
%evaluate functions if necessary (it supports matrix notation)
if ~isempty(g_grind) && isfield(g_grind, 'funcnames')&&~isempty(g_grind.funcnames.names) ...
      && isoverlap(g_l_symvar_afun, g_grind.funcnames.names) && ~isoverlap(g_l_symvar_afun,perms)
    if g_grind.statevars.vector
      error('Error: auxilary variables are no longer supported for output \nof vector/matrix models, define them instead \nas permanent variables using "defpermanent"');
%       %  i_evalfuncs; %update g_grind.funcnames
%       global g_func;
%       i_update_g_func;
%       for g_l_i = 1:length(g_grind.funcnames.names)
%          eval(sprintf('%s = g_func(:,%d:%d);',g_grind.funcnames.names{g_l_i}, g_grind.funcnames.dims{g_l_i}.from,g_grind.funcnames.dims{g_l_i}.to));
%       end

    else
       i_evalarray(g_grind.funcs);
    end

end

%catch
%end

g_l_s1 = g_afun;
g_l_i=strfind(g_l_s1, '''');
if ~isempty(g_l_i) && (length(g_l_i) == 1)
   g_N0 = zeros(size(g_Y, 1), g_grind.statevars.dim);
   for g_l_i = 1:g_grind.statevars.dim
      g_N0(:, g_l_i) = g_Y(:, g_l_i);
   end
   NRes = zeros(size(g_Y, 1), g_grind.statevars.dim);
   g_grind.hfun.curr=i_getodehandle(1,'');
   for g_l_i = 1:size(g_Y, 1)
      NRes(g_l_i, :) = transpose(feval(g_grind.hfun.curr, 1, g_N0));
   end

   if ~g_grind.statevars.vector
   for g_l_i = 1:g_grind.statevars.dim
      eval([char(g_grind.statevars.names{g_l_i}) '_ddif = NRes(:,' num2str(g_l_i) ');']);
   end

   else
     error('GRIND:outfun:NotImplemented','_ddif not implemented')
   end

   g_l_s1=strrep(g_l_s1,'''','_ddif');
end

%set all operands to array operands
%g_l_s1=['g_res1=' g_l_s1 ';'];
g_res1=i_evalarray(g_l_s1);
if size(g_res1, 1) == 1
   g_res1 = g_res1 * ones(size(g_Y, 1), 1);
end



function [g_afun1, g_afun2] = checkvecfuncs(g_afun)
global g_grind
g_afun2 = g_afun;
g_afun1 = g_afun;
if g_grind.statevars.vector
   i=length(g_afun1);
   fr=[];
   fc=[];
   fcom=[];
   while i>0
      if g_afun1(i)==''''
          i=i-1;
          while (i>0) && g_afun1(i)~=''''
              i=i-1;
          end

      end

      if g_afun1(i)==')'
         fr=i;
      end

      if g_afun1(i)==':'
         fc=i;
      end

      if g_afun1(i)==','
         fcom=i;
      end

      if (g_afun1(i)=='(')
         fl=i;
         i=i-1;
         while (i>0)&&(isletter(g_afun1(i))||((g_afun1(i)>='0')&&(g_afun1(i)<='9'))||(g_afun1(i)=='_'))
            i=i-1;
         end

         i=i+1;
         if i<fl 
           inam=i_getno(g_afun1(i:fr));
           if ~isempty(inam.no)
              if (isempty(fc)||(fc>fr))&&~inam.ispar
                 if ~isempty(fcom)&&(fcom<fr)
                    g_afun1=sprintf('%s(:,%d)%s',g_afun1(1:fl-1),inam.ndx,g_afun1(fr+1:end));
                 else
                    g_afun1=[g_afun1(1:fl) ':,' g_afun1(fl+1:end)];
                 end

              end

              g_afun2=[g_afun2(1:fl-1) g_afun2(fr+1:end)];
           end

        end

     end

     i=i-1;
  end

end


function found = isoverlap(list1, list2)
%found=isempty(intersect(list1,list2));
 found = 0;
 for i = 1:length(list1)
 %   for j = 1:length(list2)
%       if strcmp(list1{i}, list2{j})
       if any(strcmp(list1{i}, list2))
          found = 1;
          return;
       end

 %   end

 end


