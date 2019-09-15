function  [varlist1,otherswitches] = i_getoutlist(varlist)
global g_grind g_data;
varlist1 = cell(1, length(varlist));
otherswitches={};
k = 1;
if nargin > 0
   if ~iscell(varlist)
      varlist = {varlist};
   end

   for i = 1:length(varlist)
      if isempty(varlist{i})|| (varlist{1}(1) ~= '-')
         [varlist1, k] = AddFun(varlist{i}, varlist1, k);
      else
         outfun = getoutf(varlist{i});
         if ~isempty(outfun)
            if ~g_grind.statevars.vector
               for h = 1:g_grind.statevars.dim
                  varlist1{k} = i_statevars_names(h);
                  k = k + 1;
               end

            else
               for i1 = 1:length(g_grind.statevars.vectnames)
                  varlist1{k} = sprintf('outf(''%s'',''%s'')',outfun,g_grind.statevars.vectnames{i1});
                  k = k + 1;
               end

            end

         elseif strcmp(varlist{i}, '-extern')
            for h = 1:length(g_grind.externvars)
               varlist1{k} =  g_grind.externvars{h}.name;
               k = k + 1;
            end

         elseif strcmp(varlist{i}, '-data')
            if ~isempty(g_data)&&(isfield(g_data,'obs') && ~isempty(g_data.obs))
               for h = 1:length(g_data.varlist)
                  if ~min(isnan(g_data.obs(:, h)))
                     varlist1{k} =  sprintf('observed(''%s'')', g_data.varlist{h});
                     k = k + 1;
                  end

               end

            end

         elseif  strncmpi(varlist{i}, '-def',4)
               otherswitches= [ otherswitches varlist(i)];
         else
            [varlist1, k] = AddFun(varlist{i}, varlist1, k);
         end

      end

   end

   if k > 1
      varlist1 = varlist1(1:k - 1);
   else
      varlist1  = {};
   end

end

function [varlist1, k] = AddFun(fun, varlist1, k)
global g_grind;
if g_grind.statevars.vector
   f=strcmp(fun, g_grind.statevars.vectnames);
   if any(f)
       h=find(f,1);
       siz=[g_grind.statevars.dims{h}.dim1,g_grind.statevars.dims{h}.dim2];
       psiz=prod(siz);
       varlist1(k:k+psiz-1)=reshape(allelems(g_grind.statevars.vectnames{h},siz),1,psiz);
       k=k+psiz;
      return;
   end

end

varlist1{k} = fun;
k = k + 1;

function outfun = getoutf(s)
funcs={'_min','_max','_mean','_sum','_median','_std','_cv','_variance','_var','_perc','_shannon'};
outfun = '';
for i = 1:length(funcs)
   fun = funcs{i};
   fun(1) = '-';
   if strncmpi(s, fun, length(fun))
      if strcmpi(fun, '-perc')
         outfun = s(2:end);
      else
         outfun = fun(2:end);
      end

      return;
   end

end


