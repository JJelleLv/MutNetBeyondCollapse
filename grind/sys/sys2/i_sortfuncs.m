function res = i_sortfuncs(funcs1)
global g_grind;
if isempty(funcs1)
   res = funcs1;
   return;
end

if ischar(funcs1)
    funcs1=str2cell(funcs1);
end

funcs = removepoints(funcs1, 0);
permnames = cell(1, length(g_grind.permanent));
for i = 1:length(g_grind.permanent)
   permnames{i} = g_grind.permanent{i}.name;
end

if g_grind.statevars.vector
   filter = [g_grind.pars, g_grind.statevars.vectnames,permnames];%,'g_grind','t'};
else
   filter = [g_grind.pars, g_grind.statevars.names,permnames];%,'g_grind','t'};
end


funnames =  cell(1, length(funcs));
funs = cell(1, length(funcs));
hasif = 0;
for i = 1:length(funcs)
    hasc=strncmp(funcs{i},'   for(int g_ii=0;',18); 
   if ~hasc&&strncmp(funcs{i},'if ',3)||strncmp(funcs{i},'end;',4)||strncmp(funcs{i},'for',3)
      hasif = 1;
   end

   hasbox=strncmp(funcs{i},'   boxcar_result ',17); 
   f=strfind(funcs{i}, '=');
   if ~isempty(f)
      s = funcs{i};
      if hasc
          f1=strfind(funcs{i},'{');
          fun.name=strtrim(s(f1(1)+1:f(2)-1));
          f2=strfind(fun.name,'(');
          if ~isempty(f2)
             fun.name=fun.name(1:f2(1)-1);
          end

      elseif hasbox
          fun.name = strtrim(s(17:f(1) - 1));
      else
          fun.name = strtrim(s(1:f(1) - 1));
      end

      % if ~isempty(setdiff(fun.name,filter)) %not assignment to permanent or statevar
      funnames{i} = fun.name;
      fun.comm = funcs{i};
      if hasc
          s=s(f(2)+1:end);
          s=strrep(s,'(g_ii)','');
      else
      s = s(f(1) + 1:end);
      end

      fun.references = i_symvar(s, filter);
      funs{i} = fun;
      %  end

   else
      if ~hasc&&~(strcmp(funcs{i},'else')||strncmp(funcs{i},'elseif ',7)||strncmp(funcs{i},'if ',3)||strncmp(funcs{i},'end;',4)||strncmp(funcs{i},'for',3))
         warning('GRIND:model:incompleteformula','"%s" incomplete formula?\n', funcs{i});
      end

      fun.references = {};
      funs{i} = fun;
   end

end

if hasif
   res = funcs1;
   return;
end

%remove references to other than functions:
for i = 1:length(funs)
   if ~isempty(funs{i})
      ref = cell(1, length(funs{i}.references));
      l_k = 0;
      for j = 1:length(funs{i}.references)
         found = 0;
         j1 = 1;
         while ~found && (j1<=length(funnames))
            if strcmp(funnames{j1}, funs{i}.references{j})
               found = 1;
            end

            j1 = j1 + 1;
         end

         if found
            l_k = l_k + 1;
            ref{l_k} = funs{i}.references{j};
         end
      end

      funs{i}.references = ref(1:l_k);
   end

end


%Sort

sortedfuncs = cell(1, length(funs));
refs = cell(1, length(funs));
j = 0;
changedsome = 1; %safety for errorous funcs;
while (j < length(funs)) && changedsome
   changedsome = 0;
   for i = 1:length(funs)
      if ~isempty(funs{i}.name) && allunion(funs{i}.references, refs)
         changedsome = 1;
         j = j + 1;
         sortedfuncs{j} = funs{i}.comm;
         refs{j} = funs{i}.name;
         funs{i}.name = '';
      end

   end

end

res =sortedfuncs;

function res = allunion(s1, s2)
res = 1;
for i = 1:length(s1)
   found = 0;
   j = 1;
   while ~found && (j<=length(s2))
      if strcmp(s1{i}, s2{j})
         found = 1;
      end

      j = j + 1;
   end

   if ~found
      res = 0;
      return;
   end

end



function alist = removepoints(alist, docompound)
h = 1;
i = 1;
siz = length(alist);
s = stripcomments(alist);
sfull = alist;
while h <= length(alist)
   alist{i} = '';
   h2=stackedfindend(s,'[',']',h);
   if isempty(h2) || (h2==h)
      h2=stackedfindend(s,'{','}',h);
      if docompound && (isempty(h2) || (h2==h))
         h2=stackedfindend(s,{'for','if','while','switch','case'},'end',h,1);
      end

   end

   %merge h2 lines
   if ~isempty(h2)
      for j = h:h2 - 1
         alist{i} = [char(alist{i}) char(sfull{j}) sprintf('\n') ];  %#ok<SPRINTFN>
         h = h + 1;
      end

   end

   %merge ... lines
   while (h < siz) && (length(strfind(s{h}, '...')) == 1)
      alist{i} = [char(alist{i}) char(sfull{h}) sprintf('\n')];  %#ok<SPRINTFN>
      h = h + 1;
   end

   alist{i} = [char(alist{i}) char(sfull{h})];
   h = h + 1;
   i = i + 1;
end

alist = alist(1: i - 1);

function h = stackedfindend(list, sstart, send, i1,wholewords)
if nargin < 4
   i1 = 1;
end

if nargin < 5
   wholewords = 0;
end

h = i1;
stack = 0;
if isempty(findstrings(sstart, list{h}, wholewords))
   h = [];
else
   found = 0;
   while (h < length(list)) && ~found
      stack = stack + length(findstrings(sstart, list{h},wholewords)) - length(findstrings(send, list{h},wholewords));
      found=stack == 0;
      if ~found
         h = h + 1;
      end

   end

end



function s = stripcomments(s)
%global g_grind;
if iscell(s)
   for j = 1:length(s)
      i = strfind(s{j}, '%'); %strip off comments
      if ~isempty(i)
         s{j} = s{j}(1:i(1) - 1);
      end

   end

else
   i = strfind(s, '%'); %strip off comments
   if ~isempty(i)
      s = s(1:i(1) - 1);
   end

end


