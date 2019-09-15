function [ggrind, parsedmodel, parsedspec] = i_parsedmodel(amodel)
%analyse the model as typed by the user to extract all parameters
%state variables, auxiliaries etc.
%The auxiliary variables are sorted based on dependency if necessary
%The odefile is constucted and the g_grind stucture filled.
%first parse is rather slow, the subsequent analysis is very fast and efficient.

%constants
vr = parsed_equation.vr;
keywords={'t', 'pi',  'inf','Inf', 'nan','NaN','true','false','elseif','else','if','end','while','for',...
   'switch','otherwise','case','global','function','return','definepars','defextern','defpermanent',...
   'lag', 'dwiener',    'rednoise', 'definespace', 'rand',     'randn',      'randlogn',  'randi',   'randperm',...
   'boxcartrain','setevent','implicitdisperse'};
typkw = [vr.t, vr.kw, vr.kw, vr.kw, vr.kw, vr.kw, vr.kw, vr.kw, vr.kw, vr.else, vr.if, vr.block2, vr.block1, vr.block1, ...
   vr.block1,   vr.kw,    vr.kw,vr.kw, vr.kw,     vr.kw,vr.definepars,vr.defextern,vr.defpermanent, ...
   vr.lag,vr.dwiener,vr.rednoise,vr.definespace, vr.stochfun,vr.stochfun,vr.stochfun,vr.stochfun, vr.stochfun, ...
   vr.boxcartrain, vr.setevent, vr.implicitdisperse];
ggrind.model = amodel;
ggrind.locfunc = {};
ggrind.locfuncname = {};
ggrind.solver.nonautonomous  = false;
ggrind.solver.isdiffer  = false;
ggrind.solver.haslags = false;
ggrind.statevars.names = {};
ggrind.statevars.vectnames = {};
ggrind.statevars.dims = {};
ggrind.statevars.vector = false;
ggrind.externvars = {};
ggrind.permanent = {};
ggrind.funcnames.names = {};
ggrind.funcs = '';
ggrind.parsedmodel = {};
ggrind.theodefile = {};
ggrind.solver.name = 'ode45';
ggrind.solver.isstochastic = false;
ggrind.boxcar.names = {};
ggrind.boxcar.gcycl = [];
aparsedmodel = cell(size(amodel));
aparsedspec = cell(size(amodel));
maxlen = 0;
infunction = 0;
ggrind.locfunc = {};
ggrind.locfuncname = {};
ggrind.errors = {};
ggrind.solver.haslags = true;
k = 1;
k1 = 0;
%parse and remove local functions
for i = 1:length(amodel)
   %parsing, slowest step:
   obj = parsed_equation(strtrim(amodel{i}));
   aparsedmodel{i} = obj.fs;
   aparsedspec{i} = obj.types;
   %   [aparsedmodel{i}, aparsedspec{i}] = parseeq(strtrim(amodel{i}));
   if ~isempty(aparsedmodel{i})&&strcmp(aparsedmodel{i}{1}, 'function')
      infunction = true;
      f=find(strcmp(aparsedmodel{i}, '='));
      k1 = k1 + 1;
      if ~isempty(f)
         ggrind.locfuncname{k1, 1} = aparsedmodel{i}{f + 1};
      else
         ggrind.locfuncname{k1, 1} = '';
      end

   end

   if ~infunction
      if length(aparsedspec{i}) > maxlen
         maxlen = length(aparsedspec{i});
      end

   else
      ggrind.locfunc{k, 1} = amodel{i};
      k = k + 1;
      if strcmp(aparsedmodel{i}{1}, 'return')
         infunction = false;
      end

      aparsedmodel{i} = {''};
      aparsedspec{i} = vr.empty;
   end

end

%merge '...' and remove empty lines
for i = length(aparsedmodel):-1:1
   if (length(aparsedmodel{i})==1 && isempty(aparsedmodel{i}{1})) ||isempty(aparsedmodel{i})
      aparsedmodel(i) = [];
      aparsedspec(i) = [];
   elseif i < length(aparsedmodel) && strcmp('...', aparsedmodel{i}{end})
      aparsedmodel{i} = [aparsedmodel{i}(1:end - 1) aparsedmodel{i + 1}];
      aparsedspec{i} = [aparsedspec{i}(1:end - 1) aparsedspec{i + 1}];
      if length(aparsedspec{i}) > maxlen
         maxlen = length(aparsedspec{i});
      end

      aparsedmodel(i + 1) = [];
      aparsedspec(i + 1) = [];
   end

end

%make a single matrix of the model (for efficiency)
parsedmodel = cell(length(aparsedmodel), maxlen);
parsedmodel(:) = {''};
parsedspec = zeros(size(parsedmodel)) + vr.empty;
for i = 1:length(aparsedmodel)
   parsedmodel(i, 1:length(aparsedmodel{i})) = aparsedmodel{i};
   parsedspec(i, 1:length(aparsedmodel{i})) = aparsedspec{i};
end

[index] = makeindex(parsedmodel, parsedspec, vr);
parsedspec(findstrcmpndx(parsedmodel, ':',index)) = vr.colon;
parsedspec(findstrcmpndx(parsedmodel, '(',index)) = vr.brack1;
parsedspec(findstrcmpndx(parsedmodel, ')',index)) = vr.brack2;
parsedspec(findstrcmpndx(parsedmodel,',',index)) = vr.comma;
parsedspec(findstrcmpndx(parsedmodel, ' ', index)) = vr.space;
parsedspec(findstrcmpndx(parsedmodel, ';', index)) = vr.semicolon;
parsedspec(findstrcmpndx(parsedmodel, 'for', index)) = vr.block1;

%find elements before equal sign
iseq=strcmpndx(parsedmodel, '=',index);
parsedspec(iseq) = vr.assign;
ieq = sum(iseq, 2);
beforeeq = false(size(parsedspec));
for i = 1:length(ieq)
   if ieq(i) > 0
      fj = find(iseq(i, :));
      beforeeq(i, 1:fj(1)) = true;
      if length(fj) > 1
         fsemicolon=find(parsedspec(i, :) == vr.semicolon);
         for j = 2:length(fj)
            beforeeq(i, fsemicolon(j - 1):fj(j)) = true;
         end

      end

   end

   %f1 = find(iseq(f(i), :), 1);
end

index.fbeforeeq = beforeeq(index.fndx);
%if definepars is used in command style make the arguments variables
[fi, fj] = findstrcmpndx(parsedmodel, 'definepars', index);
[fi2, fj2] = findstrcmpndx(parsedmodel, 'global', index);
fi = [fi; fi2];
fj = [fj; fj2];

for i = 1:length(fi)
   for j = fj(i) + 1:size(parsedmodel, 2)
      if ~isempty(strtrim(parsedmodel{fi(i),j}))&&~any(strcmp(parsedmodel{fi(i),j},{'(',',',')'}))
         parsedspec(fi(i), j) = vr.variable;
         v = parsedmodel{fi(i), j};
         f = strfind(v,',');
         if ~isempty(f)
            f = [0 f length(v) + 1];
            for ii = 2:length(f)
               parsedmodel{fi(i), j + ii - 2} = v(f(ii - 1) + 1:f(ii) - 1);
               parsedspec(fi(i), j + ii - 2) = vr.variable;
            end

            index = makeindex(parsedmodel, parsedspec, vr);
            index.fbeforeeq = beforeeq(index.fndx);
            break;
         end

      end

   end

end

%replace 'dN','/','dt' by  '','N','''
fdt = findstrcmpndx(parsedmodel, 'dt', index, index.fbeforeeq);
if ~isempty(fdt)
   [is, js] = ind2sub(size(parsedspec), fdt);
   for i = 1:length(is)
      var = parsedmodel{is(i), js(i) - 2};
      var = var(2:end);
      parsedmodel(is(i),js(i)-2:js(i))={'',var,''''};
      parsedspec(is(i), js(i) - 2:js(i)) = [vr.empty, vr.variable, 0];
   end

end


fisacc1=findstrcmpndx(parsedmodel,'''',index,index.fbeforeeq);
fisacc=findstrcmpndx(parsedmodel,'''',index);

%isacc=strcmp(parsedmodel,'''');
parsedspec(fisacc) = vr.accent;
parsedspec(fisacc1) = vr.diff;

isdiff = ~isempty(fisacc1);

ft = findstrcmpndx(parsedmodel, 't',index);
ggrind.solver.nonautonomous = false;
if ~isempty(ft)
   parsedspec(ft) = vr.t;
   fdiff = findstrcmpndx(parsedmodel, 't', index, index.fbeforeeq);
   ggrind.solver.isdiffer = ~isempty(fdiff);
   if isdiff&&ggrind.solver.isdiffer
      ggrind.errors{end + 1} = 'Cannot combine difference equation with differential equation';
   end

   if isdiff
      ggrind.solver.nonautonomous = true;
   else
      [lineswt,~] = findstrcmpndx(parsedmodel, 't', index);
      lineswt = unique(lineswt);
      for i = 1:length(lineswt)
         %remove (t+1) or (t-1) only if they are after statevariables
         f1 =  strfind(parsedspec(lineswt(i), :), [vr.brack1 vr.t vr.oper vr.number vr.brack2]);
         for j = 1:length(f1)
            parsedspec(lineswt(i), f1(j):f1(j) + 4) = vr.empty;
         end

         %remove (t)
         f1 =  strfind(parsedspec(lineswt(i), :), [vr.brack1 vr.t  vr.brack2]);
         for j = 1:length(f1)
            parsedspec(lineswt(i), f1(j):f1(j) + 2) = vr.empty;
         end

      end

      parsedspec(fdiff) = vr.diff;
   end

else
   ggrind.solver.isdiffer = false;
end

if ~isdiff&&~ggrind.solver.isdiffer
   ggrind.errors{end + 1} = 'Not any model equations entered';
end


lineswdiff=repmat(any(parsedspec == vr.diff, 2), 1, size(parsedspec, 2));
beforeeqdiff = beforeeq & lineswdiff;

index.utypes = parsedspec(index.fndx(index.undxs1));
fbeforeeqdiff = beforeeqdiff(index.fndx);
fbeforeeq = beforeeq(index.fndx);
varnrs=find(index.utypes >= vr.funcname);
ggrind.explvector = false;
%Main loop to assign the types of variables.
haswrongorder = false(size(index.uelems));

index.uvarfirstass = zeros(size(varnrs)) - 1;
for k = 1:length(varnrs)
   i = varnrs(k);
   %  uvar = index.uelems(i); %convenient with debugging
   andx = transpose(index.undxs1(i):index.undxs2(i));
   avarndx = index.fndx(andx);
   assigned = fbeforeeq(andx);
   if any(assigned)
      index.utypes(i) = vr.auxil;
      iv = index.fi(andx(assigned));
      abeforeeqdiff = fbeforeeqdiff(andx(assigned));
      ms = abeforeeqdiff;
      if any(abeforeeqdiff)
         if any(~abeforeeqdiff)
            disp('Warning: it is not a good practice to assign to state variables, it is better to use auxiliary variables');
         end

         dd.dim1 = 0;
         dd.dim2 = 0;
         for m = 1:length(ms)
            if ms(m)
               index.utypes(i) = vr.statevar;
               elems = parsedspec(iv(m), :);
               elems = elems(beforeeq(iv(m), :));
               selems = parsedmodel(iv(m), :);
               selems = selems(beforeeq(iv(m), :));
               f = strfind(elems, [vr.brack1 vr.number vr.brack2]);
               colons=find(elems == vr.colon);
               isvector =  ~isempty(colons);
               if isvector
                  d.dim1 = str2double(selems{colons(1) + 1});
                  f1=find(elems == vr.brack1);
                  f2=find(elems == vr.brack2);
                  if length(f1) == 1
                     parsedspec(iv(m), f1:f2) = vr.empty;
                  end

                  if length(colons) > 1
                     d.dim2 = str2double(selems{colons(2) + 1});
                  else
                     d.dim2 = 1;
                  end

                  if d.dim1 * d.dim2 == 1
                     isvector = false;
                  end

               elseif ~isempty(f)
                  ggrind.explvector = true;
                  d.dim1 = str2double(selems{f(1) + 1});
                  d.dim2 = 1;
               else
                  d.dim1 = 1;
                  d.dim2 = 1;
               end

               if dd.dim1 < d.dim1
                  dd.dim1 = d.dim1;
               end

               if dd.dim2 < d.dim2
                  dd.dim2 = d.dim2;
               end

            end

         end

         ggrind.statevars.dims{end + 1} = d;
         ggrind.statevars.vector = ggrind.statevars.vector||isvector||ggrind.explvector;
         index.uvarfirstass(i) = iv(1);
      else
         iva = index.fi(andx(~assigned));
         firstuse = min(iva);
         index.uvarfirstass(i) = min(iv);
         parsedspec(avarndx) = vr.auxil;
         if ~isempty(firstuse)&&(firstuse <= index.uvarfirstass(i))
            haswrongorder(i) = true;
         end

      end

   else
      kw = strcmp(index.uelems{i}, keywords);
      if any(kw)
         f = find(kw, 1);
         index.utypes(i) = typkw(f);
         parsedspec(avarndx) = typkw(f);
      else
         if all(parsedspec(avarndx) == vr.funcname)
            if strcmp(index.uelems{i}, ggrind.locfuncname)
               index.utypes(i) = vr.locfun;
            elseif ~isempty(which(index.uelems{i}))
               index.utypes(i) = vr.funcname;
            else
               index.utypes(i) = vr.parameter;
               parsedspec(avarndx) = vr.parameter;
               
            end

         else
            index.utypes(i) = vr.parameter;
            parsedspec(avarndx) = vr.parameter;
         end

      end

   end

end

firststateass=index.uvarfirstass(index.utypes == vr.statevar);
ggrind.statevars.names=transpose(index.uelems(index.utypes == vr.statevar));
[~, ndx1] = sort(firststateass);
if ggrind.statevars.vector
   ggrind.statevars.vectnames = ggrind.statevars.names(ndx1);
   ggrind.statevars.names = {};
   ggrind.statevars.dims = ggrind.statevars.dims(ndx1);
   offset = 0;
   for i = 1:length(ggrind.statevars.dims)
      d = ggrind.statevars.dims{i};
      d.from = offset + 1;
      d.to = offset + d.dim1 * d.dim2;
      offset = d.to;
      ggrind.statevars.dims{i} = d;
   end
   ggrind.statevars.dim = offset;
else
   ggrind.statevars.names = ggrind.statevars.names(ndx1);
   ggrind.statevars.vectnames = {};
   ggrind.statevars.dims = {};
   ggrind.statevars.dim = length(ggrind.statevars.names);
end


if any(index.utypes == vr.defextern)
   [fi, fj]=find(parsedspec == vr.defextern);
   ggrind.externvars = cell(length(fi), 1);
   for i = 1:length(fi)
      pars = getfunctionpars(parsedmodel,parsedspec, fi(i),fj(i),vr);
      f =  strcmp(index.uelems, pars{1});
      if ~isempty(f)
         index.utypes(f) = vr.extern;
      end

      e.name = pars{1};
      if length(pars) > 1
         e.default = pars{2};
         if length(pars) > 2
            e.options = pars(3:end);
         end

      else
         e.default = '0';
      end
      ggrind.externvars{i} = e;
   end

end


if any(index.utypes == vr.defpermanent)
   [fi, fj]=find(parsedspec == vr.defpermanent);
   for i = 1:length(fi)
      pars =  getfunctionpars(parsedmodel,parsedspec, fi(i),fj(i),vr);
      f =  strcmp(index.uelems, pars{1});
      indexcode = index.uelemnr(f);
      parsedspec(index.elems == indexcode)=vr.permanent;
      if ~isempty(f)
         index.utypes(f) = vr.permanent;
         perm.name = pars{1};
         if length(pars) > 1
            perm.currentval = pars{2};
         end

         ggrind.permanent{i} = perm;
         haswrongorder(f) = false;
      end

   end

end

ggrind.definespace = {};
if any(index.utypes == vr.definespace)
   [fi, ~]=find(parsedspec == vr.definespace);
   ggrind.definespace = cell(length(fi), 1);
   for i = 1:length(fi)
      ggrind.definespace{i} = sprintf('%s', parsedmodel{fi(i), :});
      parsedspec(fi(i), :) = vr.empty;
   end

end

%make an unique number for each block of code (only level 1)
startblock=sum(parsedspec == vr.block1 | parsedspec==vr.if, 2);
endblock=sum(parsedspec == vr.block2, 2);
endblock = [0; endblock(1:end - 1)];
blocklevel = cumsum(startblock - endblock);
blocknr = zeros(size(startblock));
bnr = 1;
for i = 1:length(startblock)
   if startblock(i) > 0||blocklevel(i)==0
      bnr = bnr + 1;
   end

   blocknr(i) = bnr;
end

diffblocks=sum(parsedspec == vr.diff, 2);
diffblocks = unique(blocknr(diffblocks > 0));

for i = 1:length(diffblocks)
   blocknr(diffblocks(i) == blocknr)=diffblocks(i) + vr.mindiff;
end

deleterows=sum(parsedspec==vr.definepars | parsedspec==vr.defextern | parsedspec==vr.defpermanent, 2) > 0;
blocknr(deleterows) = 0;
ggrind.pars=index.uelems(index.utypes == vr.parameter);
if size(ggrind.pars, 1) > 1
   ggrind.pars = transpose(ggrind.pars);
end

p = lower(ggrind.pars);
[~, ndx1] = sort(p);
ggrind.pars = ggrind.pars(ndx1);
ggrind.funcnames.names=index.uelems(index.utypes == vr.auxil);
if ~isempty(ggrind.funcnames.names)
   firstuse=index.uvarfirstass(index.utypes == vr.auxil);
   [~, ndx] = sort(firstuse);
   ggrind.funcnames.names = transpose(ggrind.funcnames.names(ndx));
end

%if ~isempty(ggrind.funcnames.names) %for testing
if any(haswrongorder)
   %sort auxvars
   disp('sorting auxiliary variables');
   selaux=parsedspec == vr.auxil;
   %    auxblocks=auxblocks(auxblocks<vr.mindiff);
   assignedaux = index.elems.*(selaux & beforeeq);
   assignedaux=assignedaux(:, sum(assignedaux, 1) ~= 0);
   doubleaux=~([true; diff(assignedaux(:, 1))] | assignedaux(:, 1)==0);
   auxblocks = unique(blocknr(any(selaux, 2) & blocknr < vr.mindiff));
   usedaux = index.elems.*(selaux & ~beforeeq);
   usedaux = sort(usedaux, 2);
   usedaux=usedaux(:, sum(usedaux, 1) ~= 0); %auxilaries that are used per block
   assignedset = 0;
   block = zeros(size(auxblocks));
   k = 0;
   usedblocks = false(size(auxblocks));
   %loop to sort the blocks, loop till nothing is changed
   ischanged = true;
   while ischanged
      ischanged = false; %this is also safe if not all blocks can be assigned
      i = 1;
      while i  <= length(auxblocks)
         usedau=usedaux(auxblocks(i) == blocknr, :);
         if size(usedau, 1) > 1
            %more than one lines in the block, it can be that it is
            %self-assigning (don't check details)
            assignedau=assignedaux(auxblocks(i) == blocknr, :);
            usedau = usedau(~ismember(usedau(:), assignedau(:)));
         end

         usedau = usedau(usedau > 0);
         % a block of code can be added if all used auxiliaries are already
         % defined
         if ~usedblocks(i)&&all(ismember(usedau, assignedset))
            while 1  % if there are two the same assigned auxiliary variables both lines are added
               usedblocks(i) = true;
               ischanged = true;
               k = k + 1;
               block(k) = auxblocks(i);
               assignedset=union(assignedset, assignedaux(auxblocks(i) == blocknr, :));
               if (i==length(auxblocks))||usedblocks(i+1)||~any(doubleaux(auxblocks(i+1)==blocknr))
                  break;
               end

               i = i + 1;
            end

         end

         i = i + 1;
      end

   end

   if all(usedblocks)
      %only if all blocks are assigned we can trust the result
      [~, ndx] = sort(block);
      auxblocks1 = auxblocks(ndx);
      blockndx = any(selaux, 2) & blocknr < vr.mindiff;
      if sum(blockndx) == length(auxblocks) % no blcoks  > 1
         blocknr(blockndx) = auxblocks1;
      else
         blocknr1 = blocknr;
         for i = 1:length(auxblocks)
            blocknr1(blocknr == auxblocks(i)) = auxblocks1(i);
         end

         blocknr = blocknr1;
      end

   else
      %write the partly sorted list of auxiliaries to help finding the
      %circular dependence
      fprintf(2, 'Cannot sort all auxiliary variables:\nThe following auxiliary variables are sorted in sequential dependence:\n');
      f = find(selaux);
      for i = 1:length(f)
         parsedmodel{f(i)}=sprintf('<a href="matlab:%%" style="font-weight:bold">%s</a>', parsedmodel{f(i)});
      end

      usedb=block(block ~= 0);
      if isempty(usedb)
         disp('none');
      else
         for i = 1:length(usedb)
            ndx=find(blocknr == usedb(i));
            for j = 1:length(ndx)
               disp(writeline(ndx(j), parsedspec, parsedmodel, vr));
            end

         end

      end

      fprintf(2,'\nThese cannot be sorted, probably because of circular dependence:\n');
      unusedb = auxblocks(~usedblocks);
      for i = 1:length(unusedb)
         ndx=find(blocknr == unusedb(i));
         for j = 1:length(ndx)
            disp(writeline(ndx(j), parsedspec, parsedmodel, vr));
         end

      end

      ggrind.errors{end + 1} = 'Cannot sort auxuliary variables, possibly due to circular dependence or double assignments or two statements on one row';
   end

end


% if nargin > 1
%    if g_grind.statevars.vector
%       statenames1 = g_grind.statevars.vectnames;
%       statenames2 = ggrind.statevars.vectnames;
%    else
%       statenames1 = g_grind.statevars.names;
%       statenames2 = ggrind.statevars.names;
%    end

%    haserror=(length(g_grind.pars)~=length(ggrind.pars))|| any(~strcmp(g_grind.pars, ggrind.pars));
%    haserror= haserror || (length(statenames1)~=length(statenames2))||any(~strcmp(statenames1, statenames2));
%    if 0 % the current grind has too many errors to check this
%       uoldfuncs = unique(g_grind.funcnames.names); %old names are not unique
%       unewfuncs = sort(ggrind.funcnames.names);
%       if ~isempty(uoldfuncs)&&~isempty(unewfuncs)
%          haserror= haserror || (length(uoldfuncs)~=length(unewfuncs))||any(~strcmp(uoldfuncs, unewfuncs));
%       end

%    end
%    if haserror
%       error('test:analysemodel','params not correct');
%    end

% end


% ndx1 = ndx(b < vr.mindiff);
% funcs = cell(size(ndx1));
% for i = 1:length(ndx1)
%   funcs{i} = printline(parsedmodel,parsedspec,ndx1(i),vr);
% end

% ggrind.funcs = sprintf('%s\n', funcs{:});
parsedspec(parsedspec == vr.diff)=vr.empty;
%ggrind.parsedmodel = parsedmodel(ndx, :);
[b, fndx] = sort(blocknr);
fndx=fndx(b ~= 0 & b < vr.mindiff);
ggrind.funcs = cell(length(fndx), 1);
kf = 1;
for i = 1:size(fndx, 1)
   s1 = writeline(fndx(i), parsedspec, parsedmodel, vr);
   if ~isempty(s1)
      ggrind.funcs{kf} = s1;
      kf = kf + 1;
   end

end

ggrind.funcs = sprintf('%s\n', ggrind.funcs{1:kf - 1});
%Adapt the code of parsedmodel such that we can write an odefile:
ggrind = makeodefile(ggrind, parsedmodel, parsedspec, blocknr, index, vr, beforeeqdiff);

function ggrind = makeodefile(ggrind, parsedmodel, parsedspec, blocknr, index, vr, beforeeqdiff)
fbeforeeqdiff = beforeeqdiff(index.fndx);
ggrind.solver.name = 'ode45';
odefile = cell(length(parsedmodel) + 20, 1);
fparsedspec = parsedspec(index.fndx);
felemnr = index.elems(index.fndx);
k = 1;
odefile{k} = '%function created by GRIND'; k = k + 1;
odefile{k}='function g_X2=curr_ode0(t,g_X1)';kfun=k; k=k + 1;
odefile{k} = ''; kglob = k; k = k + 1;
sglobals = {};
if ~isempty(ggrind.externvars) ||~isempty(ggrind.permanent)||ggrind.statevars.vector
   sglobals = {'g_grind'};
   if ~isempty(ggrind.permanent)
      odefile{k} = 'i_updatepermanent(t);';k=k + 1;
      for i = 1:length(ggrind.permanent)
         p = ggrind.permanent{i};
         iuelem = strcmp(p.name, index.uelems);
         elnr = index.uelemnr(iuelem);
         parsedmodel(index.fndx(felemnr == elnr)) ={sprintf('g_grind.permanent{%d}.currvalue', i)};
         %        permndx = strcmpndx(parsedmodel, p.name, index);
         %        parsedmodel(permndx) = {sprintf('g_grind.permanent{%d}.currvalue', i)};
      end

   end
   for i = 1:length(ggrind.externvars)
      e = ggrind.externvars{i};
      odefile{k} = sprintf('%s=externvar(%d,%s,t);',e.name,i,e.default);k = k + 1;
   end

end

%analyse lags
flagndx=fparsedspec == vr.lag;
ggrind.solver.haslags = any(flagndx);
if ggrind.solver.haslags
   if ggrind.statevars.vector
      ggrind.errors{end + 1} = 'Lags are not yet supported in vector models';
   end

   if ggrind.solver.isdiffer
      ggrind.errors{end + 1} = 'Lags are not supported in difference equations';
   end

   %[lagi, lagj] = findstrcmpndx(parsedmodel, 'lag', index);
   lagi = index.fi(flagndx);
   lagj = index.fj(flagndx);
   ggrind.dde.lags = {};
   for i = 1:length(lagi)
      [pars,js] =  getfunctionpars(parsedmodel,parsedspec, lagi(i),lagj(i),vr);
      avar = pars{1};
      isstrng = false;
      if  strcontains(avar,'''')
         avar = avar(2:end - 1);
         isstrng = true;
      end

      ivar =  find(strcmp(avar, ggrind.statevars.names));
      parsedmodel{lagi(i), js(1,1)} = int2str(ivar);
      if ~isstrng
         parsedspec(lagi(i), js(1,1))  = 0;
         index.elems(lagi(i), js(1, 1)) =  max(max(index.elems)) + 1;
         felemnr = index.elems(index.fndx);
      end

      ggrind.dde.lags{end + 1} = sprintf('%s', pars{2});
      parsedmodel{lagi(i), js(2,1)} = int2str(length(ggrind.dde.lags));
      parsedmodel{lagi(i), lagj(i)} = 'Lags';
   end

   ggrind.solver.name = 'ddesol';
   if ~any(strcmp('g_grind', sglobals))
      sglobals = [sglobals {'g_grind'}];
   end

   odefile{kfun}='function g_X2=curr_ode0(t,g_X1,Lags)';
   odefile{k}=sprintf('if nargin==2\n  Lags=repmat(g_X1,length(g_grind.dde.lags));\nend;');k=k + 1;
end

%boxcartrain
%example:
% juv_=boxcartrain(juv,devjuv,sd);
% [juv_,juv]=boxcartrain(1,juv,devjuv,sd);
fbox=fparsedspec == vr.boxcartrain;
if any(fbox)
   is = index.fi(fbox);
   js = index.fj(fbox);
   ggrind.boxcar.names = cell(1, length(is));
   ggrind.boxcar.gcycl = zeros(1, length(is));
   for i = 1:length(is)
      pars =  getfunctionpars(parsedmodel,parsedspec, is(i),js(i),vr);
      j = 1;
      while j < size(parsedspec, 1)&&parsedspec(is(i), j)~=vr.auxil
         j = j + 1;
      end

      parsedmodel{is(i),j} = sprintf('[%s,%s]',parsedmodel{is(i),j},pars{1});
      parsedmodel{is(i),js(i) + 1} = sprintf('(%d,',i);
      ggrind.boxcar.names{i} = pars{1};
   end

end

%setevent
if any(fparsedspec == vr.setevent)
   ggrind.errors{end + 1} = '"setevent" should be in the lower parameters panel';
end

%implicitdisperse
if ggrind.statevars.vector
   fimplicit=fparsedspec == vr.implicitdisperse;
   if any(fimplicit)
      is = index.fi(fimplicit);
      js = index.fj(fimplicit);
      ggrind.implicdisp = cell(1, length(js));
      for i = 1:length(is)
         [pars,js1] =  getfunctionpars(parsedmodel,parsedspec, is(i),js(i),vr);
         parsedmodel{is(i), js1(1, 1)} = int2str(i);
         parsedmodel(is(i), js1(1, 1) + 1:js1(1, 2)) = {''};
         parsedmodel{is(i), js1(2, 1)}  = pars{1};
         parsedmodel(is(i), js1(2, 1) + 1:js1(2, 2)) = {''};
         p.Name = pars{1};
         p.D = pars{2};
         ggrind.implicdisp{i} = p;
      end

      ggrind.solver.name = 'euler';
   end

end

%dwiener
fdwiener=fparsedspec == vr.dwiener;
if any(fdwiener)
   is = index.fi(fdwiener);
   js = index.fj(fdwiener);
   for i = 1:length(is)
      [~,js1] =  getfunctionpars(parsedmodel,parsedspec, is(i),js(i),vr);
      parsedmodel{is(i),js1(1,1)} = sprintf('t,%s',parsedmodel{is(i),js1(1,1)});
   end

end

globalpar = true;
if globalpar
   sglobals = [ggrind.pars sglobals];
else
   for i = 1:length(ggrind.pars)
      fndx = findstrcmpndx(parsedmodel, ggrind.pars{i}, index);
      parsedmodel(fndx) = {sprintf('g_grind.parvalues{%d}', i)};
   end
   %     ppars=index.fndx(fparsedspec == vr.parameter);
   %     for i = 1:length(ppars)
   %        parsedmodel{ppars(i)} = sprintf('g_grind.parvalues.%s', parsedmodel{ppars(i)});
   %     end

   if ~any(strcmp('g_grind', sglobals))
      sglobals = [sglobals {'g_grind'}];
   end
end

%stochastic functions
if ~isempty(intersect(fparsedspec,[vr.dwiener,vr.rednoise, vr.stochfun, vr.boxcartrain]))
   ggrind.solver.isstochastic = true;
   ggrind.solver.name = 'euler';
end

if ggrind.solver.isdiffer
   ggrind.solver.name = 'i_differ';
end


if ggrind.statevars.vector
   odefile{k}='g_X2=zeros(g_grind.statevars.dim,1);'; k=k + 1;
   for i = 1:length(ggrind.statevars.vectnames)
      if  ggrind.statevars.dims{1}.dim2 == 1
         odefile{k}= sprintf('%s = g_X1(g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to);',ggrind.statevars.vectnames{i}, i, i);
      else
         odefile{k}= sprintf('%s = reshape(g_X1(g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to),g_grind.statevars.dims{%d}.dim1,g_grind.statevars.dims{%d}.dim2);', ...
            ggrind.statevars.vectnames{i}, i, i, i ,i);
      end

      k = k + 1;
      iuelem = strcmp(ggrind.statevars.vectnames{i}, index.uelems);
      elnr = index.uelemnr(iuelem);
      if ggrind.explvector
         andx=felemnr==elnr & fbeforeeqdiff;
         iv = index.fi(andx);
         jv = index.fj(andx);
         for j = 1:length(iv)
            f = strfind(parsedspec(iv(j),:), [vr.brack1 vr.number vr.brack2]);
            if ~isempty(f)
               parsedspec(iv(j), f(1):f(1) + 2) = vr.empty;
               aval = str2double(parsedmodel(iv(j), f(1) + 1));
               parsedmodel{iv(j),jv(j)} = sprintf('g_X2(g_grind.statevars.dims{%d}.from+%d,1)',i,aval - 1);
            end

         end

      else
         parsedmodel(index.fndx(felemnr==elnr & fbeforeeqdiff)) ={sprintf('g_X2(g_grind.statevars.dims{%d}.from:g_grind.statevars.dims{%d}.to)', i, i)};
      end

   end

else
   for i = 1:length(ggrind.statevars.names)
      iuelem = strcmp(ggrind.statevars.names{i}, index.uelems);
      elnr = index.uelemnr(iuelem);
      parsedmodel(index.fndx(felemnr==elnr & fbeforeeqdiff)) ={sprintf('g_X2(%d,1)',i)};
      parsedmodel(index.fndx(felemnr==elnr & ~fbeforeeqdiff)) ={sprintf('g_X1(%d)', i)};
   end

end

[b, ndx] = sort(blocknr);
ndx=ndx(b ~= 0);
for i = 1:size(ndx, 1)
   s1 = writeline(ndx(i), parsedspec, parsedmodel, vr);
   if ~isempty(s1)
      odefile{k} = s1;
      k = k + 1;
   end

end

if ~isempty(sglobals)
   odefile{kglob} = sprintf('global %s;', strtrim(sprintf('%s ',sglobals{:})));
end

if ~isempty(ggrind.permanent)
   odefile{k} = 'i_updatepermanent;';k=k + 1;
end


ggrind.theodefile = odefile(1:k - 1);

function s = writeline(i, parsedspec, parsedmodel, vr)
s2 = parsedspec(i, :);
ndxs2=(s2~=vr.empty & s2~=vr.comment);
f = find(ndxs2, 1, 'last');
if s2(ndxs2(f)) == vr.semicolon
   ndxs2(f) = false;
end

s1 = parsedmodel(i, ndxs2);
if ~isempty(s1)
   s2 = s2(ndxs2);
   %  if ~isempty(intersect([vr.else,vr.if],s2))]
   if any(s2==vr.else | s2==vr.if)
      s=sprintf('%s',sprintf('%s',s1{:}));
   else
      s=sprintf('%s;',sprintf('%s',s1{:}));
   end

else
   s = '';
end

function index = makeindex(parsedmodel, parsedspec, vr)
[b, ndx] = sort(parsedmodel);
if size(b, 1) == 1
   b = transpose(b);
   ndx = transpose(ndx);
end

nempty=parsedspec(ndx) ~= vr.empty & parsedspec(ndx) ~= vr.comment;
b = b(nempty);
index.fndx = ndx(nempty);
[index.fi, index.fj] = ind2sub(size(parsedspec), index.fndx);
ddiff = ~strcmp(b(1:end - 1), b(2:end));
sameelems = cumsum([1; ddiff(1:end)]);
index.uelems = b([ddiff; true]);
index.undxs1 = find([true; ddiff]);
index.undxs2 = find([ddiff; true]);
index.uelemnr = transpose(1:length(index.undxs1));
index.elems = zeros(size(parsedmodel));
index.elems(index.fndx) = sameelems;

function [ndx1, ndx2] = findstrcmpndx(~, elem, index, afilter)
iuelem = find(strcmp(elem, index.uelems));
if ~isempty(iuelem)
   if nargin == 4
      endx = false(size(index.fndx));
      endx(index.undxs1(iuelem):index.undxs2(iuelem)) = true;
      endx(~afilter) = false;
      if nargout == 2
         ndx1 = index.fi(endx);
         ndx2 = index.fj(endx);
      else
         ndx1 = index.fndx(endx);
      end

   else
      if nargout == 2
         ndx1 = index.fi(index.undxs1(iuelem):index.undxs2(iuelem));
         ndx2 = index.fj(index.undxs1(iuelem):index.undxs2(iuelem));
      else
         ndx1 = index.fndx(index.undxs1(iuelem):index.undxs2(iuelem));
      end

   end

else
   ndx1 = [];
   ndx2 = [];
end



function [ndx] = strcmpndx(parsedmodel, elem, index)
iuelem = strcmp(elem, index.uelems);
elemnr = index.uelemnr(iuelem);
if ~isempty(elemnr)
   ndx=index.elems == index.uelemnr(iuelem);
else
   ndx = false(size(parsedmodel));
end


function [pars,js] =  getfunctionpars(parsedmodel,parsedspec, i,j,vr)
if parsedspec(i, j + 1) == vr.space
   %command style function
   specs = parsedspec(i, j + 1:end);
   js=find(specs == vr.space);
   jend=find(specs == vr.semicolon);
   if isempty(jend)
      jend = length(specs);
      while jend > 0&&specs(jend)==vr.empty
         jend = jend - 1;
      end

      jend = jend + 1;
   end

   jend = jend + j;
   js = js + j;
   js = [js'+1,[js(2:end)'; jend] - 1];
   pars = cell(size(js, 1), 1);
   for k = 1:size(js, 1)
      pars{k} = sprintf('%s', parsedmodel{i, js(k, 1):js(k, 2)});
      f= strfind(pars{k},'''');
      if length(f) == 2
         pars{k} = pars{k}(f(1) + 1:f(2) - 1);
      end

   end

else
   %brackets
   j = j + 2;
   poplevel = 1;
   pars = cell(5, 1);
   js = zeros(5, 2);
   ip = 1;
   while (j < size(parsedspec, 2))&&poplevel > 0
      k = 0;
      while (j+k < size(parsedspec, 2)) && poplevel~=0 && ~(poplevel==1&&any(parsedspec(i, j+k)== vr.comma))
         if parsedspec(i, j + k) == vr.brack1
            poplevel = poplevel + 1;
         end

         k = k + 1;
         if any(parsedspec(i, j + k) == [vr.semicolon, vr.brack2])
            poplevel = poplevel - 1;
         end

      end

      js(ip, :) = [j, j + k - 1];
      pars{ip} = sprintf('%s', parsedmodel{i, js(ip, 1):js(ip, 2)});
      ip = ip + 1;
      j = j + k + 1;
   end

   pars = pars(1:ip - 1);
   js = js(1:ip - 1, :);
end



