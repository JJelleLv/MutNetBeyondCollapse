%CSOLVERS - create and use fast c++ solvers (ode45 Euler and RK4)
%   CSOLVERS can translate GRIND models to C++ and can compile them to a program including 
%   some solvers (ode45, euler, RK4 and DIFFER for difference equations), it also offers a fast interface
%   to these solvers. The resulting solver can be up to 100 times faster than the MATLAB implementations, 
%   but note that the original MATLAB model can be faster (simple models with big matrices).
%   (depending on the model details).
%  
%    Tips:
%    * On compilation the resulting c++ is tested on consistency. 
%    If this test fails, there is a bug, and the resulting model cannot be used
%    * Not all models are compilable: delay differential equations, DAE and user defined functions cannot 
%    be compiled.
%    * User defined functions can be added as a c++ function, Just use the same name as the MATLAB function
%    (but with .c or .cpp as extension). Only functions on scalars are supported (use double in c++)     
%    * Vector notation is more difficult to compile and the result is not always faster than MATLAB.
%    * In vector/matrix models use dot multiplications wherever possible.
%    * Not all MATLAB functions are (fully) supported (esp. vector/matrix models).
%    * Test the speed/consistency of the resulting model (CSOLVERS -T)
%
%    Usage:
%    CSOLVERS('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'xtra' [general] - extra parameter for option (see '-t')
%    CSOLVERS('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-?' - give information about the current model (read from inifile.exe).
%     '-c' - translate and compile the model. You need to install a C++ compiler to be able to do this 
%	   (see <a href="matlab:commands c-compiler">c-compiler</a>). The resulting executable has the same name as the inifile 
%	   (inifile.exe).
%     '-isok' - quick tests if the MATLAB solver gives the same results.
%     '-l' - translate single line to c++.
%     '-r' - recompile the translated model.
%     '-s' - silent, suppress messages
%     '-t' XTRA - compare the speed of the c.solver with the MATLAB implementation by running for t 
%	   time units (default is defined by <a href="matlab:help simtime">simtime</a>).
%     '-xfuncs' - used internally
%
%  
%    See also c.ode45, c.euler, c.rk4, solver
%
%   Reference page in Help browser:
%      <a href="matlab:commands('csolvers')">commands csolvers</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function [ts,ys,comm] = csolvers(ode, tspan, y0, ~, varargin)
global g_grind g_permanent;
if nargin == 0
   error('csolvers:grind','Not enough input arguments');
end
if ischar(ode) %options including compiling
   switch nargin
       case 1
          argsin={ode};
       case 2
           argsin={ode,tspan};
       case 3
           argsin={ode,tspan,y0};
   end
   fieldnams={'xtra','','extra parameter for option (see ''-t'')',[]}';
   args=i_parseargs(fieldnams,'xtra','-isok,-xfuncs,-s,-r,-t,-?,-l,-c',argsin);
   if any(strcmp(args.opts, '-r'))
      if any(strcmp(args.opts, '-s'))
         silent = 1;
      else
         silent = 0;
      end
      olddir = cd;
      cd(grindpath);
      cd csolver
      disp('compiling');
      [notok, comm] = system('compile.bat');
      if notok
         disp('Error compiling the model (compile.bat)');
         comm=str2cell(comm);
         ndx=~cellfun('isempty',comm);
         comm=comm(ndx);
         ndx=~strcontains(comm,'>REM');
         comm=comm(ndx);
         fprintf('%s\n',comm{:});
         cd(olddir);
         error('compile:csolvers','Compilation of the model unsuccessfull')
      else
         [~, filename] = fileparts(g_grind.inifile);
         filename=checkname(filename);
         if ispc
            filename = [filename '.exe'];
            copyfile('csolver.exe', [olddir filesep filename]);
         elseif ismac % app is the correct extension?
            filename = [filename '.app'];
            copyfile('csolver.app', [olddir filesep filename]);
         end
         if ~silent
            if g_grind.solver.isdiffer
               fprintf('Model compiled successfully\nuse <a href="matlab:solver(''c.differ'');solver(''?'')">solver c.differ</a>  to use this fast solver.\n');
            else
               fprintf('Model compiled successfully\nuse <a href="matlab:solver(''c.ode45'');solver(''?'')">solver c.ode45</a>, <a href="matlab:solver(''c.rk4'');solver(''?'')">solver c.rk4</a> or <a href="matlab:solver(''c.euler'');solver(''?'')">solver c.euler</a> to use these fast solvers.\n');
            end
         else
            disp('Model compiled successfully');
         end
      end
      cd(olddir)
      return;
   elseif any(strcmp(args.opts, '-xfuncs'))  %information of the compiled program
      exefile = getexefile(g_grind.inifile);
      if ~exist(exefile, 'file')
         return;
      end
      [ok, comm] = system(['"' exefile '" ?']);
      if ok == 0
         comm = str2cell(comm);
         nxfunc = str2double(comm{2});
         g_grind.xfuncnames = {};
         g_grind.xfuncs = comm(3:3 + nxfunc - 1);
         for j = 1:nxfunc
            g_grind.xfuncnames{j} = sprintf('g_aux%03d', j);
         end
      end
      return;
   elseif any(strcmp(args.opts, '-isok'))  %check the if curr_ode gives the same result as the c program
       exefile = getexefile(g_grind.inifile);
       comm='';
       if exist(exefile, 'file')
           [ok, comm] = system(['"' exefile '" ?']);
           if ok == 0
               comm = str2cell(comm);
               nxfunc = str2double(comm{2});
               if isnan(nxfunc)
                   nxfunc = 0;
               end
               s = transpose(comm(2 + nxfunc + 1:end));               
               s1=sprintf('%s\n',g_grind.model{:});
               s1=transpose(str2cell(s1));               
               %remove comments
               s=regexp(s, '^[^%]*','match','once');
               s1=regexp(s1, '^[^%]*','match','once');


               if length(s1) ~= length(s)
                   ts = 0;
                   ys = nan;
                   return;
               end
               difcount = ~strcmp(strtrim(s1), strtrim(s));
               if any(difcount)
                   ts = 0;
                   ys = nan;
                   return;
               end
           end
       end          
      N0 = i_initvar;
      ts = 1;
      ys = nan;
      try
         dydt = csolvers([], 1, N0);
      catch err
         disp('This model is not yet compiled');
         disp(err.message)
         ts = 0;
         return;
      end
      dydt1=i_runsinglestep(1,N0',true)';
      same=length(dydt) == length(dydt1);
      if same
         err = sqrt(sum((dydt - dydt1).^2)) / length(dydt);
         if any(isnan(dydt))&&~any(isnan(dydt1))
            fprintf('The c++ version of the model is not working properly: it returns nan''s\n');
            ts = 0;
         elseif (err > 1E-7&&~g_grind.solver.isstochastic)
            fprintf('The model is probably changed (err=%g)\n', err);
            ts = 0;
         elseif g_grind.solver.isstochastic
            disp('Warning: cannot tell if a stochastic model is correct');
         end
         ys = err;
      end
      return
   elseif any(strcmp(args.opts, '-t'))  %test the speed
      global g_Y g_t;  %#ok<TLEV>
      if strcmp(g_grind.solver.name, 'csolvers')
         msolver = g_grind.solver.csolver;
      else
         msolver = g_grind.solver.name;
      end
      csolver = ['c.' msolver];
      if ~isfield(args,'xtra')
         ndays = g_grind.ndays;
      else
         ndays = i_checkstr(args.xtra);
      end
      solver(csolver);
      csolvers('-?')
      solver(msolver);
      fprintf('Running %d time units\n', ndays);
      tic;
      time('-r', ndays);
      t1 = toc;
      y1 = g_Y;
      tt1 = g_t;
      fprintf('Elapsed time with %s is   %.5g seconds\n', msolver, t1);
      solver(csolver);
      tic;
      time('-r', ndays);
      t2 = toc;
      fprintf('Elapsed time with %s is %.5g seconds\n', csolver, t2);
      if t2 < t1
         fprintf('%s is %.4g times faster\n', csolver, t1 / t2);
      else
         fprintf('%s is slower! (%.4g times)\n', csolver, t1 / t2);
      end
      if length(g_t) ~= length(tt1)
         disp('Note: the sizes of the results are different, interpolated results');
         y1 = interp1(tt1, y1, g_t);
      end
      dist = sqrt(sum((g_Y - y1).^2, 2));
      i_figure(100,g_grind.figopts{:}); 
      plot(g_t, dist);
      ylabel('Euclidean distance between runs');
      xlabel('Time (t)');
      fprintf('maximum distance between the runs: %g\n', max(dist));
      return;
   elseif (isfield(args,'xtra')&&strncmp(args.xtra, '?')) ||any(strcmp(args.opts, '-?')) %information of the compiled program
      exefile = getexefile(g_grind.inifile);
      if ~exist(exefile, 'file')
         disp('Executable does not exist, use <a href="matlab:csolvers(''-c'')">csolvers -c</a> to compile model.')
         return;
      end
      [isok, err,comm] =  csolvers('-isok');     
      if ~isempty(comm)
         nxfunc = str2double(comm{2});
         if isnan(nxfunc)
            nxfunc = 0;
         end
         fprintf('%s\n', comm{2 + nxfunc + 1:end});
         disp(' ');
         disp(comm{1});
         s = transpose(comm(2 + nxfunc + 1:end));
         if length(g_grind.model) ~= length(s)
            disp('"Warning: the current inifile is different from this model');
            return;
         end
         difcount = ~strcmp(strtrim(g_grind.model), s);
         if ~any(difcount)
            disp('The compiled model is based on exactly the same definition as the current model');
         else
            fprintf('Warning: the definition of the current model is not the same: %d difference(s)\nline(s) %s\n',sum(difcount),sprintf('%d ',find(difcount)));
         end
      end
      if isok
         fprintf('It gives the same results (err=%g)\n', err);
      end
      return;
   elseif any(strcmp(args.opts, '-l')) % translate a single line to c +  +  (for debugging and userdefined funcs)
      if ~isfield(args,'xtra')
         cmodel = varpaste('cell');
      else
          cmodel=args.xtra;
      end
      if ~iscell(cmodel)
         cmodel = {cmodel};
      end
      userfuncs = getuserfuncs(cmodel);
      includefun = cell(size(userfuncs));
      for i = 1:length(userfuncs)
         includefun{i} = which([userfuncs{i} '.c']);
         if isempty(includefun{i})
            includefun{i} = which([userfuncs{i} '.cpp']);
            if isempty(includefun{i})
               error('grind:csolver','Unknown function in the model: to compile this function please define "%s.c"',userfuncs{i});
            end
         end
      end
      cfuncs = {};
      if g_grind.statevars.vector
         g_grind.xfuncs = {};
         g_grind.xfuncnames = {};
         for i = 1:length(cmodel)
            cmodel{i} = analysevecline(cmodel{i});
         end
         nx2 = length(g_grind.xfuncs);
         cfuncs = cell(1, nx2);
         nx = 0;
         while nx < nx2
            for i = nx + 1:nx2
               cfuncs{i} = analysevecline(g_grind.xfuncs{i});
            end
            nx = nx2;
            nx2 = length(g_grind.xfuncs);
         end
      end
      hasswitch = 0;
      for i = 1:length(cmodel)
         [cmodel{i}, hasswitch] = analysenormline(cmodel{i}, hasswitch);
      end
      hasswitch = 0;
      for i = 1:length(cfuncs)
         [cfuncs{i}, hasswitch] = analysenormline(cfuncs{i}, hasswitch);
      end
      cmodel = [cfuncs cmodel];
      fprintf('%s\n', cmodel{:});
      if nargout == 1
         ts = cmodel;
      end
      return;
   elseif any(strcmp(args.opts, '-c')) %compile odefile
      olddir = cd;
      try
         cd(grindpath);
         cd csolver
         if g_grind.statevars.vector
            warning('vector:csolvers:grind','Vector notation is not always compilable: see <a href="matlab:help csolvers">help csolvers</a> for tips');
            %   in csolvers');
         end
         if g_grind.solver.haslags
            cd(olddir);
            error('lags:csolvers:grind','Time lags not yet supported in csolvers');
         end
         if g_grind.solver.isimplicit
            cd(olddir);
            error('implicit:csolvers:grind','Implicit models (DAEs) not supported in csolvers');
         end
         if g_grind.solver.iscomplex
            cd(olddir);
            error('complex:csolvers:grind','Complex numbers not supported in csolvers');
         end
         fprintf('Translating the model %s to c++\n', g_grind.inifile);
         cmodel = g_grind.model(g_grind.modelndx);
         %remove dots
         for i = 1:length(cmodel)
            cmodel{i}=strrep(cmodel{i},sprintf('...\n'),'');
            cmodel{i}=strrep(cmodel{i},' = ','=');
         end
         cfuncs =  transpose(str2cell(g_grind.funcs));
         %remove externvar(
         for i = length(cfuncs):-1:1
            if ~isempty(strfind(cfuncs{i}, 'externvar('))
               cfuncs(i) = [];
            end
         end
         nmatch = 0;
         cmodel1 = cmodel;
         for i = length(cmodel):-1:1
            %add semicolons
            if isempty(strfind(cmodel{i}, ';'))
               cmodel{i} = sprintf('%s;', cmodel{i});
            end
            %remove funcs (to avoid double funcs)
            if any(strcmp(strtrim(cmodel{i}), cfuncs))
               cmodel(i) = [];
               nmatch = nmatch + 1;
            end
         end
         if nmatch == length(cfuncs)
            cmodel = [cfuncs, cmodel];
         else
            cmodel = cmodel1;
         end
         userfuncs = getuserfuncs(g_grind.model);
         includefun = cell(size(userfuncs));
         newdir = cd;
         cd(olddir);
         for i = 1:length(userfuncs)
            includefun{i} = which([userfuncs{i} '.c']);
            if isempty(includefun{i})
               includefun{i} = which([userfuncs{i} '.cpp']);
               if isempty(includefun{i})
                  error('grind:csolver','Unknown function in the model: to compile this function please define "%s.c"',userfuncs{i});
               end
            end
         end
         cd(newdir);
         if g_grind.solver.isdiffer
            if g_grind.statevars.vector
               if g_grind.statevars.dims{1}.dim2 == 1
                  nam1 = sprintf('%s(1:%d)', g_grind.statevars.vectnames{1}, g_grind.statevars.dims{1}.dim1);
               else
                  nam1 = sprintf('%s(1:%d,1:%d)',g_grind.statevars.vectnames{1},g_grind.statevars.dims{1}.dim1,g_grind.statevars.dims{1}.dim2);
               end
            else
               nam1 = g_grind.statevars.names{1};
            end
            k = strfind(cmodel, [nam1 '(']);
            for i = 1:length(k)
               if ~isempty(k{i})
                  differchar = cmodel{i}(k{i}(1)+length(nam1) + 1);
                  break;
               end
            end
            differ1 = sprintf('(%s+1)', differchar);
            differ2  = sprintf('(%s)', differchar);
            k = strfind(cmodel, differ1);
            found = 0;
            for i = 1:length(k)
               if ~isempty(k{i})
                  found = 1;
               end
            end
            if ~found
               differ1 = sprintf('(%s)', differchar);
               differ2  = sprintf('(%s-1)', differchar);
               k = strfind(cmodel, differ1);
            end
            for i = 1:length(k)
               if g_grind.statevars.vector
                  for j = 1:length(g_grind.statevars.vectnames)
                     cmodel{i} = strrep(cmodel{i},sprintf('%s%s', g_grind.statevars.vectnames{j}, differ2), g_grind.statevars.vectnames{j});
                  end
               else
                  for j = 1:length(g_grind.statevars.names)
                     cmodel{i} = strrep(cmodel{i},sprintf('%s%s', g_grind.statevars.names{j}, differ2), g_grind.statevars.names{j});
                  end
               end
               if ~isempty(k{i})
                  f=strfind(cmodel{i}, '=');
                  avar = strtrim(cmodel{i}(1:k{i} - 1));
                  fbr = strfind(avar, '(');
                  if ~isempty(fbr)
                     avar = avar(1:fbr(1) - 1);
                  end
                  v = i_getno(avar);
                  if v.isvar
                     cmodel{i}=sprintf('%s'' = %s', avar, strtrim(cmodel{i}(f(1) + 1:end))); %add ' to replace later with the dydt
                  end
               end
            end
         end
         if g_grind.statevars.vector
            g_grind.xfuncs = {};
            g_grind.xfuncnames = {};
            for i = 1:length(cmodel)
               cmodel{i} = analysevecline(cmodel{i});
            end
            nx2 = length(g_grind.xfuncs);
            cfuncs = cell(1, nx2);
            nx = 0;
            while nx < nx2
               for i = nx + 1:nx2
                  cfuncs{i} = analysevecline(g_grind.xfuncs{i});
               end
               nx = nx2;
               nx2 = length(g_grind.xfuncs);
            end
         end
         hasswitch = 0;
         for i = 1:length(cmodel)
            [cmodel{i}, hasswitch] = analysenormline(cmodel{i}, hasswitch);
         end
         hasswitch = 0;
         for i = 1:length(cfuncs)
            [cfuncs{i}, hasswitch] = analysenormline(cfuncs{i}, hasswitch);
         end
         for i = length(cmodel):-1:1
            if isempty(cmodel{i})
               cmodel(i) = [];
            end
         end
         ieq=[];
         if isfield(g_grind, 'xfuncs')&&~isempty(cfuncs)
            for i = 1:length(cmodel)
               if ~isempty(strfind(cmodel{i}, '.v_dydt['))
                  ieq = i - 1;
                  break
               end
            end
            if ~isempty(ieq)
               cfuncs = [cfuncs cmodel(1:ieq)];
            end
            for i = length(cfuncs):-1:1
               s = strtrim(cfuncs{i});
               if isempty(s)||strcmp(s, ';')
                  cfuncs(i) = [];
               end
            end
            for i = 1:length(cfuncs)
               if ~isempty(strfind(cmodel{i}, '.v_dydt['))
                  ieq = i - 1;
                  break
               end
            end
          %  cfuncs = i_sortfuncs(cfuncs);
            cmodel = [cfuncs cmodel(ieq + 1:end)];
            %             f = strfind(cmodel, 'g_aux');
            %             xnames = g_grind.xfuncnames;
            %             for i = 1:length(f)
            %                if ~isempty(f(i))&&isempty(strfind(cmodel{i}, '.v_dydt['))
            %                   for j = length(cfuncs):-1:1
            %                      if ~isempty(strfind(cmodel{i}, xnames{j}))
            %                         cmodel = [cmodel(1:i - 1) cfuncs(j) cmodel{i:end}];
            %                         xnames(j) = [];
            %                         cfuncs(j) = [];
            %                      end
            %                   end
            %                end
            %             end
            %             cmodel1 = cmodel;
            %             cmodel = cell(1, length(cmodel1) + length(cfuncs));
            %             h = 0;
            %             for i = 1:length(cmodel1)
            %                if (h==0)&&~isempty(strfind(cmodel1{i}, '.v_dydt['))
            %                   for j = 1:length(cfuncs)
            %                      cmodel{i + j - 1} = cfuncs{j};
            %                      h = h + 1;
            %                   end
            %                   break;
            %                end
            %                cmodel{i + h} = cmodel1{i};
            %             end
            funcnames = unique([g_grind.funcnames.names g_grind.xfuncnames]);
            xfuncs = g_grind.xfuncs;
         else
            funcnames = unique(g_grind.funcnames.names);
            xfuncs = {};
         end
         fid=fopen('model.tmp','w');
         if g_grind.statevars.vector
            fprintf(fid, '#define vectormodel\n');
         end
         for i = 1:length(includefun)
            fprintf(fid, '#include "%s"\n',includefun{i});
         end
         fprintf(fid, '//generated by GRIND on %s\nvoid initialize (TModel & g_model)\n', date);
         
         fprintf(fid, '{\n    g_model.setstatevars(%d);\n    g_model.setparameters (%d);\n', iif(g_grind.statevars.vector, length(g_grind.statevars.vectnames), g_grind.statevars.dim), length(g_grind.pars));
         if ~isempty(g_grind.externvars)
            fprintf(fid, '    g_model.externvars.resize(%d);\n', length(g_grind.externvars));
         end
         if g_grind.solver.isdiffer
            fprintf(fid, '    g_model.setdiffer(true);\n');
         end
         if g_grind.statevars.vector
            fprintf(fid, '    g_model.setauxil(%d); \n',  length(funcnames));
         end
         if isfield(g_grind, 'boxcar') &&~isempty(g_grind.boxcar.names)
            s1 = sprintf('"%s",',g_grind.boxcar.names{:});
            fprintf(fid, '    g_model.boxcar_name={%s};\n',s1(1:end - 1));
         end
         fprintf(fid, '    g_model.setpermanent(%d); \n}\nstring helpmodel()\n', length(g_grind.permanent));
         
         moddescr=strrep(g_grind.model,'\','\\'); %\ is used for escape sequence
         moddescr = sprintf('%s%s\\n',sprintf('%s\\n\n', moddescr{1:end-1}), moddescr{end});
         moddescr=strrep(moddescr,'"','\"'); %" becomes \" in c++
         moddescr=strrep(moddescr,sprintf('\n'),sprintf('"\n  "')); %#ok<SPRINTFN>
         fprintf(fid,'{\n  return "The model ''%s'' was compiled on %s.\\n%d\\n"\n  "%s%s";\n}\n',strrep(g_grind.inifile,'\','\\'),...
            date, length(xfuncs), sprintf('%s\\n', xfuncs{:}), moddescr);
         
         %       if g_grind.solver.isdiffer
         %          for i = 1:length(g_grind.statevars.names)
         %             fprintf(fid,'#define %s(%s) g_X1[%d]\n', g_grind.statevars.names{i}, differchar, i - 1);
         %          end
         %       else
         if g_grind.statevars.vector
            for i = 1:length(g_grind.pars)
               fprintf(fid,'#define %s g_model.v_param[%d]\n', g_grind.pars{i}, i - 1);
            end
            for i = 1:length(g_grind.statevars.vectnames)
               fprintf(fid,'#define %s g_model.v_y0[%d]\n', g_grind.statevars.vectnames{i}, i - 1);
            end
            for i = 1:length(g_grind.permanent)
               fprintf(fid,'#define %s g_model.v_permanent[%d]\n', g_grind.permanent{i}.name, i - 1);
            end
         else
            for i = 1:length(g_grind.pars)
               fprintf(fid,'#define %s g_model.param[%d]\n', g_grind.pars{i}, i - 1);
            end
            for i = 1:length(g_grind.statevars.names)
               fprintf(fid,'#define %s g_X1[%d]\n', g_grind.statevars.names{i}, i - 1);
            end
            for i = 1:length(g_grind.permanent)
               fprintf(fid,'#define %s g_model.permanent[%d]\n', g_grind.permanent{i}.name, i - 1);
            end
         end
         if g_grind.statevars.vector
            for i = 1:length(cmodel) %beautfy the code
               s=strrep(cmodel{i},'}',sprintf('\n   }'));
               s=strrep(s,';;',';');
               cmodel{i}=strrep(s,'{',sprintf('\n   {    \n       '));
            end
            for i = 1:length(funcnames)
               if isempty(strfind(funcnames{i},','))
                  fprintf(fid,'#define %s g_model.v_auxil[%d]\n', funcnames{i}, i - 1);
               end
            end
         end
         %       end
         for i = 1:length(g_grind.externvars)
            fprintf(fid,'#define %s g_model.externvalues[%d]\n', g_grind.externvars{i}.name, i - 1);
            %            fprintf(fid,'#define %s g_model.externvars[%d].valueat(t)\n', g_grind.externvars{i}.name, i - 1);
         end
         
         fprintf(fid,'void odefile (double t, double * g_X1,  double * g_X2,  TModel & g_model)\n{\n');
         %declare local variables
         for i = 1:length(g_grind.permanent)
            p = strcmp(funcnames, g_grind.permanent{i}.name);
            funcnames(p) = [];
         end
         if ~g_grind.statevars.vector&&~isempty(funcnames)
            vars = sprintf('%s,',funcnames{:});
            vars = vars(1:end - 1); %remove last comma
            fprintf(fid, '   double %s;\n', vars);
         end
         ndx = strncmp('   g_X2[', cmodel, 8);
         fprintf(fid, '%s\n', cmodel{~ndx});
         fprintf(fid, '%s\n', cmodel{ndx});
         fprintf(fid, '}\n');
         
         if ~isempty(g_grind.pars)
            fprintf(fid,'#undef %s\n', g_grind.pars{:});
         end
         if g_grind.statevars.vector
            fprintf(fid,'#undef %s\n', g_grind.statevars.vectnames{:});
         else
            fprintf(fid,'#undef %s\n', g_grind.statevars.names{:});
         end
         for i = 1:length(g_grind.permanent)
            fprintf(fid,'#undef %s\n', g_grind.permanent{i}.name);
         end
         for i = 1:length(g_grind.externvars)
            fprintf(fid,'#undef %s\n', g_grind.externvars{i}.name);
         end
         if g_grind.statevars.vector
            for i = 1:length(funcnames)
               if isempty(strfind(funcnames{i},','))
                  fprintf(fid,'#undef %s\n', funcnames{i});
               end
            end
         end
         fclose(fid);
         cd(olddir);
         if nargin == 1
            csolvers('-r');
         else
            csolvers('-r', tspan);
         end
         if ~csolvers('-isok')
            error('csolvers:bugcheck','Don''t use the csolver as there a bug in the compiled model, MATLAB gives a different result');
         end
      catch err
         cd(olddir);
         rethrow(err);
      end
      return;
   end
end
%Start here for applying the compiled model
exefile = getexefile(g_grind.inifile);
i = 1;
parfile = sprintf('%s%sparam%d.tmp', grindpath, filesep, i);
while exist(parfile, 'file')
   i = i + 1;
   parfile = sprintf('%s%sparam%d.tmp', grindpath, filesep, i);
end
outfile = sprintf('%s%sout%d.tmp', grindpath, filesep, i);

%write binary parameter file
signature  = 58275485;
c_param =  1;
c_svar =  2;
c_permanent =  3;
c_auxil =  4;
c_step =  5;
c_reltol =  6;
c_abstol =  7;
c_tstart =  8;
c_tend =  9;
c_nout =  10;
c_srand =  11;
c_outf =  12;
c_solver = 13;
c_tout =  14;
c_extern_data =  15;
c_extern_cycle =  16;
c_extern_active = 17;
c_iters =  18;
c_backw =  19;
c_stats =  20;
c_refine =  21;
c_extern_default = 22;
c_ndays = 23;
c_gcycl = 24;
c_nonnegative = 25;
fid = fopen(parfile, 'wb');
fwrite(fid, signature, 'int32');
if g_grind.solver.isdiffer
   g_grind.solver.opt.StepSize = 1;
end
if isempty(g_grind.solver.opt.StepSize)
   g_grind.solver.opt.StepSize = 0.1;
end
tspan = tspan(:);
if length(tspan) < 2 %this occurs only for singlerun
   fwrite(fid, [c_solver;0], 'int32'); %singlerun
   fwrite(fid, c_tstart, 'int32'); %tstart
   if length(tspan) == 1
      fwrite(fid, tspan, 'double');
   else
      fwrite(fid, 1, 'double');
   end
else
   
   %   if length(tspan) == 2
   %      tspan = (tspan(1):g_grind.solver.opt.StepSize:tspan(2))';
   %   end
   %fprintf(fid, 'outf "%s"\n', outfile);
   
   difftspan = diff(tspan);
   if max(difftspan)-min(difftspan) < 1E-8
      %   fprintf(fid, 'tstart %g\ntend %g\nnout %g\n', tspan(1), tspan(end), length(tspan) - 1);
      fwrite(fid, c_tstart, 'int32'); %tstart
      fwrite(fid, tspan(1), 'double');
      fwrite(fid, c_tend, 'int32'); %tend
      fwrite(fid, tspan(end), 'double');
      if length(tspan) > 2
         fwrite(fid, c_nout, 'int32'); %nout
         fwrite(fid, length(tspan), 'double');
      end
   else
      fwrite(fid, c_tout, 'int32');
      fwrite(fid, length(tspan), 'int32');
      fwrite(fid, tspan, 'double');
      %   fprintf(fid,'tout %s\n',strtrim(sprintf('%g ',tspan)));
   end
   
   if ~isfield(g_grind.solver, 'csolver')
      g_grind.solver.csolver = g_grind.solver.name;
   end
   if strcmp(g_grind.solver.csolver, 'ode45')
      if isempty(g_grind.solver.opt.RelTol)
         g_grind.solver.opt.RelTol = 5E-5;
      end
      if isempty(g_grind.solver.opt.AbsTol)
         g_grind.solver.opt.AbsTol = 1E-7;
      end
      if ~isempty(g_grind.solver.opt.Refine)&&~g_grind.solver.opt.Refine&&(length(tspan)==2)
         fwrite(fid, [c_refine; 0], 'int32');
      end
      fwrite(fid, c_solver, 'int32');
      fwrite(fid, 1, 'int32'); %ode45
      fwrite(fid, c_reltol, 'int32'); %reltol
      fwrite(fid, g_grind.solver.opt.RelTol, 'double');
      fwrite(fid, c_abstol, 'int32'); %abstol
      fwrite(fid, g_grind.solver.opt.AbsTol, 'double');
      %   fprintf(fid, 'solver %s\nreltol %g\nabstol %g\n', ...
      %     g_grind.solver.csolver, g_grind.solver.opt.RelTol, g_grind.solver.opt.AbsTol);
   else
      if ~g_grind.solver.isdiffer
         if strcmp(g_grind.solver.csolver, 'rk4')
            isolver = 2; %rk4
         else
            isolver = 3; %euler
         end
         fwrite(fid, c_solver, 'int32'); %solver
         fwrite(fid, isolver, 'int32');
         fwrite(fid, c_step, 'int32');
         fwrite(fid, g_grind.solver.opt.StepSize, 'double');
         %   fprintf(fid, 'solver %s\nstep %g\n', g_grind.solver.csolver,g_grind.solver.opt.StepSize);
      elseif g_grind.solver.iters > 1
         fwrite(fid, c_iters, 'int32'); %iters
         fwrite(fid, g_grind.solver.iters, 'int32');
         %   fprintf(fid, 'iters %d\n', g_grind.solver.iters);
      end
   end
end
if isfield(g_grind, 'boxcar')&&~isempty(g_grind.boxcar.gcycl)
   fwrite(fid, c_gcycl, 'int32');
   fwrite(fid, length(g_grind.boxcar.gcycl), 'int32');
   fwrite(fid, g_grind.boxcar.gcycl(:), 'double'); %outf
end
fwrite(fid, c_outf, 'int32'); %outf
fwrite(fid, length(outfile), 'int32');
fwrite(fid, outfile, 'char');
fwrite(fid, c_ndays, 'int32'); %reltol
fwrite(fid, g_grind.ndays, 'double');

if g_grind.solver.isstochastic
   fwrite(fid, c_srand, 'int32');
   fwrite(fid, rand(1), 'double');
   %   fprintf(fid, 'srand %0.16e\n', rand(1));
end
if strcmp(g_grind.solver.opt.Stats, 'on')
   fwrite(fid, c_stats, 'int32');
   fwrite(fid, 1, 'int32');
   %    fprintf(fid, 'stats 1\n');
end
if g_grind.statevars.vector
   for i = 1:length(g_grind.statevars.vectnames)
      p = evalin('base', g_grind.statevars.vectnames{i});
      s = size(p);
      fwrite(fid, [c_svar; i - 1; s(:)], 'int32');
      fwrite(fid, p(:), 'double');
      %       fprintf(fid,'svar %d %d %d%s\n', i - 1,s(:),sprintf(' %0.16e',p(:)));
   end
   if ~isfield(g_grind, 'xfuncs')
      csolvers('-xfuncs');
   end
   sizfuncts = sizefuncs([transpose(str2cell(g_grind.funcs)), g_grind.xfuncs], unique([g_grind.funcnames.names g_grind.xfuncnames]));
   for i = 1:length(sizfuncts)
      %p = zeros(sizfuncts{i});
      if prod(sizfuncts{i}) >= 1
         fwrite(fid, [c_auxil; i - 1; sizfuncts{i}(:)], 'int32');
      end
      %         fprintf(fid,'auxil %d %d %d%s\n', i - 1,sizfuncts{i}(:),sprintf(' %d',p(:)));
   end
   for i = 1:length(g_grind.pars)
      p = evalin('base', g_grind.pars{i});
      s = size(p);
      fwrite(fid, [c_param; i - 1; s(:)], 'int32');
      fwrite(fid, p(:), 'double');
      % fprintf(fid,'para %d %d %d%s\n', i - 1,s(:),sprintf(' %0.16e',p(:)));
   end
   for i = 1:length(g_grind.permanent)
      p = g_grind.permanent{i}.currvalue;
      s = size(p);
      fwrite(fid, [c_permanent; i - 1; s(:)], 'int32');
      fwrite(fid, p(:), 'double');
      %    fprintf(fid,'perm %d %d %d %0.16e\n', i - 1,s(:),p);
   end
else
   fwrite(fid, c_svar, 'int32');
   fwrite(fid, length(y0), 'int32');
   fwrite(fid, y0, 'double');
   ppar = ones(length(g_grind.pars), 1);
   for i = 1:length(g_grind.pars)
      p = evalin('base', g_grind.pars{i});
      ppar(i) = p;
   end
   fwrite(fid, c_param, 'int32');
   fwrite(fid, length(ppar), 'int32');
   fwrite(fid, ppar, 'double');
   %    fprintf(fid,'para %d %d %d%s\n', i - 1,s(:),sprintf(' %0.16e',p(:)));
   pperm = ones(length(g_grind.permanent), 1);
   for i = 1:length(g_grind.permanent)
      p = evalin('base', g_grind.permanent{i}.name);
      pperm(i) = p;
   end
   fwrite(fid, [c_permanent; length(pperm)], 'int32');
   fwrite(fid, pperm(:), 'double');
   %    fprintf(fid,'perm %d %d %d %0.16e\n', i - 1,s(:),p);
end
if ~isempty(g_grind.solver.opt.NonNegative)
   fwrite(fid, c_nonnegative, 'int32');
   fwrite(fid, length(g_grind.solver.opt.NonNegative), 'int32');
   fwrite(fid, g_grind.solver.opt.NonNegative, 'int32');
%    if strcmp(g_grind.solver.csolver, 'ode45')
%       warning('GRIND:code45:NonNegativeIgnored','C.ODE45 does not constrain solution to be non-negative. Option ''NonNegative'' will be ignored.');
%    end
end
for i = 1:length(g_grind.externvars)
   if ~isempty(g_grind.externvars{i}.data)
      s = size(g_grind.externvars{i}.data);
      fwrite(fid, [c_extern_data; i - 1; s(:)], 'int32');
      fwrite(fid, g_grind.externvars{i}.data, 'double');
      %       fprintf(fid,'extern.data %d %d %d %s\n', i - 1,s(:),strtrim(sprintf('%.16g ',g_grind.externvars{i}.data(:))));
   end
   fwrite(fid, [c_extern_default; i - 1], 'int32');
   defval = str2num(g_grind.externvars{i}.default);  %#ok<ST2NM>
   if numel(defval) > 1
      error('grind:csolvers','Matrix of vector external variables not yet supported')
   end
   fwrite(fid, defval, 'double');
   if g_grind.externvars{i}.options.cycle
      fwrite(fid, [c_extern_cycle; i - 1; 1], 'int32');
      %       fprintf(fid, 'extern.options.cycle 1\n');
   end
   if ~g_grind.externvars{i}.options.active
      fwrite(fid, [c_extern_active; i - 1; 0], 'int32');
      %       fprintf(fid, 'extern.active 0\n');
   end
end
if g_grind.solver.backwards
   if  g_grind.solver.isdiffer
      error('backw:csolver', 'Backwards not implemented for difference equations, use i_differ instead');
   end
   fwrite(fid, [c_backw; 1], 'int32');
   %     fprintf(fid, 'backw -1\n');
end
fclose(fid);
copyfile(parfile,[grindpath filesep 'csolver' filesep 'param1.tmp']);
system(sprintf('"%s" "%s"', exefile, parfile));
ts = tspan;
%ys = load(outfile);
% read of binary file is muuuch faster than load of ascii file
fid = fopen(outfile, 'r');
signature = fread(fid, 1, 'int32');
if signature == 1234567890
   siz = fread(fid, 2, 'int32');
   if siz(2) <= 0
      siz = double(siz);
      siz(2) = inf;
   end
   ys=fread(fid,transpose(siz),'double');
   ys = transpose(ys);
else
   error('read:csolver','Error reading output file');
end
fclose(fid);
if length(tspan) == 2
   ts = ys(:, 1);
   ys = ys(:, 2:end);
end
if length(tspan) < 2
   ts = transpose(ys);
elseif ~isempty(g_grind.permanent)
   g_permanent.Y = ys(:, end - g_grind.permanent{end}.to + 1:end);
   g_permanent.t = ts;
   n = length(ts);
   if n > 10
      n = 10;
   end
   g_permanent.lasti = n;
   g_permanent.nextt = length(ts) + 1;
   g_permanent.active = 0;
   g_permanent.lastt(1:n) = ts(end - n + 1:end);
   for i = n-1:-1:0
      g_permanent.lastY{i + 1} = ys(end - i, :);
   end
   defpermanent('-s', g_permanent.Y(end, :));
   ys = ys(:, 1:g_grind.statevars.dim);
end

delete(outfile);
%for debugging a copy of the parfile is made to the csolver directory
delete(parfile);

function exefile  = getexefile(inifile)
[~, exefile] = fileparts(inifile);
exefile=checkname(exefile);
if ispc
   exefile = [exefile '.exe'];
   %  s = dir(exefile);
   %  if isempty(s)
   %     gcd(exefile);
   %  end
elseif ismac
   exefile = [exefile '.app'];
   s = dir(exefile);
   if isempty(s)
      gcd(exefile);
   end
end
function exefile=checkname(exefile)
%some names could be regarded as dangerous by Windows
%change to grindexefile
if any(strcmpi(exefile,{'twopatchass'}))
    exefile='grindexefile';
end

function [userfun] = getuserfuncs(smodel)
userfun = {};
k = 1;
for i = 1:length(smodel)
   %      l_function = 30;
   h = parsed_equation(smodel{i});
   f=find(h.types == parsed_equation.vr.funcname);
   % [elems,struc] = parseeq(smodel{i});
   %    isvar = false(size(elems));
   %    for j = 1:length(elems)
   %       isvar(j)=isletter(elems{j}(1))||(elems{j}(1)=='_');
   %    end
   %    f = strcmp(elems, '(') & [0 isvar(1:end - 1)];
   %    f = find(f);
   %    if ~isempty(f)
   %       feq=find(strcmp(elems, '='));
   %       if isempty(feq)
   %          feq = 0;
   %       end
   %       k = 1;
   %       for j = 1:length(f)
   %   f=find(struc==l_function);
   for j = 1:length(f)
      afun = h.fs{f(j)};
      no = i_getno(afun);
      if isempty(no.no)
         a = which(afun);
         if ~strncmp(a, matlabroot, length(matlabroot))&&~strncmp(a, 'built-in ', 9)&&~strncmp(a, grindpath, length(grindpath))
            userfun{k} = afun;
            k = k + 1;
         end
      end
   end
end


function [smodel, hasswitch] = analysenormline(smodel, hasswitch)
smodel = strtrim(smodel);
f = strfind(smodel, '%');
if ~isempty(f)
   smodel = smodel(1:f(1) - 1);
end
hh =  parsed_equation(smodel);
elems = addbrackets(hh);
%replace strings
f=find(strncmp('''',elems,1)&~strcmp('''',elems));
for i = 1:length(f)
   elems{f(i)} = sprintf('"%s"', elems{f(i)}(2:end - 1));
end

elems=cellrep(elems,'.*','*');
elems=cellrep(elems,'./','/');
elems=cellrep(elems,'.^','^');
elems=cellrep(elems,'~','!');
elems=cellrep(elems,'~=','!=');
elems=cellrep(elems,'| ','||');
elems=cellrep(elems,'&','&&');
elems = changeop(elems, {'^','.^'},'power');
fdivs = find(strcmp(elems, '/'));
%check and repair integer divide
for i = 1:length(fdivs)
   if (fdivs(i) < length(elems))&&~isnan(str2double(elems{fdivs(i)+1}))
      if isempty(strfind(elems{fdivs(i) + 1}, '.'))
         elems{fdivs(i) + 1} = [elems{fdivs(i) + 1} '.0'];
      end
   end
end

% [elems, ischanged] = changeop(elems, {'=='},'eq');
% if ischanged
%    [elems] = parseeq(sprintf('%s',elems{:}));
% end
% [elems, ischanged] = changeop(elems, {'~='},'ne');
% if ischanged
%    [elems] = parseeq(sprintf('%s',elems{:}));
% end
isoper = 0;
f = strcmp(elems, 'switch');
if any(f)
   elems=cellrep(elems,';','');
   elems=[{ 'switch',' ', '(','int','(','round','('} elems(2:end) {')',')',')','{'}];
   %   elems{1} = 'switch (int(round(';
   %   elems{end + 1} = '))){';
   hasswitch = 1;
   isoper = true;
end
if hasswitch
   f = strcmp(elems, 'case');
   if any(f)
      elems=cellrep(elems,';','');
      elems = [{'break',';',' '},elems(1:end),{':'}];
      %sprintf('break;%s', elems{1});
      %elems{end + 1} = ':';
   end
   f = strcmp(elems, 'otherwise');
   if any(f)
      elems = {'break',';','default',':'};
   end
end
if any(strcmp(elems,'definepars'))||any(strcmp(elems,'defextern'))||any(strcmp(elems,'defpermanent'))
   elems = {''};
end
f = strcmp(elems, 'val');
while any(f)
   f = find(f, 1);
   f3 = find(strncmp(elems, '"',1));
   avar = strtrim(elems{f3(1)});
   avar = avar(2:end - 1);
   v = i_getno(avar);
   if v.isvar
      elems{f} = sprintf('g_model.y_start[%d]',  v.no - 1);
      f = f + 1;
      while (f < length(elems))&&~strcmp(elems{f}, ')')
         elems{f} = '';
         f = f + 1;
      end
      if strcmp(elems{f}, ')')
         elems(f) = '';
      end
   end
   f = strcmp(elems, 'val');
end
f = strcmp(elems, 'elseif');

if any(f)
   f = find(f, 1);
   elems{f} = '} else if(';
   elems{end + 1} = '){';
   isoper = true;
else
   f = strcmp(elems, 'if');
   if any(f)
      f = find(f, 1);
      elems=cellrep(elems,';',''); %this is not correct to put a ; after if, but MATLAB does not complain
      elems{f} = 'if(';
      elems{end + 1} = '){';
      isoper = true;
   end
   if any(strcmp(elems, 'else'))
      elems = {'} else {'};
      isoper = true;
   end
   if any(strcmp(elems, 'end'))
      elems = {'}'};
      isoper = true;
   end
end
feq = find(strcmp(elems, '='));
if ~isempty(feq)
   f =  strcmp(elems, '«');
   f =  find(~f, 1);
   avar = elems{f};
   v = i_getno(avar);
   if v.isvar
      elems{1} = sprintf('g_X2[%d]', v.no - 1);
      for j = 2:feq(1) - 1
         elems{j} = '';
      end
   end
end

if ~isempty(elems)
   f=strcmp(elems,'«')|strcmp(elems,'»');
   elems = elems(~f);
   if ~isoper
      f = strcmp(elems, ';');
      if ~any(f)
         elems{end + 1} = ';';
      end
   end
   smodel = sprintf('   %s', sprintf('%s',elems{:}));
else
   smodel = '';
end
function [elems] = changeop(elems, opers,func)
f = strcmp(elems, opers{1});
for k = 2:length(opers)
   f =  f | strcmp(elems, opers{k});
end
while any(f)
   brackets=(strcmp(elems,'(')|strcmp(elems,'«'))-(strcmp(elems,')')|strcmp(elems,'»'));
   f = find(f, 1);
   f1 = popbrack(brackets,f - 1, - 1);
   %    h = 0;
   %    while any(strcmp(elems{f+1+h},{'+',' ','-'}))
   %       h = h + 1;
   %    end
   f2 = popbrack(brackets, f + 1, 1);
   elems=[elems(1:f1-1) {func,'('} elems(f1:f-1) {','} elems(f+1:f2) {')'} elems(f2+1:end)];
   f = strcmp(elems, opers{1});
   for k = 2:length(opers)
      f =  f | strcmp(elems, opers{k});
   end
end

function fs =  cellrep(fs, sfrom, sto)
ks = find(strcmp(sfrom, fs));
for i = 1:length(ks)
   fs{ks(i)} = sto;
end

function [smodel] = analysevecline(smodel)
global g_grind;
f = strfind(smodel, '%');
if ~isempty(f)
   smodel = smodel(1:f(1) - 1);
end
%[elems] = parseeq(smodel,1);
obj = parsed_equation(smodel);
elems = addbrackets(obj, true);
[elems] = changeop(elems,{'^','.^'},'power');
elems=cellrep(elems,'[',' ');
elems=cellrep(elems,']',' ');
faccent=strcmp(elems,'''');
feq=strcmp(elems, '=');
if any(faccent)
   f = find(feq);
   f1 = find(faccent);
   if f1(1) < f(1)
      f1(1) = [];
   end
   if length(f1) >= 1
      for j = length(f1):-1:1
         brack=(strcmp(elems,'(')|strcmp(elems,'«'))-(strcmp(elems,')')|strcmp(elems,'»'));
         fstart = popbrack(brack,f1(1) - 1, -1);
         elems{f1(j)} = ')';
         elems =  [elems(1:fstart-1) {'transpose','('} elems(fstart:end) ];
         %           if isvar(f1(j) - 1)
         %             elems{f1(j) - 1} = ['transpose(' elems{f1(j) - 1} ];
         %             elems{f1(j)} = ')';
         %          else
         %             error('grind:csolvers','(xx)'' is not supported, please replace with transpose((xx))');
         %          end
      end
      f=strcmp(elems,'«')|strcmp(elems,'»');
      elems = elems(~f);
      smodel = sprintf('%s', elems{:});
      obj = parsed_equation(smodel);
      elems = addbrackets(obj, true);
      %  [elems] = parseeq(smodel,1);
      feq=strcmp(elems, '=');
   end
end

matmults = find(strcmp(elems, '*')); %find the matrix multiplications and make auxiliary vars of the argument

for j = 1:length(matmults)
   if ~((j > 0)&&isnan(str2double(elems{matmults(j) - 1})))||~((j < length(elems))&&isnan(str2double(elems{matmults(j)+1})))
      elems{matmults(j)} = '.*';
   end
end
dim2 = 2;
if any(feq)
   isvar = false(size(elems));
   for j = 1:length(elems)
      isvar(j)=(isletter(elems{j}(1))||(elems{j}(1)=='_'));
   end
   vars = elems(isvar);
   no = i_getno(vars{1});
   if no.isvar
      dim2 = g_grind.statevars.dims{no.vecno}.dim2;
   end
end
if dim2 == 1
   mtimesfun = 'mtimes1d';
else
   mtimesfun = 'mtimes';
end
[elems] = changeop(elems, {'*'},mtimesfun);
%matrix funcs, the arguments are matrices, if an argument is a matrix
%expression, the arguments should become auxiliary vars. We need
%to add g_ii before the arguments.
matrixfuncs={'mtimes','mtimes1d','transpose','leftcells','rightcells','upcells', 'downcells','numel',...
   'size','length','sum','cumsum','prod','cumprod','flipud','fliplr','sumneighbors','repmat','g_min','min','max','mean','boxcarinflow','parlookup'};
%this function has matrix as argument, but we have to add it to an
%auxiliary var (if in expression)
%we only need to add g_ii here
matrixresfunc={'eye','linspace','logspace'};


%smodel = sprintf('%s', elems{:});
%[elems,  isvar] = parseeq(smodel);
isvar = false(size(elems));
for i = 1:length(elems)
   isvar(i)=isletter(elems{i}(1))||(elems{i}(1)=='_');
end
k = 1;
while k < length(elems)
   if isvar(k)&&strcmp(elems{k+1}, '(')
      f = strcmp(elems{k}, matrixfuncs);
      f1 = strcmp(elems{k}, matrixresfunc);
      if any(f)||any(f1)
         br1 = k + 1;
         arg = {'(','g_ii',','};
         %      br2=findvarfrom(elems,isvar,br1,1);
         if any(f)
            %check whether the arguments of the matrixfuncs are equations or
            %single matrices (equations need to be replaced
            %by generated auxiliary vars)
            brack=(strcmp(elems,'(')|strcmp(elems,'«'))-(strcmp(elems,')')|strcmp(elems,'»'));
            ncom = ',';
            fcomma1 = br1 + 1;
            while strcmp(ncom,',')
               fcomma = popbrack(brack, fcomma1, 1);
               ncom = elems{fcomma + 1};
               arg1 = elems(fcomma1:fcomma);
               if length(arg1) > 1
                  nx = length(g_grind.xfuncs) + 1;
                  f=strcmp(arg1,'«')|strcmp(arg1,'»');
                  g_grind.xfuncnames{nx} = sprintf('g_aux%03d', nx);
                  g_grind.xfuncs{nx}=sprintf('g_aux%03d=%s;', nx,sprintf('%s',arg1{~f}));
                  arg1 = g_grind.xfuncnames(nx);
               end
               arg = [arg arg1 {ncom}]; %sprintf('%s%s%s', arg, elems{fcomma1}, ncom);
               fcomma1 = fcomma + 2;
            end
            br2 = fcomma + 1;
            for i = 1:length(arg)
               if isletter(arg{i}(1))||(arg{i}(1)=='_')
                  arg{i} = sprintf('@%s', arg{i}); %flag to avoid adding (g_ii)
               end
            end
            elems = [elems(1:br1 - 1) arg elems(br2 + 1:end)];
         else
            elems=[elems(1:br1) {'g_ii',','}, elems(br1+1:end)];
         end
         isvar = false(size(elems));
         for i = 1:length(elems)
            isvar(i)=isletter(elems{i}(1))||(elems{i}(1)=='_');
         end
      end
   end
   k = k + 1;
end
%this is because a name conflicts with c++ functions
elems=cellrep(elems,'[]','g_emptyvar');
elems=cellrep(elems,'sum','g_sum');
elems=cellrep(elems,'min','g_min');  % in non vector models min is used to compare doubles
elems=cellrep(elems,'max','g_max');  % g_min and g_max are for matrices
fbox = strcmp(elems, 'boxcartrain');
if any(fbox)
   % boxcar_result   juv_=boxcartrain(1,juv,devjuv(0),sd(0));
   fbrack1 = find(strcmp(elems, '('), 1);
   fbrack2 = find(strcmp(elems, ')'));
   fbrack2 = fbrack2(end);
   fcomma = [find(strcmp(elems,',')), fbrack2];
   boxvar = elems(fbrack1 + 1:fcomma(1) - 1);
   resvar = elems{1};
   devvar = elems(fcomma(1) + 1:fcomma(2) - 1);
   if any(strcmp(devvar{1}, '«'))
      devvar = devvar(2:end - 1);
   end
   devvar = makecbrack(devvar, 1, '0');
   if length(fcomma) > 2
      cvvar = elems(fcomma(2) + 1:fcomma(3) - 1);
      if any(strcmp(cvvar{1}, '«'))
         cvvar = cvvar(2:end - 1);
      end
      cvvar = makecbrack(cvvar, 1, '0');
   else
      cvvar = {};
   end
   boxnr = find(strcmp(boxvar, g_grind.boxcar.names)) - 1; % - 1 for zero based c +  +
   elems=[{'boxcar_result',' ',resvar,'=','boxcartrain','(',num2str(boxnr),','},boxvar,{','},devvar,{','},cvvar,{')'}];
elseif any(feq)
   isvar = false(size(elems));
   for j = 1:length(elems)
      isvar(j)=isletter(elems{j}(1))||(elems{j}(1)=='_')||(elems{j}(1)=='.');
   end
   vars = elems(isvar);
   pvar = find(isvar);
   varno = [];
   for i1 = length(vars):-1:1
      no = i_getno(vars{i1});
      if i1 == 1
         varno = no;
      end
      isxfun = any(strcmp(vars{i1}, g_grind.xfuncnames))||~isempty(strfind(vars{i1},'.flow'));
      if isxfun||no.ispar||no.isvar||no.isfun||no.isperm
         ndx = pvar(i1);
         if ~isempty(ndx)
            elems = makecbrack(elems, ndx, 'g_ii');
         end
      end
   end
   if ~isempty(varno)&&varno.isvar
      f= find(strcmp(elems, '='), 1);
      for n = 1:f(1) - 1
         elems{n} = '';
      end
      elems{1} = sprintf('g_model.v_dydt[%d](g_ii)', varno.vecno - 1);
   end
   if ~strcmp(elems{end}, '; ')
      elems{end + 1} = ';';
   end
   elems=[{sprintf('for (int g_ii=0; g_ii<%s.numel; g_ii++){',vars{1})},elems,{'}'}];
end
fat = find(strncmp(elems, '@', 1)); %flag to avoid adding (g_ii)
for i = 1:length(fat)
   elems{fat(i)} = elems{fat(i)}(2:end);
end
f=strcmp(elems,'«')|strcmp(elems,'»');
elems = elems(~f);
smodel = sprintf('%s', elems{:});

function j = popbrack(brack, i, ddir)
% use parseeq to add brackets to elems
%create brack as:
%brack=(strcmp(elems,'(')|strcmp(elems,'«'))-(strcmp(fs,')')|strcmp(elems,'»'));
if ddir < 0
   brack = -brack;
end
if brack(i) == 0
   j = i;
else
   b = brack(i);
   j = i;
   while b > 0&&j > 0&&j<=length(brack)
      j = j + ddir;
      if j > 0
         b = b + brack(j);
      end
   end
end

function fs = makecbrack(fs, i, s)
if i > length(fs)
   i = length(fs);
end
if (i==length(fs))||~strcmp(fs{i+1}, '(')
   fs=[fs(1:i) {'(',s,')'} fs(i+1:end)];
else
   brack=(strcmp(fs,'(')|strcmp(fs,'«'))-(strcmp(fs,')')|strcmp(fs,'»'));
   j = popbrack(brack, i + 1, 1);
   fs=[fs(1:j-1) {'-','1'} fs(j:end)];
end

% function f = findvarfrom(s,istart, istep, precedence)
% %s and ivar should be created with parseq
% %use istart=find(fpos>=istart,1) to create an istart from a position in the
% %original string
% if nargin < 4
%    precedence = {''};
% end
%
% if length(s) < 2
%    f = istart;
%    return;
% end
% isvar = false(size(s));
% for i = 1:length(s)
%    isvar(i)=isletter(elems{i}(1))||(elems{i}(1)=='_');
% end
% f = istart;
% found = 0;
% brack = 0;
% if istep > 0
%    brack1 = '(';
%    brack2 = ')';
% else
%    brack1 = ')';
%    brack2 = '(';
% end
% while (f > 0)&&(f<=length(s))&&strcmp(s{f}, ' ')
%    f = f + istep;
% end
% while ~found
%    while ((f > 0)&&(f<=length(s)))&&(isvar(f)||~isnan(str2double(s{f}))||any(strcmp(precedence, s{f})))
%       f = f + istep;
%    end
%    if (f<=0)||(f > length(s))
%       found = true;
%    else
%       s1 = s(f);
%       prec = false;
%       if strcmp(s1, brack1)
%          brack = brack + 1;
%       elseif strcmp(s1, brack2)
%          brack = brack - 1;
%          if brack == 0 &&~prec
%             prec = f < length(s)&&any(strcmp(precedence, s{f+1}));
%             if ~prec
%                f = f + istep;
%             end
%          end
%       end
%       if brack <= 0 &&~prec
%          found = true;
%          f = f - istep;
%       else
%          f = f + istep;
%       end
%    end
% end
% if f == 0
%    f = 1;
% end
% if f > length(s)
%    f = length(s);
% end
% %if (istep==-1)&&(f > 1)&&(s(f)=='(')&&isvarletter(s(f-1)) %exp(x)^exp(x) left side solved
% if (istep==-1)&&(f > 1)&&(strcmp(s{f}, '(')&&isvar(f-1)) %exp(x)^exp(x) left side solved
%    f = f - 1;
%    while (f > 0)&&isvar(f)
%       f = f - 1;
%    end
%    f = f + 1;
% end
%
%
% faccent=find(strcmp('''',fs));
% %check whether there are strings in the equation, a single ' is interpreted as transpose,
% %a string cannot occur before =
% if ~isempty(faccent)
%    feq=find(strcmp('=', fs), 1);
%    if ~isempty(feq)&&(faccent(1) < feq)
%       faccent(1) = [];
%    end
%    nfaccent = length(faccent);
%    if nfaccent > 1 && (mod(nfaccent, 2)==0)%iseven
%       for i = 1:2:nfaccent
%          fs{faccent(i)} = sprintf('%s', fs{faccent(i):faccent(i + 1)});
%          for k = faccent(i + 1):-1:faccent(i) + 1
%             fs{k} = '';
%          end
%       end
%    end
%    fs = fs(~strcmp(fs, ''));
% end

function C = mtimes1d(A, B)  %#ok<DEFNU>
C = A * B;

function g_ss = sizefuncs(g_funcs, g_funcnames)
[N0] = i_initvar(1);
evalin('base',sprintf('%s\n',g_funcs{:}));

g_ss = cell(size(g_funcnames));
for g_ii1 = 1:length(g_funcnames)
   if isempty(strfind(g_funcnames{g_ii1},','))
      g_ss{g_ii1} = size(evalin('base', g_funcnames{g_ii1}));
   else
      g_ss{g_ii1} = [0, 0];
   end
end
i_keep(N0);






