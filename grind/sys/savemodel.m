%SAVEMODEL   Save current parameters and model
%    Save the model and the default parameter settings. After changing a model
%    with the command <a href="matlab:help model">model</a> this is done automatically.
%
%    Usage:
%    SAVEMODEL- prompts for a filename and saves the model to the 
%    that file.
%    SAVEMODEL FILENAME - saves the model as FILENAME.
%    SAVEMODEL FILENAME 1 - saves the model as FILENAME and overwrites
%    existing files without confirmation.
%    SAVEMODEL('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'file' [file name] - Filename for saving the model
%     'overwrite' [logical] - If true, the file is overwritten without notice
%     'type' [coco | ini or empty] - Type of m-file (coco|ini|empty)
%
% 
%    See also savepar, model
%
%   Reference page in Help browser:
%      <a href="matlab:commands('savemodel')">commands savemodel</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function savemodel(varargin)
%afilename, overwrite, filetype
global g_grind;
fieldnams={'file', 'F', 'Filename for saving the model','';...
   'overwrite', 'l', 'If true, the file is overwritten without notice',false;...
   'type', 'e[coco|ini]#E', 'Type of m-file (coco|ini|empty)',''}';
args=i_parseargs(fieldnams,'if nargs>1&&~argtype(2,''l''),deffields=''file,type'';else,deffields=''file,overwrite,type'';end', '',varargin);
if ~isfield(args,'type')
    args.type='';
end
if ~isfield(args,'overwrite')
    args.overwrite = false;
end

if ~isfield(g_grind, 'inifile')
    i_errordlg('No model selected');
    error('GRIND:savemodel:NoModel','No model to save');
end
if ~isfield(args,'file')
    if isempty(g_grind.inifile)
        args.file=uiputfile('*.ini', 'Save model as');
    else
        if isempty(fileparts(g_grind.inifile))
            args.file = uiputfile(fullfile(cd,g_grind.inifile), 'Save model as');
        end
    end
    args.overwrite = 1;
end

%if length(g_grind.model) == 0
%   errordlg('No model specified');
%   error('savemodel');
%end

if iscell(args.file)
    args.file = args.file{1};
end
if isempty(strfind(args.file, '.'))
    if strcmp(args.type,'coco')
        args.file = [args.file '.m'];
    else
        args.file = [args.file '.ini'];
    end
end
if exist(args.file, 'file') && ~args.overwrite
    but=questdlg(['Ok to overwrite ' ,args.file, '?'],'Saving model','Yes','No','Cancel','Yes');
else
    but = 'Yes';
end
LF = sprintf('\n'); %#ok<SPRINTFN>
if strcmp('Yes', but)
    [path, name, ext] = fileparts(args.file);
    if exist(args.file, 'file')
        ext = ['.~', ext(2:length(ext))];
        copyfile(args.file, fullfile(path,[name ext]));
    end
    if ~isempty(strfind(args.file, '.all'))
        saveall(args.file);
        return;
    end
    if ~isempty(strfind(args.file, '.m'))
        saveas_m(args.file, args.type);
        return;
    end
    fid = fopen(args.file, 'w');
    try
        fwrite(fid, ['%model', LF]);
        for i = 1:length(g_grind.model)
            if ~isempty(g_grind.model{i})
                fwrite(fid, [char(g_grind.model{i}) LF]);
            end
        end
        fwrite(fid, ['%commands' LF]);
        for i = 1:length(g_grind.commands)
            if ~isempty(g_grind.commands{i})
                fwrite(fid, [char(g_grind.commands{i}) LF]);
            end
        end
        if isfield(g_grind, 'scheme') && ~isempty(g_grind.scheme)
            fwrite(fid, ['%scheme' LF]);
            for i = 1:length(g_grind.scheme)
                if ~isempty(g_grind.scheme{i})
                    fwrite(fid, [char(g_grind.scheme{i}) LF]);
                end
            end
        end
        fclose(fid);
    catch err
        fclose(fid);
        rethrow(err);
    end
    disp(['Model saved to ', args.file]);
    g_grind.inifile = args.file;
elseif strcmp('No', but)
    args.file=uiputfile('*.ini', 'Save model as');
    savemodel(args.file, 1);
end

function saveall(afilename)
global g_grind;
%g_version = 3;
%who global
%eval(i_globalstr(who('global')));
save(afilename, g_grind, '-mat');
disp(['Model and current state saved to ', afilename]);


function saveas_m(afilename, filetype)
global g_grind t g_cont;
fid = fopen(afilename, 'w');
for i1 = 1:length(g_grind.model)
    f = strfind(g_grind.model{i1}, '%');
    if ~isempty(f)&&min(f)==1
        fprintf(fid, '%s\n', g_grind.model{i1});
    else
        fprintf(fid, '%%%s\n', g_grind.model{i1});
    end
end
fprintf(fid, '%%\n');
[~, name] = fileparts(afilename);
if ~strcmp(filetype,'coco')
    fprintf(fid, '%%Usage:\n');
    fprintf(fid,'%%[tvalues,yvalues] = %s(trange,nsteps, N0)\n',name);
    fprintf(fid, '%%trange = range for time\n');
    fprintf(fid, '%%nsteps = number of steps for output (NAN=maximum)\n');
    fprintf(fid, '%%N0 = vector with initial values\n');
    fprintf(fid,'function [tvalues,yvalues] = %s(trange,nsteps, N0)\n',name);
    fprintf(fid, 'if nargin < 1\n');
    fprintf(fid, '   trange = [%g %g];\n', t, g_grind.ndays);
    fprintf(fid, 'end\n');
    fprintf(fid, 'if nargin < 2\n');
    fprintf(fid, '   nsteps = %g;\n', g_grind.tstep);
    fprintf(fid, 'end\n');
    fprintf(fid, 'if nargin < 3\n');
    fprintf(fid, '   N0 = %s;\n', mat2str(i_initvar));
    fprintf(fid, 'end\n');
    fprintf(fid, 'if length(trange) == 1\n');
    fprintf(fid, '   trange = [0 trange];\n');
    fprintf(fid, 'end\n');
    fprintf(fid, 'if ~isnan(nsteps)\n');
    fprintf(fid,'   ts = (trange(1):(trange(2) - trange(1)) / nsteps:trange(1) + trange(2))'';\n');
    fprintf(fid, '   if length(ts) == 2\n');
    fprintf(fid, '      ts = [ts(1); ts(1) + (ts(2) - ts(1)) / 2; ts(2)];\n');
    fprintf(fid, '   end\n');
    fprintf(fid, 'else\n');
    fprintf(fid, '   ts = trange;\n');
    fprintf(fid, 'end\n');
    
    fprintf(fid, '\n%%run model:\n');
    fprintf(fid, 'opt.AbsTol=%g;\n', g_grind.solver.opt.AbsTol);
    if ~isempty(g_grind.solver.opt.StepSize)
        fprintf(fid, 'opt.StepSize=%g;\n', g_grind.solver.opt.StepSize);
    end
    fprintf(fid, 'opt.RelTol=%g;\n', g_grind.solver.opt.RelTol);
    fprintf(fid,'funhandles=get_model;\n');
    fprintf(fid,'[ts, Ys] = %s(funhandles.odehandle, ts, N0, opt);\n',g_grind.solver.name);
    fprintf(fid, 'if nargout == 2\n');
    fprintf(fid, '   tvalues = ts;\n');
    fprintf(fid, '   yvalues = Ys;\n');
    fprintf(fid, 'elseif nargout == 1\n');
    fprintf(fid,'   tvalues = [ts, Ys];\n');
    fprintf(fid, 'else\n');
    fprintf(fid, '   figure;\n');
    fprintf(fid,'   plot(ts,Ys);\n');
    fprintf(fid,'   xlabel(''%s'');\n',g_grind.outt{1});
    fprintf(fid,'   ylabel(''state'');\n');
    fprintf(fid, 'end\n\n');
    if any(strcmp(g_grind.solver.name,{'i_differ','euler'}))
        filecopy(fid, g_grind.solver.name);
    end
else
    fprintf(fid,'function [bd1] = %s\n',name);
    if ~isempty(g_cont)
        fprintf(fid, 'mod=get_model;\n%% coco settings\nprob = coco_prob();\n');
        fprintf(fid, 'prob = coco_set(prob, ''ode'', ''vectorized'',''%s'');\n',  g_grind.solver.opt.Vectorized);
        fprintf(fid, 'prob = coco_set(prob,''EP'', ''NSA'',true,''BTP'',true);\n');
        tbs=fieldnames(g_cont.settings.coco);
        settings=g_cont.eval_settings;
        for i=1:length(tbs)
            tb=tbs{i};
            f=fieldnames(settings.coco.(tb));
            fprintf(fid,'prob = coco_set(prob,''%s''',tb);
            for j=1:length(f)
                v=settings.coco.(tb).(f{j});
                if isempty(v)
                    fprintf(fid,',''%s'',[]',f{j});
                else
                    fprintf(fid,',''%s'',%g',f{j},v);
                end
            end
            fprintf(fid,');\n');
        end
        fprintf(fid, '\n%%run coco 1D\nprob = ode_isol2ep(prob, '''',mod.coco_ode, mod.coco_jac, mod.coco_jacp,%s, mod.parnames, mod.p0);\n',mat2str(i_initvar));
        fprintf(fid,'bd1 = coco(prob,''1'', [], 1,  ''%s'', %s);\n\n',g_cont.settings.freepars{1},mat2str(g_cont.settings.parranges{1}));
    end
end
% if isfield(g_grind, 'permanent')&&~isempty(g_grind.permanent)
%    fprintf(fid, 'global g_permanent;\n');
%    for i = 1:length(g_grind.permanent)
%       fprintf(fid,'g_permanent.permanent{%d}.name=''%s'';\n',i,g_grind.permanent{i}.name);
%       fprintf(fid, 'g_permanent.permanent{%d}.currvalue=%g;\n', i, g_grind.permanent{i}.currvalue);
%       fprintf(fid, 'g_permanent.permanent{%d}.initiate=%d;\n', i, g_grind.permanent{i}.initiate);
%    end
% end
fprintf(fid,'\nfunction g_l_res=get_model\n%%parameters:\n');
pp = par; %grind command par
fprintf(fid, '%s\n', pp{:});
fprintf(fid,'g_l_res.odehandle=%s;\n',func2str(i_getodehandle('normal')));
fprintf(fid,'g_l_res.coco_ode=%s;\n',func2str(i_getodehandle('coco')));
fprintf(fid,'g_l_res.coco_jac=%s;\n',func2str(i_getodehandle('coco_jac')));
fprintf(fid,'g_l_res.coco_jacp=%s;\n',func2str(i_getodehandle('coco_jacp')));
s=sprintf('%s;',g_grind.pars{:});
fprintf(fid,'g_l_res.p0=[%s];\n\n', s(1:end-1));
s=sprintf('''%s'',',g_grind.pars{:});
fprintf(fid,'g_l_res.parnames={%s};\n\n', s(1:end-1));
m1 = i_analysemodel(g_grind.model);
m1.theodefile=strrep(m1.theodefile,'=curr_ode0(',['=' g_grind.odefile '(']);
fprintf(fid, '%s\n', m1.theodefile{:});
fprintf(fid,'\nfunction m=i_reshapevcat(siz,varargin)\nif siz(3)>1\n    flen1=find(cellfun(''size'',varargin,2)==1);\n    if ~isempty(flen1)\n');
fprintf(fid,'for i=1:length(flen1)\n            varargin{flen1(i)}=varargin{flen1(i)}(ones(1,siz(3)));\n        end\n    end\nend\nm=reshape(vertcat(varargin{:}),siz);\n');

fclose(fid);
function filecopy(fid, name, rep1, to1)
fid2 = fopen(which(name), 'r');
if fid2 < 0
    error('grind:savemodel:error','Error opening file %s',name);
end
while ~feof(fid2)
    line = fgetl(fid2);
    if nargin == 4
        line = strrep(line, rep1, to1);
    end
    fprintf(fid, '%s\n', line);
end
