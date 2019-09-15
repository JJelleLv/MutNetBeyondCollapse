%VAL   initial and final values
%   Show initial and final values of state variables. Also to be used in a model
%   if you need to know the initial value of a state variable.
%  
%   Usage:
%   VAL  - displays the initial and final values of a state variable (after a run).
%   VAL('AVAR') - gets the initial value(s) of AVAR during a simulation.
%   VAL('AVAR(i,j)') - gets the initial value of AVAR(i,j) during a simulation.
%   VAL('argname',argvalue,...) - Valid argument <a href="matlab:commands func_args">name-value pairs</a> [with type]:
%     'output' [string] - output "full" or "partial" or name of state variable to output the initial value
%     'title' [string] - title for the input box (option -?)
%     'var' [state variable] - state variable
%   VAL('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-?' - edit the values of state variables using a dialog box
%     '-vector' - outputs the whole vector of initial conditions   
%   See also par, stat, model    
%
%   Reference page in Help browser:
%      <a href="matlab:commands('val')">commands val</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function v = val(varargin)
%(full, atitle)
global g_grind g_Y;

%use in a model to get the initial value (bypass i_parseargs for speed)
if nargout>0&&nargin==1&&ischar(varargin{1})
    if strcmp(varargin{1},'-vector')
       v=i_initvar;
    else
      v = evalin('base', varargin{1});
    end
    return;
end
fieldnams={'var', 'v', 'state variable';...
   'output', 's', 'output "full" or "partial" or name of state variable to output the initial value';...
   'title', 's', 'title for the input box (option -?)'}';
args=i_parseargs(fieldnams,'if(hasoption(''-?'')),deffields=''title'';else,deffields=''output,title'';end','-?,-vector',varargin);

if ~isfield(args,'output')
    args.output = '';
elseif any(strcmp(args.output,g_grind.statevars.names))||any(strcmp(args.output,g_grind.statevars.vectnames))
    v = evalin('base', args.output);
    return;
end
if isfield(args, 'var')
    v=evalin('base',args.var);
    return;
end
if ~isfield(g_grind, 'statevars')|| ~isfield(g_grind.statevars, 'dim')||(g_grind.statevars.dim==0)
    disp('No state variables defined');
    v=[];
else
    n = g_grind.statevars.dim;
    if n > 50&&~strcmp(args.output, 'full')
        n = 50;
        shortened = 1;
    else
        shortened = 0;
    end
    if any(strcmp(args.opts,'-?'))&&(n < 20)
        vals = cell(n, 1);
        N0 = i_initvar;
        for i = 1:n
            vals{i} = num2str(N0(i));
        end
        if ~isfield(args,'title')
            args.title = 'Set initial values';
        end
        answer = inputdlg(g_grind.statevars.names, args.title, 1, vals);
        if ~isempty(answer)
            for i = 1:n
                N0(i) = str2double(answer{i});
            end
            i_keep(N0);
        end
    end
    if nargout>0
        if ~g_grind.statevars.vector       
            N0=i_initvar;
            v=cell(1,g_grind.statevars.dim);
            for i=1:g_grind.statevars.dim
                v{i}=sprintf('%s=%s;',g_grind.statevars.names{i},num2str_complex(N0(i)));
            end
        else
            v=cell(1,length(g_grind.statevars.vectnames));
            for i=1:length(g_grind.statevars.vectnames)
                mat=evalin('base',g_grind.statevars.vectnames{i});
                if min(mat(:))==max(mat(:))
                    v{i}=sprintf('%s=zeros(%d,%d)+%.15g;',g_grind.statevars.vectnames{i},size(mat,1),size(mat,2),mean(mat(:)));
                else
                    v{i}=sprintf('%s=%s;',g_grind.statevars.vectnames{i},mat2str(mat));
                end
            end
        end
    elseif isempty(g_Y)
        disp('Initial values of state variables:');
        N0 = i_initvar;
        maxparlen = par('-maxparlen');
        s=['%-' num2str(maxparlen) 's = %s\n'];
        s1=['%-' num2str(maxparlen) 's = %s      %%%s,%s\n'];
        for i = 1:n
            avar = i_statevars_names(i);
            [adescr,aunit] = par('-d', avar);
            if ~isempty(aunit)||~isempty(adescr)
                fprintf(s1, i_statevars_names(i), num2str_complex(N0(i),'%0.6g'), adescr, aunit);
            else
                fprintf(s, i_statevars_names(i), num2str_complex(N0(i),'%0.6g'));
            end
        end
    else
        disp('Initial and final values of state variables in last run:');
        maxparlen = par('-maxparlen');
        s=['%-' num2str(maxparlen) 's = %s / %s\n'];
        s1=['%-' num2str(maxparlen) 's = %s / %s      %%%s,%s\n'];
        statevars=i_statevars_names;
        for i = 1:n
            avar = statevars{i};
            [adescr,aunit] = par('-d', avar);
            if ~isempty(aunit)||~isempty(adescr)
                fprintf(s1,i_statevars_names(i),num2str_complex(g_Y(1, i),'%0.5g'),num2str_complex(g_Y(end,i),'%0.5g'),adescr,aunit);
            else
                fprintf(s,i_statevars_names(i),num2str_complex(g_Y(1, i),'%0.5g'),num2str_complex(g_Y(end,i),'%0.5g'));
            end
        end
    end
    if shortened
        fprintf('Only the first 50 of %d elements shown\n', g_grind.statevars.dim);
    end
end
if isfield(g_grind, 'permanent')&&~isempty(g_grind.permanent)
    disp(' ');
    defpermanent('-l');
end

function s=num2str_complex(N0,format)
if nargin<2
    format='%g';
end
if isreal(N0)
    s=sprintf(format,N0);
else
    formatc=sprintf('%s%%+%si',format, format(2:end));
    s=sprintf(formatc,real(N0),imag(N0));
end


