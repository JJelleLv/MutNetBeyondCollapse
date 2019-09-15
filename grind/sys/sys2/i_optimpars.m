function J = i_optimpars(x)
global t g_data g_t ;
% set parameters
if g_data.options.PosOnly && max(x < 0)
    J = 9999 + rand;
    return;
end

for i = 1:size(g_data.pars, 2)
    multassignin('base', g_data.pars{i}, x(i));
end

ndat = size(g_data.t, 1);
i_ru(t, g_data.t(ndat,1)-g_data.t(1,1), i_initvar,0);
YY=outfun(g_data.varlist);
g_data.pred=interp1(g_t,YY,g_data.t);

%penalty function: sumofsquares
if isa(g_data.options.penaltyfun,'function_handle')
    J=feval(g_data.options.penaltyfun,g_data);
else
    switch g_data.options.penaltyfun
        case 'normalized SS'
            J = 0;
            for k = 1:numel(g_data.varlist)
                if ~isnan(g_data.minobs(k))
                    %standardize different state variables between 0 and 1
                    pred = (g_data.pred(:, k) - g_data.minobs(k)) / (g_data.maxobs(k) - g_data.minobs(k));
                    obs = (g_data.obs(:, k) - g_data.minobs(k)) / (g_data.maxobs(k) - g_data.minobs(k));
                    for i = 1:ndat
                        if ~isnan(obs(i))
                            %sum of squares
                            J =  J + (pred(i) - obs(i))^2;
                        end
                        
                    end
                end
                
            end
            
        case 'SS'
            J = 0;
            for k = 1:size(g_data.obs,2)
                if ~isnan(g_data.minobs(k))
                    %or without standardization:
                    pred=g_data.pred(:,k);
                    obs=g_data.obs(:,k);
                    for i = 1:ndat
                        if ~isnan(obs(i))
                            %sum of squares
                            J =  J + (pred(i) - obs(i))^2;
                        end
                        
                    end
                end
                
            end
            
    end
end
g_data.iter = g_data.iter + 1;
if g_data.stopped > 0
    if g_data.stopped == 1
        error('GRIND:optimpars:Stopped','Stopped');
    else
        error('GRIND:optimpars:Cancelled','Cancelled, original parameter values are reset');
    end
    
end

if isempty(g_data.fval) || (g_data.fval > J)
    g_data.fval = J;
    g_data.X = x;
end

f=findobj(gcf,'Tag','ProgressEdit');
if ~isempty(f)
    v = [{sprintf('Param\tBest \tCurrent')}, g_data.pars];
    for i = 1:size(g_data.pars, 2)
        v{i + 1}=sprintf('%s\t%0.5g\t%0.5g',v{i + 1},g_data.X(i), x(i));
    end
    
    v=[v, {'',sprintf('Method: %s',g_data.optmethod) ,sprintf('Penalty function: %s',char(g_data.options.penaltyfun)), ...
        sprintf('Iteration   \t%0d',g_data.iter),sprintf('Current fit   \t%0.5g',J),...
        sprintf('Best fit    \t%0.5g',g_data.fval)}];
    set(f, 'String', v);
end

drawnow;

function multassignin(ws, name, V)
fbrack = strfind(name, '(');
if isempty(fbrack)
    assignin(ws, name, V);
else
    temp = evalin(ws, name(1: fbrack(1) - 1));
    s = ['temp' name(fbrack(1):length(name)) ' = ' num2str(V) ';' ];
    eval(s);
    assignin(ws, name(1, fbrack(1) - 1), temp);
end



