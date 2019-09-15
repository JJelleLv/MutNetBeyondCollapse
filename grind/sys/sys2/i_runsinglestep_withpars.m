function res = i_runsinglestep_withpars(t, N0, update,parno,P0)
%run a single step with the current odefile
global g_grind;
if isempty(parno)
    res = i_runsinglestep(t, N0, update);
else
    if update
        %update parameters!
        g_grind.hfun.curr = i_getodehandle(1,sprintf(',%s',g_grind.pars{parno}));
    end
    if isempty(N0)
        res=[];
    else
        siz1 = size(N0);
        if g_grind.statevars.dim<siz1(end)
            siz1(end+1)=g_grind.statevars.dim;
        end

        switch length(siz1)
            case 3
                N0 = reshape(N0, siz1(1) * siz1(2), siz1(3));
                P0 = reshape(P0, siz1(1) * siz1(2), length(parno));
            case 4
                N0 = reshape(N0, siz1(1) * siz1(2) * siz1(3), siz1(4));
                P0 = reshape(P0, siz1(1) * siz1(2), length(parno));
        end

        res = ones(size(N0));
        if length(t)==1
            t=t+zeros(size(N0,1),1);
        end

        if strcmp(g_grind.solver.opt.Vectorized,'on')
            if length(parno)==1
                yy = feval(g_grind.hfun.curr, transpose(t), transpose(N0),transpose(P0(:,1)));
            elseif length(parno)==2
                yy = feval(g_grind.hfun.curr, transpose(t), transpose(N0),transpose(P0(:,1)),transpose(P0(:,2)));
            end

            res = transpose(yy);
        else  
            if length(parno)==1
                for i = 1:size(N0, 1)
                    yy = feval(g_grind.hfun.curr, t(i), transpose(N0(i,:)),P0(i,1));
                    res(i, :) = transpose(yy);
                end

            elseif length(parno)==2
                for i = 1:size(N0, 1)
                    yy = feval(g_grind.hfun.curr, t(i), transpose(N0(i,:)),P0(i,1),P0(i,2));
                    res(i, :) = transpose(yy);
                end
            end

        end

        if length(siz1) > 2
            res = reshape(res, siz1);
        end

    end

end

