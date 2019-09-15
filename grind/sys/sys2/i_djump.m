function res = i_djump(no,t,varargin)
global g_grind
args=g_grind.solver.djump.args(no);
if ~isfield(args,'timing_fun')
    djump('-i');
    args=g_grind.solver.djump.args(no);
end

if g_grind.statevars.vector
    siz=[g_grind.statevars.dims{args.var}.dim1,g_grind.statevars.dims{args.var}.dim2];
    res=zeros(siz);
else
    siz=[];
    res=0;
end
%synchronize the timing of different terms?
if args.synchronize>0&& isfield(g_grind.solver.djump,'nextt')&&length(g_grind.solver.djump.nextt)>=args.synchronize
    sync_nextt=g_grind.solver.djump.nextt{args.synchronize};
else
    %no synchronized value available
    sync_nextt=[];
end
if isempty(args.nextt)
    if ~isempty(sync_nextt)
        args.nextt=sync_nextt;
    else
        args.nextt=t+args.timing_fun(varargin{1:length(args.timing_pars)});
        if args.synchronize>0
            g_grind.solver.djump.nextt{args.synchronize}=args.nextt;
        end
    end
else
    while t>args.nextt
        if ~isempty(siz)
            if strcontains(char(args.size_fun),'g_size_var')
                res=res+args.size_fun(varargin{length(args.timing_pars)+1:end},siz)/solver('step');
            else
                for i=1:siz(1)
                    for j=1:siz(2)
                        res(i,j)=res(i,j)+args.size_fun(varargin{length(args.timing_pars)+1:end})/solver('step');
                    end
                end
            end
        else
            res=res+args.size_fun(varargin{length(args.timing_pars)+1:end})/solver('step');
        end
        if isempty(sync_nextt)||t>sync_nextt
            %if t is also larger than sync_nextt update sync_next
            args.nextt=args.nextt+args.timing_fun(varargin{1:length(args.timing_pars)});
            if args.synchronize>0
                g_grind.solver.djump.nextt{args.synchronize}=args.nextt;
                sync_nextt=args.nextt;
            end
        else
            args.nextt=sync_nextt;
        end
    end
end
g_grind.solver.djump.args(no)=args;
