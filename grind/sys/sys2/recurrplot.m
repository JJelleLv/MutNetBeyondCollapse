function recurrplot(varargin)
global g_grind g_Y t
fieldnams={'npoints', '', ' ';...
   'fracpoints', 'n>=0&&n<=1', ' '}';
arg = i_time_options(i_parseargs(fieldnams,'fracpoints',...
    {'-r','-a','-n'},varargin));
if ~isfield(arg,'npoints')
    arg.npoints=g_grind.tstep;
end

if isnan(arg.npoints)
    arg.npoints=1000;
end

if ~isfield(arg,'fracpoints')
    arg.fracpoints=0.04;
end


oldtstep = g_grind.tstep;
if isnan(oldtstep)
    g_grind.tstep = args.ndays;
end

try
    N0=i_initvar;
    if arg.rerun
        i_ru(t, arg.npoints, N0, 0);
    end

    g_grind.tstep = oldtstep;
catch err
    %   err=lasterror;
    g_grind.tstep = oldtstep;
    rethrow(err);
end

sumdiff=zeros(size(g_Y,1));
for i=1:size(g_Y,2)
    [x,y]=meshgrid(g_Y(:,i));
    sumdiff=sumdiff+(x-y).^2;
end

sumdiff=sumdiff.^0.5;
figure
ncritpoints=round(numel(sumdiff)*arg.fracpoints);
critdiff=fzero(@(x)getnumpoints(x,sumdiff,ncritpoints),1);
pcolor(double(sumdiff<critdiff))
i_plotdefaults;
xlabel('Time (t)');
ylabel('Time (t)');

shading flat
colormap(flipud(gray))

function res=getnumpoints(critdiff,sumdiff,ncritpoints)
res=sum(sum(double(sumdiff<critdiff)))-ncritpoints;
