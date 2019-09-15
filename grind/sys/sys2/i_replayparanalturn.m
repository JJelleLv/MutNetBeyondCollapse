function [ok,buttontext]=i_replayparanalturn
global g_paranal;
ok=0;
buttontext='PLAY >';
if isfield(g_paranal,'prevrun')&&~isempty(g_paranal.prevrun)
    ok=1;
    %exchange Y and prevY etc.
    h1=g_paranal.prevrun;
    g_paranal.prevrun=g_paranal.run;
    g_paranal.run=h1;
    if g_paranal.run.parvalues(1)>g_paranal.run.parvalues(end)
         buttontext='Backward >';
    else
         buttontext='Forward >';      
    end
end

