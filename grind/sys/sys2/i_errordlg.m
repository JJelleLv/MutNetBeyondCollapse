function i_errordlg(varargin)
global g_grind;
h=errordlg(varargin{:});
if isfield(g_grind,'figopts')&&~isempty(g_grind.figopts)
    set(h,g_grind.figopts{:});
end


