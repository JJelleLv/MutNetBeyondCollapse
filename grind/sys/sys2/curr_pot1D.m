%function created by GRIND
function g_X2 = curr_pot1D(g_X1, ~)
global g_grind;
g_X2 = -g_grind.hfun.curr(1 + 1E-5, g_X1);

