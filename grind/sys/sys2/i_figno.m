function H = i_figno(afig)
global g_grind;
if ~isfield(g_grind, 'timevars')
   ntime = 8;
else
   ntime = max(8, length(g_grind.timevars));
end

switch afig
 case 'time'
   H = 1;
 case 'phase2'
   H = 1 + ntime;
 case 'phase3'
   H = 2 + ntime;
 case 'phase1'
   H = 3 + ntime;
 case 'phase1a'
   H = 4 + ntime;
 case 'funplot'
   H = 5 + ntime;
 case 'dirfield'
   H = 6 + ntime;
 case 'eigen'
   H = 7 + ntime;
 case 'paranal'
   H = 11 + ntime;
 case 'poinsec'
   H = 14 + ntime;
 case 'potent1'
   H = 15 + ntime;
 case 'potent2'
   H = 16 + ntime;
 case 'potent3'
   H = 17 + ntime;
 case 'poinmap'
   H = 18 + ntime;
 case 'lyap1'
   H = 19 + ntime;
 case 'lyap2'
   H = 20 + ntime;
 case 'itermap'
   H = 21 + ntime;
 case 'takens'
   H = 22 + ntime;
 case 'lorenzmap'
   H = 23 + ntime;
 case 'combfig'
   H = 24 + ntime;
 case 'autocorr'
   H = 25 + ntime;
 case 'torus'
   H = 26 + ntime;
 case 'attrbasin'
   H = 27 + ntime;
 case 'trdet'
   H = 28 + ntime;
 case 'dialog'
   H = 29 + ntime;
 case 'setmat'
   H = 30 + ntime;
 case 'mcarlo'
   H = 31 + ntime; %4 figs
 case  {'conteq2d','contbif'}
   H = 35 + ntime;
 case  'conteq'
   H = 36 + ntime; 
case  'growths'
   % for each state variable one growth plot
   H = 37 + ntime;
 case 'obspred'
   H = 38 + ntime+ nstatevars;
 case 'vectplot'
   H = 39 + ntime + 2 * nstatevars;
 case 'varcontour' %obsolete, synonym of vectplot
   H = 39 + ntime + 2 * nstatevars;
 case 'uncertain'
   H = 39 + ntime + 3 * nstatevars;
 case 'paranal2d'
   H = 39 + ntime + 4 * nstatevars;
 case  'viewcells'
   H = 39 + ntime + 4 * nstatevars + 1 * nvectors;
 case 'fokkerplanck1'
   H = 41 + ntime + 4 * nstatevars + 2 * nvectors;
  case 'fokkerplanck2'
   H = 42 + ntime + 4 * nstatevars + 2 * nvectors;
  case 'fokkerplanck3'
   H = 43 + ntime + 4 * nstatevars + 2 * nvectors;
  case 'fokkerplanck4'
   H = 44 + ntime + 4 * nstatevars + 2 * nvectors;
 case 'maxno'
   H = 45 + ntime + 4 * nstatevars + 2 * nvectors;
 case 'quasipot'
   H = 46 + ntime + 4 * nstatevars + 2 * nvectors;
 case 'lyapspect'
   H = 47 + ntime + 4 * nstatevars + 2 * nvectors;
 otherwise
   error('GRIND:figno:UnknownType','Unknown figure: %s', afig);
end

function n = nstatevars
global g_grind
if ~isfield(g_grind, 'statevars')
   n = 0;
else
   n = g_grind.statevars.dim;
end

function n = nvectors
global g_grind;
if ~isfield(g_grind, 'statevars')
   n = 0;
elseif ~isempty(g_grind.statevars) && g_grind.statevars.vector
   n = length(g_grind.statevars.dims);
else
   n = 0;
end

