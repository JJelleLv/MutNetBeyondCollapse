%E2R   Erase and redraw 2D phase plane without nullclines
%   Combination of the commands <a href="matlab:help era">era</a>, <a href="matlab:help phas">phas 2</a> and <a href="matlab:help ru">ru</a>.
%
%   See also e2p, e2n, e3r
%
%   Reference page in Help browser:
%      <a href="matlab:commands('e2r')">commands e2r</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
%start
function e2r(varargin)
i_parseargs('','','',varargin);
era();
phas(2);
ru();
