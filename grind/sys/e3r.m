%E3R   Erase and redraw 3D phase plane without nullclines
%   Combination of the commands <a href="matlab:help era">era</a>, <a href="matlab:help phas">phas 3</a> and <a href="matlab:help ru">ru</a>.
%
%   See also e2p, e2r, e2n
%
%   Reference page in Help browser:
%      <a href="matlab:commands('e3r')">commands e3r</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
%start
function e3r(varargin)
i_parseargs('','','',varargin);
era();
phas(3);
ru();
