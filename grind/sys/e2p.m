%E2P   Erase and redraw 2D Phase plane
%   Combination of the commands <a href="matlab:help era">era</a>, <a href="matlab:help null">null</a> and <a href="matlab:help ru">ru </a>
%
%   See also e2n, e2r, e3r
%
%   Reference page in Help browser:
%      <a href="matlab:commands('e2p')">commands e2p</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
%start
function e2p(varargin)
i_parseargs('','','',varargin);
era();
null();
ru();
