%E2N   Erase and redraw 2D Nullclines
%   Combination of the commands <a href="matlab:help era">era</a> and <a href="matlab:help null">null</a>
%
%   See also e2p, e2r, e3r
%
%   Reference page in Help browser:
%      <a href="matlab:commands('e2n')">commands e2n</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
%start
function e2n(varargin)
i_parseargs('','','',varargin);
era();
null();

