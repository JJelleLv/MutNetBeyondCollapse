%CONTBIF   Continue a bifurcation in two dimensions with COCO
%   This command is no longer necessary, it does the same thing as <a href="matlab:help conteq">conteq</a>, but it selects by
%   default the 2 parameters in the user interface.
%
%   Usage:
%   See <a href="matlab:help conteq">conteq</a>
%  
%
%   See also conteq
%
%   Reference page in Help browser:
%      <a href="matlab:commands('contbif')">commands contbif</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:39 $
function contbif(varargin)
if nargin>0
    conteq('-contbif',varargin{:});
else
    conteq('-contbif');
end

