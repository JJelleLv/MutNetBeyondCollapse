%UPDATEGRIND   Download newest version from website and install this
% version, overwriting the current version. It can also update MATCONT and COCO
% if a new version is available
%
%  Usage: 
%  updategrind
%
%  See also setupgrind
%
%   Reference page in Help browser:
%      <a href="matlab:commands('updategrind')">commands updategrind</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function updategrind
curdir=pwd;
try
   cd(grindpath)
   cd ..
   setupgrind('-update');
   cd(curdir)
catch err
   cd(curdir)
   rethrow(err)
end