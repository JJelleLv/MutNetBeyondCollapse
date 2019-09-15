%GRINDPATH   Full path of the grind/sys directory
%   Display the path of the GRIND program files
%
%   Reference page in Help browser:
%      <a href="matlab:commands('grindpath')">commands grindpath</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function p = grindpath(sys)
if nargin==0
   sys=1;
end
if sys==1
   p=fileparts(which('grind.m'));
elseif sys ==0
%  the root of the path of the grind ini files, you may change these lines
  root=fileparts(which('grind.m'));
  root=root(1:length(root)-4);
  grindroot=root(1:length(root)-6);
  if isoctave
     userp=[];
  else
    userp=regexp(userpath,';','split');
    userp=userp(~cellfun('isempty',userp));
  end
  p = [{root,grindroot},userp];
elseif sys==2
  root=fileparts(which('grind.m'));
  p=root(1:length(root)-4);
end
