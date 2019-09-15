%GCD1 - grind CD
%   Handy way of navigating in grind subdirectories. Search order:
%     1. Current directory
%     2. Grind directory and subdirectories
%     3. MATLAB directory and subdirectories
%
%   GCD1 filename - searches filename and goes to directory.
%   GCD1 pathname - goes to subdirectory
%
function gcd1(varargin)
% if nargin == 2
%    d=which('gcd','-all');
%    if length(d) > 1
%       curr = pwd;
%       cd(fileparts(d{length(d)}));
%       if nargout == 0
%          gcd(varargin{:})
%       elseif nargout == 1
%          Ans = gcd(varargin{:});
%          varargout = {Ans};
%       else
%          [a1, a2, a3] = gcd(varargin{:});
%          varargout = {a1, a2, a3};
%       end
%       cd(curr);
%       return;
%    end
% end
if isempty(varargin)
   [~, adir] = uigetfile('*.*', 'Change dir');
   cd(adir);
   return;
end
adir = strtrim(sprintf('%s ', varargin{:}));
afil = findgrindfile(adir);
if ~isempty(afil)
   cd(afil)
else
   error('GRIND:gcd:cannotfinddir','Cannot find directory/file %s',adir);
end

