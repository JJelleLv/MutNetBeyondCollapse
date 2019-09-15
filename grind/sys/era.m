%ERA   Erase and close all figures
%   Close all open GRIND figures.
%       
% ERA('-opt1','-opt2',...) - Valid command line <a href="matlab:commands func_args">options</a>:
%     '-a' - close all open figures.
%
%   See also e2n, e3r, e2r, e2p
%
%   Reference page in Help browser:
%      <a href="matlab:commands('era')">commands era</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:40 $
function era(varargin)
args=i_parseargs('','','-a',varargin);
allfigures=any(strcmp(args.opts, '-a'));
if allfigures
    close all hidden;
else
    ch = findobj(get(0,'Children'), 'flat','type','figure');
    max = i_figno('maxno');
    i2 = i_figno('combfig');
    i3=i_figno('dialog');
    for i = 1:length(ch)
        if ~isnumeric(ch(i))&&isprop(ch(i),'Number')
            no=ch(i).Number;
        else
            no=ch(i);
        end
        if ~isempty(no)&&isnumeric(no)&&(no <= max) && (no ~= i2)&& (no ~= i3)
            close(ch(i));
        end
    end
end
