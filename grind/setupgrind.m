%SETUPGRIND   Setup GRIND
%   Setup script for GRIND, adds the grind system directory to
%   the default path. If the grind directories are not found they are downloaded from
%   the download websites: http://content.alterra.wur.nl/webdocs/internet/aew/downloads/ or
%   http://www.sparcs-center.org/uploads/files/. After successful installation, setupgrind needs not to be run 
%   again.
%
%Usage:
%SETUPGRIND - standard setup
%   SETUPGRIND -d - runs setupgrind also if grind was already installed.
%   SETUPGRIND -s - creates a startup.m file that runs setupgrind each time that MATLAB
%      is opened. Used in older MATLAB versions if the path is not saved correctly.
%   SETUPGRIND -u - uninstalls grind (without removing files)
%   SETUPGRIND -up - update grind from the website (this is the same as calling 
%   <a href="matlab:help updategrind">updategrind</a>).
%   SETUPGRIND -restore - if grind is installed in an old version of MATLAB (< 7.12) some syntax needs to be changed.
%      'SETUPGRIND -restore' restores the original files.
%
%
%
%   See also Installing, updategrind
%
%   Reference page in Help browser:
%      <a href="matlab:commands('setupgrind')">commands setupgrind</a>

%   Copyright 2019 WUR
%   Revision: 1.2.1 $ $Date: 15-Jul-2019 21:00:41 $
function setupgrind(option)
%urls for downloading
if ~isoctave&&verLessThan('matlab','7.6')
    error('grind:setup','MATLAB version R2008a or newer is required for GRIND');
end
if nargin > 0&&strcmp(option, '-restore')
    checkversions(true);
    return;
end
urls = {'http://content.alterra.wur.nl/webdocs/internet/aew/downloads/','http://www.sparcs-center.org/uploads/files/'};
url = urls{1};
for i = 1:length(urls)
    %find first valid url
    [dummy,status]=urlread(sprintf('%s%s',urls{i},'grindversion.txt')); %#ok<ASGLU>
    if status == 1
        url = urls{i};
        break;
    end
end
doanyway = 0;
uninstall = 0;
addstartup = 0;
warning('off','backtrace');
if nargin == 1
    if strncmpi(option, '-d', 2)
        doanyway = 1;
    end
    if strncmpi(option, '-up', 3)
        updategrindfromweb(url);
        if isempty(which('grind_coco'))
            addpath(fullfile(grindpath,'sys2'));
        end
        %update coco from the grind website (only if coco is available)
        if ~isempty(grind_coco.findupdate(false,url))
            grind_coco.findupdate(true,url);
        end
        if ~isempty(grind_matcont.findupdate(false,url,true))
            grind_matcont.findupdate(true,url,true);
        end
        if ~isempty(grind_matcont.findupdate(false,url,false))
            grind_matcont.findupdate(true,url,false);
        end
    elseif strncmpi(option, '-u', 2)
        uninstall = 1;
    end
    if strncmpi(option, '-s', 2)
        addstartup = 1;
    end
end
k = which('grind.m');
if ~isempty(k) && ~doanyway
    thegrindpath = fileparts(k);
    if uninstall
        rmpath(thegrindpath);
        savepath; %adds permanently
        fprintf('Grind uninstalled from %s (no files removed)\n', thegrindpath);
    else
        checkversions(false); %if version < 7.12 [~, ] is not supported (replaced with [DUMM]
        if exist_newversion(url) < 1
            warning('GRIND:setupgrind:newerversion','There is a newer version available, use <a href="matlab: updategrind">updategrind</a> to update.');
        end
        fprintf('Grind was already installed in %s\n', thegrindpath);
        if addstartup
            doaddstartup;
        end
    end
    return;
elseif isempty(k)&&exist(fullfile('.','sys','sys2'),'dir')~=7
    %install grind using only the setupgrind.m file
    answer=inputdlg({'Grind needs to be downloaded, enter the base folder:'},...
        'Setup GRIND for MATLAB',1,{'d:\matlab\'});
    folder=answer{1};
    warning off MATLAB:MKDIR:DirectoryExists
    mkdir(folder);
    warning on MATLAB:MKDIR:DirectoryExists
    cd(folder);
    fprintf('Downloading GRIND...\n')
    unzip([url 'grind.zip']);
    disp('Successfully downloaded');
    cd grind
    setupgrind
    return
elseif uninstall
    fprintf('Grind was not installed\n');
    return;
end
try
    try
        mroot = OCTAVE_PATH;
    catch
        mroot = matlabroot;
    end
catch
    mroot = '';
end
thegrindpath=fullfile(mroot,'work','grind','sys');
k = which(fullfile(thegrindpath, 'grind.m'));
if isempty(k)
    thegrindpath = fullfile(pwd, 'sys');
end
addpath(thegrindpath);
try
    savepath; %adds permanently
catch
    path2rc;
end
checkversions(false); %if version < 7.12 [~, ] is not supported (replaced with [DUMM]
w = 0;
if ~canwriteto(grindpath)
    fileattrib(grindpath,'+w');
    if ~canwriteto(grindpath)
        warning('GRIND:setupgrind:writepermission','No permission to write to %s %s %s\n\n',grindpath,...
            'If you have permission in another directory it is possible to run GRIND.', ...
            'Preferrably set write permission to the full GRIND directory.');
        fileattrib(grindpath,'+w');
        w = 1;
    end
end
p=fullfile(mroot,'toolbox','local','pathdef.m');
if exist(p,'file') %if in future MATLAB versions pathdef.m does not exist there
    [~,aa]=fileattrib(p);
    if ~aa.UserWrite
        fileattrib(p,'+w');
        [~,aa]=fileattrib(p);
        if ~aa.UserWrite
            warning('GRIND:setupgrind:writepermission','The path cannot be stored as there is no write permission in %s. %s\n\n',...
                p, 'Each session GRIND needs to be reinstalled (or the permission should be changed).');
            w = 1;
        end
    end
end
if addstartup
    doaddstartup;
end
%if exist_newversion(url)<1
%    warning('GRIND:setupgrind:newerversion','There is a newer version available, use <a href="matlab: updategrind">updategrind</a> to update.');
%end
if w
    disp('GRIND installed, but there are some problems (see above)');
else
    disp('Grind installed successfully');
end

function doaddstartup
thegrindpath = grindpath;
u = userpath;
if length(u) < 4  %R11 gives not the workdirectory
    u = fullfile(matlabroot, 'work');
end
f = strfind(u, '; ');
if ~isempty(f)
    u = u(1:f(1) - 1);
end
u = fullfile(u, 'startup.m');
if exist(u, 'file') == 2
    %   lines = {};
    fid = fopen(u, 'r');
    hassetup = 0;
    while ~feof(fid)
        line = fgetl(fid);
        %      lines = {lines, line};
        if ~isempty(strfind(line, 'setupgrind'))  %#ok<STREMP>
            hassetup = 1;
        end
    end
    fclose(fid);
    if ~hassetup
        fid = fopen(u, 'a');
        fprintf(fid,'\ncd ''%s''\nsetupgrind\n',thegrindpath(1:end-4));
        fclose(fid);
        disp('setupgrind added to startup.m');
    else
        disp('setupgrind was already in startup.m');
    end
else
    fid = fopen(u, 'w');
    fprintf(fid,'\ncd ''%s''\nsetupgrind\n',thegrindpath(1:end-4));
    disp('startup.m created');
    fclose(fid);
end
function ok = canwriteto(path)
ok = 0;
oldpath = pwd;
try
    cd(path);
    fid=fopen('t.tmp','wt');
    ok =fid >= 0;
    if ok
        fclose(fid);
    end
    delete('t.tmp');
    cd(oldpath);
catch
    cd(oldpath);
end

function res=isoctave
res=exist('OCTAVE_VERSION', 'builtin') ~= 0;



function updategrindfromweb(url)

% if ~exist('urlread','file') %exists in R2008
%    error('GRIND:updategrind:MATLABversion','This function works only for newer versions of MATLAB');
% end
[nonewer, web, current] = exist_newversion(url);
if isempty(web)
    error('GRIND:updategrind:noconnection','Cannot download GRIND, is internet connected?')
end
if nonewer
    warning('GRIND:updategrind:nonewversion','There is no newer version of GRIND available');
else
    if exist('i_use', 'file')==2
        model('-c','1');
    end
    fprintf('Removing current version (%s)...\n', current.date)
    gp = fileparts(which('grind.m'));
    if ~isempty(gp)
        cd(gp);
    else
        gp = pwd;
    end
    delete('*.m');
    delete('*.mat');
    delete('*.exe');
    delete('*.chm');
    delete('*.chw');
    delete('*.h');
    cd('sys2')
    delete('*.m');
    cd ..
    if exist('csolver','dir')==7
        cd('csolver')
        delete('*.cpp');
        delete('*.h');
        copyfile('compile.bat','compile1.bat')
        cd ..
    end
    cd ..
    cd ..
    fprintf('Downloading new version (%s)...\n', web.date)
    unzip([url 'grind.zip']);
    cd(gp)
    cd('csolver')
    if exist('compile1.bat','file')==2
        copyfile('compile1.bat','compile.bat');
        delete('compile1.bat');
    end
    cd ..
    cd ..
    cd ..
    disp('Successfully updated');
end
fprintf('GRIND version %s\nRelease date: %s\n', web.version, web.date);

function [existnew, web, current] = exist_newversion(url)
existnew = 1;
current = [];
web = [];
try
    [s, status] = urlread([url 'grindversion.txt']);
    if status == 1
        h = str2cell(s);
    else
        return;
    end
catch
    return;
end
web.version = h{1};
web.date = h{2};
fid=fopen('use.m','r');
while ~feof(fid)
    line = fgetl(fid);
    f = strfind(line, 'Revision:');
    if ~isempty(f)
        f1 = strfind(line, '$');
        current.version = strtrim(line(f + 9:f1(1) - 1));
        f = strfind(line, 'Date:');
        current.date = strtrim(line(f + 5:f1(end) - 1));
    end
end
fclose(fid);
try
    existnew = datenum(web.date)-datenum(current.date) < 1E-4;
catch
    return
end

function checkversions(back)
if back||verLessThan('MATLAB','7.12')
    olddir = cd;
    %replace [~] with [DUMM1]
    if isempty(which('fileregexp'))
        addpath(fullfile(grindpath,'sys2'));
    end
    files={fullfile(grindpath,'*.m'),fullfile(grindpath,'sys2','*.m')};
    oldpath=cd;
    if ~back
        %change to dumm1
        nchanged=fileregexp(files,...
            struct('subdir',false,'action','replace', 'addcounter',true, 'noquote',true),'[~]+(?=[ ]*[,\]\)])','DUMM');
        cocopath=grind_coco.findupdate;
        cd(cocopath)
        nchanged=nchanged+fileregexp({'*.m'},...
            struct('subdir',true,'action','replace', 'addcounter',true, 'noquote',true),'[~]+(?=[ ]*[,\]\)])','DUMM');
    else
        %change dumm1 back to ~
        nchanged=fileregexp(files,...
            struct('subdir',false,'action','replace', 'addcounter',false, 'noquote',true),'DUMM[0-9]*(?![a-zA-Z_0-9])','~');
        cocopath=grind_coco.findupdate;
        cd(cocopath)
        nchanged=nchanged+fileregexp({'*.m'},...
            struct('subdir',true,'action','replace', 'addcounter',false, 'noquote',true),'DUMM[0-9]*(?![a-zA-Z_0-9])','~');
    end
    cd(oldpath);
    if nchanged > 0
        fprintf('Made %d changes in the code for this matlab version\n', nchanged)
    end
    cd(olddir);
end

