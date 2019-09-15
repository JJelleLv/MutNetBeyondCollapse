% FUNCTION FILEREGEXP Apply regexprep to files in the directory
%
% fileregexp(files,options,pattern,reptext) 
% files = string or cell list with filenames including wildcards
% options = structure with fields:
%    option.subdir = true/false search subdirectories
%    options.action = 'show' (show only) 'replace' (replace all) or 'query'
%    (show and ask for each line)
%     options.addcounter = true/false add counter to the replacement text
%     options.noquote = true/false - skip quotes or comments
% pattern,reptext = see REGEXPREP
% 
% for example change all [~,]=function a(~,a) to
% [dummy1,]=function a(dummy2,a)
%
%fileregexp('*.m',struct('subdir',true,'action','show', 'addcounter',true, 'noquote',true),'[~]+(?=[ ]*[,\]\)])','dummy')
% change back:
% fileregexp('*.m',struct('subdir',true,'action','show', 'addcounter',false, 'noquote',true),'dummy[0-9]*(?![a-zA-Z_0-9])','~')
%
%or temporary remove % in all files
%fileregexp('*.m',struct('subdir',true,'action','replace'),'[%][#]ok','%_#ok','ignorecase')
%fileregexp('*.m',struct('subdir',true,'action','replace'),'[%]_[#]ok','%#ok','ignorecase')
function nchanged=fileregexp(files,opt,varargin)
if ischar(files)
    files={files};
end

if ~isfield(opt,'subdir')
    opt.subdir=0;
    %search subdirs
end

if ~isfield(opt,'action')
    opt.action='query';
    %action = show replace query
end

if ~isfield(opt,'noquote')
    opt.noquote=false;
    %skip if in quotes or in comment
end

if ~isfield(opt,'addcounter')
    opt.addcounter=false;
    %add counter to replacement text
end

nchanged=0;
for i=1:length(files)
    file=files{i};
    if strcontains(file,'*')
        %wildcard in filename
        thepath=fileparts(file);
        if isempty(thepath)
            thepath=cd;
        end

        oldpath=cd;
        cd(thepath);
        d=dir(file);
        ndx=~[d(:).isdir];
        d=d(ndx);
        for j=1:length(ndx)
            [opt,changed]=runfile(d(j).name,opt,varargin{:});
            if changed
                nchanged=nchanged+1;
            end

        end

        if opt.subdir
            subdirs=dir('*.');
            ndx=[subdirs(:).isdir]&~strcmp({subdirs(:).name},'.')&~strcmp({subdirs(:).name},'..');
            subdirs=subdirs(ndx);
            for j=1:length(subdirs)
                cd(subdirs(j).name);
                nchanged=nchanged+fileregexp(file,opt,varargin{:});
                cd(oldpath)
            end

        end

        cd(oldpath);
    else
        %no wildcard
        [opt,changed]=runfile(file,opt,varargin{:});
        if changed
            nchanged=nchanged+1;
        end

    end

end


function [opt,changed]=runfile(file,opt,varargin)
fid=fopen(file,'r');
try
    lines=textscan(fid,'%s', 'delimiter','\n', 'whitespace','');
catch err
    if strcontains(err.identifier,'BufferOverflow')
        lines=textscan(fid,'%s', 'delimiter','\n', 'whitespace','', 'bufsize',1000000); %#ok<BUFSIZE>
    else
        rethrow(err);
    end

end
fclose(fid);
changed=false;
lines=lines{1};
if length(varargin)>2&&strcmp(varargin{3},'ignorecase')
    [sstart,send]=regexp(lines,varargin{1},'start','end','ignorecase');
else
    [sstart,send]=regexp(lines,varargin{1},'start','end');
end
ndx=find(~cellfun('isempty',sstart));
if opt.noquote
    %remove if in quotes or in a comment %
    [qstart,qend]=regexp(lines(ndx),'[%].*|[''][^'']*['']','start','end');
    for i=1:length(ndx)
        sstart1=sstart{ndx(i)};
        send1=send{ndx(i)};
        ndx1=false(size(sstart1));
        if ~isempty(qstart{i})
            for n=1:length(qstart{i})
                ndx1=ndx1|sstart1>=qstart{i}(n)&send1<=qend{i}(n);
            end

            sstart{ndx(i)}=sstart1(~ndx1);
            send{ndx(i)}=send1(~ndx1);
        end

    end

    ndx=find(~cellfun('isempty',sstart));
end

if ~isempty(ndx)
    for i=1:length(ndx)
        j=ndx(i);
        if opt.addcounter
            newline=lines{j};
            sstart1=sstart{ndx(i)};
            send1=send{ndx(i)};
            for k=length(sstart1):-1:1
                newline=sprintf('%s%s%d%s',newline(1:sstart1(k)-1),varargin{2},k,newline(send1(k)+1:end));
            end

        else
            newline=regexprep(lines{j},varargin{:});
        end

        
        switch opt.action
            case 'show'
                disp(fullfile(cd,file));
                disp(lines{j});
                disp('replace with')
                disp(newline);
                disp('')
            case 'replace'
                lines{j}=newline;
                changed=true;
            case 'query'
                disp(fullfile(cd,file));
                disp(lines{j});
                disp('replace with')
                disp(newline);
                disp('')
                s=' ';
                while ~any(s=='RSAX')
                    s=upper(input('Enter choice [R]replace [S]skip [A]replace all[X]show all','s'));
                    if length(s)~=1
                        s=' ';
                    end

                end

                if s=='R'
                    lines{j}=newline;
                    changed=true;
                elseif s=='A'
                    lines{j}=newline;
                    changed=true;
                    opt.action='replace';
                elseif s=='X'
                    opt.action='show';
                end

            otherwise
                error('grind:fileregexp','Unknown action')
        end

    end

    if changed
        if strcmp(opt.action,'query')
            s=upper(input('OK to replace file?(Y/N)','s'));
            if isempty(s)||s=='N'
                changed=false;
            end

        end

        if changed
            fid = fopen(file, 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
        end

    end

end
