function ys=readbin(outfile)
%help function for debugging
if nargin==0
    outfile='out1.tmp';
end
fid=fopen(outfile,'r');
signature=fread(fid,1,'int32');
if signature==1234567890
    siz=fread(fid,2,'int32');
    ys=fread(fid,siz','double');
    ys=ys';
end
fclose(fid);
