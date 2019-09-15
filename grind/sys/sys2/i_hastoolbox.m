function yes=i_hastoolbox(toolname)
%can simulate here that toolbox is not available
yes=exist(fullfile(matlabroot,'toolbox', toolname),'dir')==7;
%yes=false;

