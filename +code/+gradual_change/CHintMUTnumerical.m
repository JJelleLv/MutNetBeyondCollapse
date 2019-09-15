function CHintMUTnumericalFile = CHintMUTnumerical(NETnr,SETnr,replaceCHintMUTnumericalFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gradually change interaction strengts %%%%
%%%% from initial network to final network %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load initial network
initialNetFolder = sprintf('data%sINITIAL_NETWORKS',filesep);
initialNetFile = sprintf('%s%s%d_INIT_NET',initialNetFolder,filesep,NETnr);

initialNetDATA=load(initialNetFile);

%% load final network
finalSetsFolder = sprintf('data%sFINAL_NETWORKS%s%d_INIT_NET',filesep,filesep,NETnr);
finalSetFile = sprintf('%s%s%d_FINAL_SET',finalSetsFolder,filesep,SETnr);

finalSetDATA=load(finalSetFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Output Folder and File %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% determine file and folder to store in
CHintMUTnumFolder = sprintf('data%sCHintMUTnumerical%s%d_INIT_NET',filesep,filesep,NETnr);
CHintMUTnumericalFile = sprintf('%s%s%d_CHintMUTnumerical',CHintMUTnumFolder,filesep,SETnr);

%% skip if this set already exists
if (~replaceCHintMUTnumericalFile)
    if exist(sprintf('%s.mat',CHintMUTnumericalFile),'file')==2
        fprintf('\nTIME SERIES LOADED FROM FILE:\n')
        fprintf('%s.mat\n',CHintMUTnumericalFile)
        return
    end
end

%% mk folder if doesnt exist
if (~exist(CHintMUTnumFolder, 'dir'))
    mkdir(CHintMUTnumFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% load GRIND (numerical) model %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load model
eval(sprintf('cd +code%s+gradual_change',filesep));
use mut_typeII_20D_10PL_10AN_M_dw
global P A rp ra dp da cp ca hp ha lap lpa dwp dwa ...
    M SIMTIME Peq Peq_null Peq_change_M Aeq Aeq_null Aeq_change_M ...
    RELMUTNETLPA RELMUTNETLPA_null RELMUTNETLPA_change_M RELMUTNETLAP RELMUTNETLAP_null RELMUTNETLAP_change_M ...
    cPeq LFULLPeq LPeq cAeq LFULLAeq LAeq

%% noise level
dw=0.1;

%% time series length
simtime 0 19999 19999
SIMTIME=2e4;

%% solver (Euler step length
solver euler 0.01

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% assign parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initial network
NRPl=initialNetDATA.NRPl;
NRAn=initialNetDATA.NRAn;
rp=initialNetDATA.rp;
ra=initialNetDATA.ra;
dp=initialNetDATA.dp;
da=initialNetDATA.da;
cp=initialNetDATA.cp;
ca=initialNetDATA.ca;
hp=initialNetDATA.hp;
ha=initialNetDATA.ha;

Peq_null=initialNetDATA.Peq_all;
Aeq_null=initialNetDATA.Aeq_all;

RELMUTNETLPA_null=initialNetDATA.RELMUTNETLPA;
RELMUTNETLAP_null=initialNetDATA.RELMUTNETLAP;

%% final network
Peq_end=initialNetDATA.Peq_all;
Aeq_end=initialNetDATA.Aeq_all;

RELMUTNETLPA_end=finalSetDATA.RELMUTNETLPA;
RELMUTNETLAP_end=finalSetDATA.RELMUTNETLAP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAKE CHANGE MAT/VECT %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% change in abundance of NTReq depending on M
Peq_change_M=(Peq_end-Peq_null);
Aeq_change_M=(Aeq_end-Aeq_null);

%% change in relative interaction strengths
RELMUTNETLPA_change_M=(RELMUTNETLPA_end-RELMUTNETLPA_null);
RELMUTNETLAP_change_M=(RELMUTNETLAP_end-RELMUTNETLAP_null);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAKE ALL GLOBAL IN GRIND %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=0;
Peq=Peq_null+Peq_change_M.*M;
RELMUTNETLPA=RELMUTNETLPA_null+RELMUTNETLPA_change_M.*M;
Aeq=Aeq_null+Aeq_change_M.*M;
RELMUTNETLAP=RELMUTNETLAP_null+RELMUTNETLAP_change_M.*M;
cPeq=cp*Peq;
LFULLPeq=-(-rp+dp+cPeq)./(-1-rp.*hp+dp.*hp+cPeq.*hp);
LPeq=RELMUTNETLPA.*(LFULLPeq*ones(1,10));
lpa=LPeq.*(ones(10,1)*(ones(1,10)./Aeq'));
cAeq=ca*Aeq;
LFULLAeq=-(-ra+da+cAeq)./(-1-ra.*ha+da.*ha+cAeq.*ha);
LAeq=RELMUTNETLAP.*(LFULLAeq*ones(1,10));
lap=LAeq.*(ones(10,1)*(ones(1,10)./Peq'));

%%%%%%%%%%%%%%%%%%%%%%
%%%% APPLY CHANGE %%%%
%%%%%%%%%%%%%%%%%%%%%%

%% initial abundances
P=Peq_null;
A=Aeq_null;

%% noise
dwp=ones(NRPl,1).*dw;
dwa=ones(NRAn,1).*dw;

%% run grind
tic
time -s
ke
toc

%% timeseries
Mtime=outfun('M');
Ptime=outfun('P');
Atime=outfun('A');

%% check M=1
if M>1.01 || M<0.99
    M
    error('M is not 1')
end

%%%%%%%%%%%%%%%%%%%%%%
%%%% Store output %%%%
%%%%%%%%%%%%%%%%%%%%%%

%% return to main folder
eval(sprintf('cd ..%s..',filesep));

save(CHintMUTnumericalFile,'Mtime','Ptime','Atime','dwp','dwa','dw', ...
    'RELMUTNETLPA_null','RELMUTNETLPA_end','RELMUTNETLPA_change_M','RELMUTNETLAP_null','RELMUTNETLAP_end','RELMUTNETLAP_change_M', ...
    'Peq_null', 'Peq_end','Peq_change_M','Aeq_null', 'Aeq_end','Aeq_change_M','SIMTIME')

fprintf('CHintMUTnumerical (dw=%.2f):, %d - %d \n',dw,NETnr,SETnr);

%%%%%%%%%%%%%%%%%%%%%%
%%%% CLEAR GLOBAL %%%%
%%%%%%%%%%%%%%%%%%%%%%

clearvars -global