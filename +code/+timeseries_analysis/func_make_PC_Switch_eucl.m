function [Switch_vect]=func_make_PC_Switch_eucl(DATA_vect)

NR_steps=length(DATA_vect(1,:));

%% determine when to switch PC1 based on Euclidian distance.
Switch_vect=nan(NR_steps,1);
Switch_vect(1,1)=1;
for lstepNR_noise=2:NR_steps
    
    %% Eucledian distance between this and previous step
    EUCLdist_nonreversed=(sum((DATA_vect(:,(lstepNR_noise-1))-DATA_vect(:,(lstepNR_noise))).^2)).^0.5;
    EUCLdist_reversed=(sum((DATA_vect(:,(lstepNR_noise-1))+DATA_vect(:,(lstepNR_noise))).^2)).^0.5;
    
    %% assign direction based on euclidian distance
    Switch_vect(lstepNR_noise,1)=Switch_vect((lstepNR_noise-1),1);
    if EUCLdist_reversed<=EUCLdist_nonreversed
        Switch_vect(lstepNR_noise,1)=-Switch_vect((lstepNR_noise-1),1);
    end
    
end