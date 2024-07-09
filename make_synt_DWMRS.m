%% Make synthetic DW-MR spectra  
% This code was used to generate synthetic DW-MR spectra for
% Genovese and Marjanska MRM (2024) doi:xxx
%
%
% by Guglielmo Genovese, PhD
% readpted for github on July 9, 2024
%%---------------------------------------------------------------

%% general parameters

load('steam_basis.mat')  
n_subj = 90; %number of datasets per diffusion schemes
metabolite = {'Cr','Glu','GPC','Ins','sNAA'};
n_av = 64; %number of averages per shell


%% parameter for mono-exponential signal decay 
% S/S0 = exp(-bADC)

mADC_metabolite = [1.2, 1.4, 1.1, 1.2, 1.2]*1e-4;
mSD_ADC_metabolite = [1.0, 2.0, 1.0, 1.0, 1.0]*1e-5; 

for i = 1:length(metabolite)
    ADC_GT_tmp = mADC_metabolite(i)+mSD_ADC_metabolite(i)*rand(n_subj,1);
    eval(['mADC_GT.' metabolite{i} '= ADC_GT_tmp;']);
end
ADC_GT_tmp = [];

%% parameter for bi-exponential signal decay 
% S/S0 = A*exp(-bADC)+(1-A)*exp(-bB)

bADC_metabolite = [1.8, 2.4, 1.3, 2.5, 1.5]*1e-4;
bSD_ADC_metabolite = [1.0, 3.0, 1.0, 2.0, 1.0]*1e-5; 
A_metabolite = [0.5, 0.4, 0.4, 0.5, 0.4];
B_metabolite = [4.0, 4.0, 3.0, 3.0, 1.0]*1e-5; 

for i = 1:length(metabolite)
    ADC_GT_tmp = bADC_metabolite(i)+bSD_ADC_metabolite(i)*rand(n_subj,1);
    eval(['bADC_GT.' metabolite{i} '= ADC_GT_tmp;']);
end
ADC_GT_tmp = [];

%% make ground-truth fids for expSTEAM (7 b-values, 1 direction)

b_value = [0, 529, 1058, 1586, 2115, 2643, 3172];
for i = 1:length(metabolite)
    eval(['ADC_GT_tmp = mADC_GT.' metabolite{i} ';']); 
    for k = 1:n_subj
        decay_tmp(k,:) = exp(-b_value'*ADC_GT_tmp(k));
    end
    eval(['expSTEAM_decay.' metabolite{i} '=decay_tmp;']);
    decay_tmp = [];
end

expSTEAM_GTfid = simGT(BasisSet,expSTEAM_decay,scale,metabolite,b_value,n_subj); 
expSTEAM_RAWfid = simRAW(expSTEAM_GTfid,b_value,n_subj,n_av); 

%% make ground-truth fids for dirSTEAM (2 b-values, 6 direction)

b_value = [0, 3172*ones(1,6)];
for i = 1:length(metabolite)
    eval(['ADC_GT_tmp = mADC_GT.' metabolite{i} ';']); 
    decay_tmp = decay_dirSTEAM(ADC_GT_tmp ,n_subj,b_value(2));
    eval(['dirSTEAM_decay.' metabolite{i} '=decay_tmp;']);
    decay_tmp = [];
end

dirSTEAM_GTfid = simGT(BasisSet,dirSTEAM_decay,scale,metabolite,b_value,n_subj); 
dirSTEAM_RAWfid = simRAW(dirSTEAM_GTfid,b_value,n_subj,n_av); 

%% make ground-truth fids for biexpSTEAM (11 b-values, 1 direction)

b_value = [0 750 1500 2250 3000 4403 6463 9487 13925 20439 30000];
for i = 1:length(metabolite)
    eval(['ADC_GT_tmp = bADC_GT.' metabolite{i} ';']); 
    for k = 1:n_subj
        decay_tmp(k,:) = A_metabolite(i)*exp(-b_value'*ADC_GT_tmp(k))+(1-A_metabolite(i))*exp(-b_value'*B_metabolite(i));
    end
    eval(['biexpSTEAM_decay.' metabolite{i} '=decay_tmp;']);
    decay_tmp = [];
end

biexpSTEAM_GTfid = simGT(BasisSet,biexpSTEAM_decay,scale,metabolite,b_value,n_subj); 
biexpSTEAM_RAWfid = simRAW(biexpSTEAM_GTfid,b_value,n_subj,n_av);    




%% functions 

function GTfid = simGT(basis_set,S,scalefactor,metab,b_value,n_subj) 

    
    for kk = 1:n_subj

        for b = 1:length(b_value)

            ftF = 0;

            for i=1:length(metab)
                eval(['ftF = ftF + basis_set.' metab{i} '.*scalefactor.' metab{i} '.*S.' metab{i} '(kk,b);']);
            end

            fid = ifft(fftshift(ftF));

            eval(['GTfid.b' int2str(b-1) '(:,kk)=fid(1:2048);']);
        end

    end

end


function RAWfid = simRAW(GTfid,b_value,n_subj,n_av) 

    kn = 3; %noise level
    nc_points = size(GTfid.b0,1);


    for j = 1:n_subj

        tmp_fid = [];
        tmp_fidT = [];

        for b=1:length(b_value)

            eval(['tmp_fid  = repmat(GTfid.b' int2str(b-1) '(:,j),1,n_av) + kn*randn(nc_points,n_av) + 1i*kn*randn(nc_points,n_av);']);
            tmp_fidT = cat(2,tmp_fidT,tmp_fid);

        end

        if j < 10
            eval(['RAWfid.subj0' int2str(j) '=tmp_fidT;']);
        else
            eval(['RAWfid.subj' int2str(j) '=tmp_fidT;']);
        end

    end

end
     

function S = decay_dirSTEAM(mADC_GT,n_subj,b_value) 

    ADCc = 0;
    epsilon = 1e-9;   %error level: 1e-9 um2/ms
    spr = 0.2; %spread level: max 20% among directions  

    for i = 1:n_subj
        target = mADC_GT(i);
        while abs(target-ADCc) > epsilon
            setV = target*((1-spr/2)+spr*rand(1,6));
            ADCc = 1/3*(sqrt(setV(1)*setV(2))+sqrt(setV(3)*setV(4))+sqrt(setV(5)*setV(6)));
        end
        S(i,:) = [1 exp(-b_value*setV)];
    end

end




