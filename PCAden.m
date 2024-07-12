%% PCA denoising at fixed rank 
% This code was used to perform denoising at fixed rank 
% using the slinding window 
% on synthetic DW-MR spectra for
% Genovese and Marjanska MRM (2024) doi:xxx
%
%
% by Guglielmo Genovese, PhD
% readpted for github on July 12, 2024
%% ---------------------------------------------------------------

function dataset_den = PCAden(dataset,n_av,n_cond,P)
%
% Input: 
% - dataset --> RAW DW-MR spectra generate with 'make_synt_DWMRS.m' 
%               (for any diffusion scheme, i.e., expSTEAM_RAWfid
%               dirSTEAM_RAWfid, and biexpSTEAM_RAWfid)
%
% - n_av --> Number of averages per shell (i.e., n_av = 64)
%
% - n_cond --> Number of shells
%
% - P --> rank at which performs denoising
%
% Output:
% - dataset_den --> Denoised DW-MR spectra
%
%% ---------------------------------------------------------------
%% set parameters
    
    a = fieldnames(dataset);
    n_subj = length(a);

    nex = 0;
    for i = 1:n_cond     
        nex(i+1) = i*n_av;         %% counter for moving along the shells
    end

%%  performing denoising on all the dataset

    for k=1:n_subj

        eval(['fids = dataset.' a{k} ';']);

%% denoising on single experiment

        Xre = real(fids);
        Xim = imag(fids);

        for i = 2:length(nex)-2      % middle shells
            Xtmp = sub_block(Xre,Xim,nex(i-1)+1,nex(i+2));      
            Xtmp_d = lam_pca(Xtmp,P);
            fids_den(:,nex(i)+1:nex(i+1)) = Xtmp_d(:,nex(2)+1:nex(3)) + 1i*Xtmp_d(:,nex(5)+1:nex(6));
        end

        Xtmp = sub_block_bord(Xre,Xim,nex,'left');  %first shell
        Xtmp_d = lam_pca(Xtmp,P);
        fids_den(:,nex(1)+1:nex(2)) = Xtmp_d(:,nex(2)+1:nex(3)) + 1i*Xtmp_d(:,nex(5)+1:nex(6));

        Xtmp = sub_block_bord(Xre,Xim,nex,'right');    %last shell
        Xtmp_d = lam_pca(Xtmp,P);
        fids_den(:,nex(end-1)+1:nex(end)) = Xtmp_d(:,nex(2)+1:nex(3)) + 1i*Xtmp_d(:,nex(5)+1:nex(6));

%%        
        eval(['dataset_den.' a{k} '= fids_den;']);

    end

end

%% functions

function Xsub_block = sub_block(Xre,Xim,st,ed)
    Xsub_block = cat(2,Xre(:,st:ed),Xim(:,st:ed));
end

function Xsub_block = sub_block_bord(Xre,Xim,nex,bord)

    if contains(bord,'left')
        Xsub_tmp_re = cat(2,Xre(:,nex(end-1)+1:nex(end)),Xre(:,nex(1)+1:nex(3)));
        Xsub_tmp_im = cat(2,Xim(:,nex(end-1)+1:nex(end)),Xim(:,nex(1)+1:nex(3)));
    elseif contains(bord,'right') 
        Xsub_tmp_re = cat(2,Xre(:,nex(end-2)+1:nex(end)),Xre(:,nex(1)+1:nex(2)));
        Xsub_tmp_im = cat(2,Xim(:,nex(end-2)+1:nex(end)),Xim(:,nex(1)+1:nex(2)));
    end

    Xsub_block = cat(2,Xsub_tmp_re,Xsub_tmp_im);

end


function X_d = lam_pca(X,P)
    [U,S_Y,V] = svd(X);
    S_Y_tmp = zeros(size(S_Y,1),size(S_Y,2));
    S_Y_tmp(1:P,1:P) = S_Y(1:P,1:P);
    X_d = U*S_Y_tmp*V';
end
    




 



