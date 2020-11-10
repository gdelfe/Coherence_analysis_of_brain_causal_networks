
function homogeneous_proj_complex(data,hitIndx,missIndx,tapers,fs,dn,fk,pad,iSubject,iSess,iCh)

% HITS_MISSES_SPECTROGRAMS_GINO
%
% Compute spectrograms of hit and miss trail for a given Subject, a given Session, and a given Channel.
% Creates directories where to store the data and figures and save each spectrogram into data file (time x frequency) and
% as image (jpg)
%
% INPUTS:
%         data: Lfp in trial x time form
%         hitIndx: array that contains the labels of the hit trails
%         missIndx: array that contains the labels of the miss trails
%         bn:
%         tapes:
%         fs:
%         dn:
%         fk:
%         pad:
%         nPerm:
%         alpha:
%         iSubject: animal number
%         iSessio: session number
%         iCh: channel number
%
% OUTPUT: spectrogram data and figures saved on files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HitSpec = []; MissSpec = [];

% Compute the complex projections onto the tapers
[HitSpec,fhit] = tfsp_proj_GINO(data(hitIndx,:), tapers,fs, dn, fk, pad);
[MissSpec,fhit] = tfsp_proj_GINO(data(missIndx,:), tapers,fs, dn, fk, pad);

HitSpec_Re = mean(real(HitSpec),3); % average the real part over the tapers
HitSpec_Im = mean(imag(HitSpec),3); % average the real part over the tapers

MissSpec_Re = mean(real(MissSpec),3); % average the imaginary part over the tapers
MissSpec_Im = mean(imag(MissSpec),3); % average the imaginary part over the tapers
% 
% figure(1); tvimage(sq(MissSpec_Re(1,:,1,:))); colorbar
% dlmwrite('homogeneous_hit_test.txt',sq(log(HitSpec(1,:,:))),'delimiter','\t')
% figure(2); tvimage(sq(MissSpec_Im(1,:,1,:))); colorbar


%% make directories for data and figures

display('Creating directories...')
% Data directories %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dir_dataHits = sprintf('Data/Hits/%d_Subject/%d_Sess/%d_Ch/homogeneous_sum',iSubject,iSess,iCh)
Dir_dataHits = sprintf('../Shaoyu_data/Data/Hits/%d_Subject/%d_Sess/%d_Ch/complex/N_015_W_25_dn_018',iSubject,iSess,iCh)

if ~exist(Dir_dataHits, 'dir')
    mkdir(Dir_dataHits)
end

% Dir_dataMisses = sprintf('Data/Misses/%d_Subject/%d_Sess/%d_Ch/homogeneous_sum',iSubject,iSess,iCh)
Dir_dataMisses = sprintf('../Shaoyu_data/Data/Misses/%d_Subject/%d_Sess/%d_Ch/complex/N_015_W_25_dn_018',iSubject,iSess,iCh)
if ~exist(Dir_dataMisses, 'dir')
    mkdir(Dir_dataMisses)
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Save data files
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HITS TRAILS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

display('Saving HIT TRAILS...')
% save spectrogram data and figure for each of the hit trail
for indx = 1:size(hitIndx,2)
    
    if mod(indx,10) == 0
        display(['Saving hit trail # ',num2str(indx)]);
    end
    % save data
    fileName_Hit_Re = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_%d_hit_real.txt',Dir_dataHits,iSubject,iSess,iCh,indx);
    dlmwrite(fileName_Hit_Re,sq(HitSpec_Re(indx,:,1,:)),'delimiter','\t')
    fileName_Hit_Im = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_%d_hit_imaginary.txt',Dir_dataHits,iSubject,iSess,iCh,indx);
    dlmwrite(fileName_Hit_Im,sq(HitSpec_Im(indx,:,1,:)),'delimiter','\t')

end

% save file with labels of the hit trails
Name_HitIndx = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_Hits_index.txt',Dir_dataHits,iSubject,iSess,iCh)
dlmwrite(Name_HitIndx,hitIndx')

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MISS TRAILS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

display('Saving MISS TRAILS...')

% save spectrogram data and figure for each of the hit trail
for indx = 1:size(missIndx,2)
    
    if mod(indx,10) == 0
        display(['Saving miss trail # ',num2str(indx)]);
    end
    % save data
    fileName_Miss_Re = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_%d_miss_real.txt',Dir_dataMisses,iSubject,iSess,iCh,indx);
    dlmwrite(fileName_Miss_Re,sq(MissSpec_Re(indx,:,1,:)),'delimiter','\t')
    fileName_Miss_Im = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_%d_miss_imaginary.txt',Dir_dataMisses,iSubject,iSess,iCh,indx);
    dlmwrite(fileName_Miss_Im,sq(MissSpec_Im(indx,:,1,:)),'delimiter','\t')

end
% save file with labels of the miss trails
Name_MissIndx = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_Misses_index.txt',Dir_dataMisses,iSubject,iSess,iCh)
dlmwrite(Name_MissIndx,missIndx')



end