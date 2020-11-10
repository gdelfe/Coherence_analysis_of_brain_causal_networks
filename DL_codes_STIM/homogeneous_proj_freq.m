
function homogeneous_proj_freq(data,hitIndx,missIndx,bn,tapers,fs,dn,fk,pad,nPerm,alpha,iSubject,iSess,iCh)

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

% Compute the spectrograms 
[HitSpec,fhit] = tfspec_GINO(data(hitIndx,:),tapers,fs,dn,fk,pad,0.05,0,1); % get the time-frequency spectrum of the hit trial. Format trial x time x freq
[MissSpec,fmiss] = tfspec_GINO(data(missIndx,:),tapers,fs,dn,fk,pad,0.05,0,1); % get the time-frequency spectrum of the miss trial. Format trial x time x freq 

HitSpec = HitSpec(:,:,1:61);
MissSpec = MissSpec(:,:,1:61);

% keyboard
% figure(1); tvimage(log(sq(HitSpec(1,:,:)))); colorbar
% figure(1); tvimage(sq(log(specA))); colorbar
% dlmwrite('homogeneous_hit_test.txt',sq(log(HitSpec(1,:,:))),'delimiter','\t')
% figure(2); tvimage(specA); colorbar

%% make directories for data and figures 

display('Creating directories...')
% Data directories %%%%%%%%%%%%%%%%%%%%%%%%%%%
Dir_dataHits = sprintf('Data/Hits/%d_Subject/%d_Sess/%d_Ch/homogeneous_sum_fmin_%d',iSubject,iSess,iCh,fk(1))
if ~exist(Dir_dataHits, 'dir')
    mkdir(Dir_dataHits)
end

Dir_dataMisses = sprintf('Data/Misses/%d_Subject/%d_Sess/%d_Ch/homogeneous_sum_fmin_%d',iSubject,iSess,iCh,fk(1))
if ~exist(Dir_dataMisses, 'dir')
    mkdir(Dir_dataMisses)
end

% Figures directories %%%%%%%%%%%%%%%%%%%%%%%%%% 
% Dir_FiguresHits = sprintf('Figures/Hits/%d_Subject_%d_Sess_%d_Ch',iSubject,iSess,iCh)
% if ~exist(Dir_FiguresHits, 'dir')
%     mkdir(Dir_FiguresHits)
% end
% 
% Dir_FiguresMisses = sprintf('Figures/Misses/%d_Subject_%d_Sess_%d_Ch',iSubject,iSess,iCh)
% if ~exist(Dir_FiguresMisses, 'dir')
%     mkdir(Dir_FiguresMisses)
% end

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
    fileName_Hit = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_%d_hit.txt',Dir_dataHits,iSubject,iSess,iCh,indx);
    dlmwrite(fileName_Hit,sq(log(HitSpec(indx,:,:))),'delimiter','\t')
    % save the figure 
%     figure('visible','off')
%     figHit = tvimage(log(sq(HitSpec(indx,:,:)))); colorbar;
%     figName = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_%d_figHit.jpg',Dir_FiguresHits,iSubject,iSess,iCh,indx);
%     saveas(figHit,figName)
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
    fileName_Miss = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_%d_miss.txt',Dir_dataMisses,iSubject,iSess,iCh,indx);
    dlmwrite(fileName_Miss,sq(log(MissSpec(indx,:,:))),'delimiter','\t')
    
    % save the figure 
%     figure('visible','off')
%     figMiss = tvimage(log(sq(HitSpec(indx,:,:)))); colorbar;
%     figName = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_%d_figMiss.jpg',Dir_FiguresMisses,iSubject,iSess,iCh,indx);
%     saveas(figMiss,figName)
end
% save file with labels of the miss trails
Name_MissIndx = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_Misses_index.txt',Dir_dataMisses,iSubject,iSess,iCh)
dlmwrite(Name_MissIndx,missIndx')

end