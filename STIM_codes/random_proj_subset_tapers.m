
function random_proj_subset_tapers(data,hitIndx,missIndx,tapers,fs,dn,fk,pad,iSubject,iSess,iCh,R)

% HITS_MISSES_SPECTROGRAMS_GINO
%
% Compute spectrograms of hit and miss trail for a given Subject, a given
% Session, and a given Channel, projected on a K (taper) subspace
% Creates directories where to store the data and figures and save each spectrogram
% into data file (time x frequency) and
% as image (jpg)
%
% INPUTS:
%         data: Lfp in trial x time form
%         hitIndx: array that contains the labels of the hit trails
%         missIndx: array that contains the labels of the miss trails
%         tapers: [N,W] form
%         fs:
%         dn:
%         fk:
%         pad:
%         iSubject: animal number
%         iSessio: session number
%         iCh: channel number
%         R: number of total random projections
%
% OUTPUT: spectrogram data and figures saved on files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = tapers(1);
w = tapers(2);
k = floor(2*n*w-1);

SpecK_hit = []; SpecK_miss = [];

% Projection of Xk on the K tapers
[SpecK_hit, fhit] = tfsp_proj_GINO(data(hitIndx,:), tapers,fs, dn, fk, pad); % Trial, Time, K (taper), frequency - complex number
[SpecK_miss, fmiss] = tfsp_proj_GINO(data(missIndx,:), tapers,fs, dn, fk, pad); % Trial, Time, K (taper), frequency - complex number


% Projection of the spectrogram on the K tapers, i.e. |X_k|^2,
% trail, time, K, frequency
Sx_Length_hit = SpecK_hit.*conj(SpecK_hit); % compute the length of the complex number
Sx_Length_miss = SpecK_miss.*conj(SpecK_miss); % compute the length of the complex number


% figure(3); tvimage(log(Sx_Length(1,:,2,:))); colorbar; % for trial 1, taper 2
% Sx_Length_mean = sq(mean(Sx_Length,3)); % average across tapers
% figure(3); tvimage(log(Sx_Length_mean(1,:,:))); colorbar; % spectrogram for trail 1

k = 2 % numb of tapers for each average projection

for rindx = 0: R % for each taper
    
    % generate random vector with k (# tapers) element with fixed sum = 1
    S = randfixedsum(k,1,1,0,1); % row, columns, tot_sum, lower bound for rand, max bound for rand
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate spectrogram with inhomogeneous random weight with fixed sum = 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for ks = 1:2:6      % for each pair of random projection 
     
        spec_hit = zeros(size(Sx_Length_hit,1),size(Sx_Length_hit,2),1,size(Sx_Length_hit,4)); % tensor for the spectrogram: trial,time,1,frequency
        spec_miss = zeros(size(Sx_Length_miss,1),size(Sx_Length_miss,2),1,size(Sx_Length_miss,4)); % tensor for the spectrogram: trial,time,1,frequency
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the sum of the random projection
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cnt = 1;
        for kindx = ks:ks+1
            ks 
            type = sprintf('random_N_015_W_15_dn_018/tapers/tapers_%d_%d/random_%d',ks,ks+1,rindx); % name of the folder
            display(['type: ',num2str(type)])
            
            spec_hit(:,:,1,:) = spec_hit(:,:,1,:) + Sx_Length_hit(:,:,kindx,:)*S(cnt);
            spec_miss(:,:,1,:) = spec_miss(:,:,1,:) + Sx_Length_miss(:,:,kindx,:)*S(cnt);
            cnt = cnt + 1;
        end
        spec_hit = sq(spec_hit); % trial,time,frequency
        spec_miss = sq(spec_miss); % trial,time,frequency
        
%       keyboard
        
        % % cntf = 0
        % % % indx = 1
        %
        % cntf = cntf + 1
        % figure(cntf); tvimage(log(spec_hit(indx,:,:))); colorbar; title('Hit - taper 56 ') % for trial 1, taper 2
        % cntf = cntf + 1
        % figure(cntf); tvimage(log(spec_miss(indx,:,:))); colorbar; title('Miss - taper 56') % for trial 1, taper 2
        %
        %
        %% make directories for data and figures
        
        display('Creating directories...')
        % Data directories %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Dir_dataHits = sprintf('Data/Hits/%d_Subject/%d_Sess/%d_Ch/%s',iSubject,iSess,iCh,type)
        if ~exist(Dir_dataHits, 'dir')
            mkdir(Dir_dataHits)
        end
        
        Dir_dataMisses = sprintf('Data/Misses/%d_Subject/%d_Sess/%d_Ch/%s',iSubject,iSess,iCh,type)
        if ~exist(Dir_dataMisses, 'dir')
            mkdir(Dir_dataMisses)
        end
       
                
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
            dlmwrite(fileName_Hit,sq(log(spec_hit(indx,:,:))),'delimiter','\t')
            
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
            dlmwrite(fileName_Miss,sq(log(spec_miss(indx,:,:))),'delimiter','\t')
            
        end
        
        % save file with labels of the miss trails
        Name_MissIndx = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_Misses_index.txt',Dir_dataMisses,iSubject,iSess,iCh)
        dlmwrite(Name_MissIndx,missIndx')
        
    end
end

end