
function orthogonal_proj_combo(data,hitIndx,missIndx,tapers,fs,dn,fk,pad,iSubject,iSess,iCh)

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
%         
% OUTPUT: spectrogram data and figures saved on files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = tapers(1);
w = tapers(2);
k = floor(2*n*w-1);

display(['The number of tapers is k = ',k]) 

SpecK_hit = []; SpecK_miss = []; 

% comput the projection of Xk on the K tapers
[SpecK_hit, fhit] = tfsp_proj_GINO(data(hitIndx,:), tapers,fs, dn, fk, pad); % Trial, Time, K (taper), frequency - complex number
[SpecK_miss, fmiss] = tfsp_proj_GINO(data(missIndx,:), tapers,fs, dn, fk, pad); % Trial, Time, K (taper), frequency - complex number


% Compute the projection of the spectrogram on the K tapers, i.e. |X_k|^2,
% trail, time, K, frequency
Sx_Length_hit = SpecK_hit.*conj(SpecK_hit); % compute the length of the complex number
Sx_Length_miss = SpecK_miss.*conj(SpecK_miss); % compute the length of the complex number


% figure(3); tvimage(log(Sx_Length(1,:,2,:))); colorbar; % for trial 1, taper 2
% Sx_Length_mean = sq(mean(Sx_Length,3)); % average across tapers
% figure(3); tvimage(log(Sx_Length_mean(1,:,:))); colorbar; % spectrogram for trail 1

keyboard 

for i = 1:(k-1) % for all the couple of tapers, i.e. 1,2 and 1,3 and 1,4 and 2,3 and so on
    for j = i+1:k
        
        % %% make directories
        
        display('Creating directories...')
        % Data directories %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Dir_dataHits = sprintf('Data/Hits/%d_Subject/%d_Sess/%d_Ch/orthogonal_%d%d',iSubject,iSess,iCh,i,j)
        if ~exist(Dir_dataHits, 'dir')
            mkdir(Dir_dataHits)
        end
        
        Dir_dataMisses = sprintf('Data/Misses/%d_Subject/%d_Sess/%d_Ch/orthogonal_%d%d',iSubject,iSess,iCh,i,j)
        if ~exist(Dir_dataMisses, 'dir')
            mkdir(Dir_dataMisses)
        end
       
        % generate tensor to store the sum of the orthogonal projections 
        spec_hit = zeros(size(Sx_Length_hit,1),size(Sx_Length_hit,2),1,size(Sx_Length_hit,4)); % tensor for the spectrogram: trial,time,1,frequency
        spec_miss = zeros(size(Sx_Length_miss,1),size(Sx_Length_miss,2),1,size(Sx_Length_miss,4)); % tensor for the spectrogram: trial,time,1,frequency
        
%         i = 1;
%         j = 2;
%         
        % compute the sum of two orthogonal projections 
        spec_hit(:,:,1,:) = 0.5*(Sx_Length_hit(:,:,i,:) + Sx_Length_hit(:,:,j,:));
        spec_miss(:,:,1,:) = 0.5*(Sx_Length_miss(:,:,i,:) + Sx_Length_miss(:,:,j,:));
        
        % squeeze the tensors
        spec_hit = sq(spec_hit); % trial,time,frequency
        spec_miss = sq(spec_miss); % trial,time,frequency
        
        
%         figure(1); tvimage(sq(log(spec_hit(1,:,:)))); title(sprintf('Hit - ')); colorbar;
%         figure(2); tvimage(sq(log(spec_miss(1,:,:)))); title(sprintf('Miss - ')); colorbar;
%         
%         
         
         
         % %%%%% Save data files
         
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
             dlmwrite(fileName_Hit,sq(log(Sx_Length_hit(indx,:,k,:))),'delimiter','\t')
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
             dlmwrite(fileName_Miss,sq(log(Sx_Length_miss(indx,:,k,:))),'delimiter','\t')
         end
         
         % save file with labels of the miss trails
         Name_MissIndx = sprintf('%s/%d_Subject_%d_Sess_%d_Ch_Misses_index.txt',Dir_dataMisses,iSubject,iSess,iCh)
         dlmwrite(Name_MissIndx,missIndx')
             
    end
end
    
  

end