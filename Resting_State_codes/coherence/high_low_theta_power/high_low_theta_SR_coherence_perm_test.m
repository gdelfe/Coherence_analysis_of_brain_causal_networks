
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code compute the theta coherence permutation test for the SR
% coherence for high and low modulator power trials
%
%    @ Gino Del Ferraro, July 2022, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data';

name_struct_input = '/session_data_lfp.mat'; % -- name file to load
filename = '.mat'; % -- filename for sess_data_info.mat

freq_band = 'theta_band';
monkey = 'Archie';

dir_high_low_theta = strcat(dir_main,sprintf('/%s/Resting_State/high_low_theta',monkey));
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

% ---- parameters for the coherence-gram
fs = 1000;
fk = 200;
pad = 2;
N = 1;
W = 5;

diff = [];
nperm = 1000; % number of permutation for each session, for each modulator


for s = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    clear send_rec sess_data_lfp
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    load(strcat(dir_Modulators,name_struct_input)); % load sess_data_lfp, structure with session modulator info
    
    if sess_data_lfp.mod_idx ~= sess_data_lfp.receiver_idx % if modulator is not receiver
        
        
        % Load sender-receiver structure
        dir_Sess_send_rec_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/send_rec/Data',Sess));
        load(strcat(dir_Sess_send_rec_data,'/send_rec_coh_for_session.mat'));
        
        % number of modulators in that session
        n_mod = size(send_rec.mod_idx,2);
        cnt_m = 1; % counter for numb of modulators checked for the outlier trials remover
        
        
        
        for m = send_rec.mod_idx
            
            display(['-- Modulator ',num2str(m),'  --- ',num2str(cnt_m),' out of ', num2str(n_mod)]);
            
            % load Sender and Receiver LFP without outliers
            lfp_S = send_rec.mod(cnt_m).lfp_S_clean;
            lfp_R = send_rec.mod(cnt_m).lfp_R_clean;
            
            
            % load low and high theta power indexes
            low_idx = send_rec.mod(cnt_m).low_pow_idx;
            high_idx = send_rec.mod(cnt_m).high_pow_idx;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Coherency Sender-Receiver (Permutation test)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            idx = [high_idx; low_idx];
            L = length(idx);
            
            for p = 1:nperm
                
                perm = idx(randperm(length(idx)));
                high_perm = perm(1:L/2);
                low_perm = perm(L/2+1:end);
                
                % -- coherence for null distribution
                [c_sr_high,f] = coherency(lfp_S(high_perm,:),lfp_R(high_perm,:),[N W],fs,fk,pad,0.05,1,1);
                [c_sr_low,f] = coherency(lfp_S(low_perm,:),lfp_R(low_perm,:),[N W],fs,fk,pad,0.05,1,1);
                
                
                d = abs(c_sr_high) - abs(c_sr_low);
                diff = [diff; d];
                
            end
            cnt_m = cnt_m + 1;        
            
        end
        
    end
    
    
end



save(strcat(dir_high_low_theta,'permutation_test/coh_diff_permutation_mav.mat'),'diff','-v7.3')




