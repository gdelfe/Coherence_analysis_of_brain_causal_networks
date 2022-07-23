
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
monkey = 'Maverick';

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

diff_mod = [];
diff_ctrl = [];
diff_mod_ctrl = [];

diff_mod_OA = [];
diff_ctrl_OA = [];
diff_mod_ctrl_OA = [];


nperm = 1000; % number of permutation for each session, for each modulator


for s = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    clear send_rec sess_data_lfp
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    load(strcat(dir_Modulators,name_struct_input)); % load sess_data_lfp, structure with session modulator info
    dir_Ctrl_SA =  strcat(dir_RS_Theta,sprintf('/Sess_%d/Controls_same_area',Sess));
    load(strcat(dir_Ctrl_SA,'/sess_controls_same_area_lfp_movie.mat')); % load controls lfp and data
    
    if sess_data_lfp.mod_idx ~= sess_data_lfp.receiver_idx % if modulator is not receiver
        
        
        % Load sender-receiver structure
        dir_Sess_send_rec_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/send_rec/Data',Sess));
        load(strcat(dir_Sess_send_rec_data,'/send_rec_coh_for_session.mat'));
        mod = send_rec;
        load(strcat(dir_Sess_send_rec_data,'/send_rec_coh_for_session_controls_SA.mat'));
        ctrl_SA = send_rec;
        load(strcat(dir_Sess_send_rec_data,'/send_rec_coh_for_session_controls_OA.mat'));
        ctrl_OA = send_rec;
        
        % number of modulators in that session
        n_mod = size(mod.mod_idx,2);
        cnt_m = 1; % counter for numb of modulators checked for the outlier trials remover
        
        
        
        for m = mod.mod_idx
            
            display(['-- Modulator ',num2str(m),'  --- ',num2str(cnt_m),' out of ', num2str(n_mod)]);
            
            % load Sender and Receiver LFP without outliers
            lfp_S = mod.mod(cnt_m).lfp_S_clean;
            lfp_R = mod.mod(cnt_m).lfp_R_clean;
            
            
            % load low and high theta power indexes for the modulator
            low_idx = mod.mod(cnt_m).low_pow_idx;
            high_idx = mod.mod(cnt_m).high_pow_idx;
            
            % load low and high theta power indexes for the control SA
            r = randperm(size(ctrl_SA.mod,2)); % select a random control to be used in comparison to the modulator 
            r = r(1);
            low_idx_SA = ctrl_SA.mod(r).low_pow_idx;
            high_idx_SA = ctrl_SA.mod(r).high_pow_idx;
            
            
            % load low and high theta power indexes for the control OA
            r = randperm(size(ctrl_OA.mod,2)); % select a random control to be used in comparison to the modulator 
            r = r(1);
            low_idx_OA = ctrl_OA.mod(r).low_pow_idx;
            high_idx_OA = ctrl_OA.mod(r).high_pow_idx;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Coherency Sender-Receiver (Permutation test)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            idx_high = [high_idx; high_idx_SA]; % concatenate high pow trial for modulator and controls SA
            idx_low = [low_idx; low_idx_SA];    % concatenate low pow trial for modulator and controls SA
            
            idx_high_mod_OA = [high_idx; high_idx_OA]; % concatenate high pow trial for modulator and controls OA
            idx_low_mod_OA = [low_idx; low_idx_OA];    % concatenate low pow trial for modulator and controls OA
            
            L_high = length(idx_high);
            L_low = length(idx_low);
            
            for p = 1:nperm
                
                perm_high = idx_high(randperm(L_high)); % permute high idx trial for control SA and modulator
                perm_low = idx_low(randperm(L_low));    % permute low idx trial for control SA and modulator
                
                perm_high_mod_OA = idx_high_mod_OA(randperm(L_high)); % permute high idx trial for control OA and modulator
                perm_low_mod_OA = idx_low_mod_OA(randperm(L_low));    % permute low idx trial for control OA and modulator
                
                % pseudo trials high and low theta, for controls SA and modulators 
                pseudo_H_mod = perm_high(1:floor(L_high/2));
                pseudo_H_ctrl = perm_high(ceil(L_high/2)+1:end);
                
                pseudo_L_mod = perm_low(1:floor(L_low/2));
                pseudo_L_ctrl = perm_low(ceil(L_low/2)+1:end);
                
                
                % pseudo trials high and low theta, for controls OA and modulators 
                pseudo_H_mod_OA = perm_high_mod_OA(1:floor(L_high/2));
                pseudo_H_ctrl_OA = perm_high_mod_OA(ceil(L_high/2)+1:end);
                
                pseudo_L_mod_OA = perm_low_mod_OA(1:floor(L_low/2));
                pseudo_L_ctrl_OA = perm_low_mod_OA(ceil(L_low/2)+1:end);
                
                
                
                % -- coherence for null distribution
                % %% modulator and control SA
                % modulator
                [c_sr_high_mod,f] = coherency(lfp_S(pseudo_H_mod,:),lfp_R(pseudo_H_mod,:),[N W],fs,fk,pad,0.05,1,1);
                [c_sr_low_mod,f] = coherency(lfp_S(pseudo_L_mod,:),lfp_R(pseudo_L_mod,:),[N W],fs,fk,pad,0.05,1,1);
                % control SA
                [c_sr_high_ctrl,f] = coherency(lfp_S(pseudo_H_ctrl,:),lfp_R(pseudo_H_ctrl,:),[N W],fs,fk,pad,0.05,1,1);
                [c_sr_low_ctrl,f] = coherency(lfp_S(pseudo_L_ctrl,:),lfp_R(pseudo_L_ctrl,:),[N W],fs,fk,pad,0.05,1,1);
                
                
                % %% modulator and control OA
                % modulator
                [c_sr_high_mod_OA,f] = coherency(lfp_S(pseudo_H_mod_OA,:),lfp_R(pseudo_H_mod_OA,:),[N W],fs,fk,pad,0.05,1,1);
                [c_sr_low_mod_OA,f] = coherency(lfp_S(pseudo_L_mod_OA,:),lfp_R(pseudo_L_mod_OA,:),[N W],fs,fk,pad,0.05,1,1);
                % control SA
                [c_sr_high_ctrl_OA,f] = coherency(lfp_S(pseudo_H_ctrl_OA,:),lfp_R(pseudo_H_ctrl_OA,:),[N W],fs,fk,pad,0.05,1,1);
                [c_sr_low_ctrl_OA,f] = coherency(lfp_S(pseudo_L_ctrl_OA,:),lfp_R(pseudo_L_ctrl_OA,:),[N W],fs,fk,pad,0.05,1,1);
                
                % Modulator control SA
                
                d_mod = abs(c_sr_high_mod) - abs(c_sr_low_mod); % coherence difference modulator
                d_ctrl = abs(c_sr_high_ctrl) - abs(c_sr_low_ctrl); % coherence difference control
                
                % difference of the difference 
                d_mod_ctrl = d_mod - d_ctrl;
                
                % store values 
                diff_mod = [diff_mod; d_mod];
                diff_ctrl = [diff_ctrl; d_ctrl];
                diff_mod_ctrl = [diff_mod_ctrl; d_mod_ctrl];
                
                % % % Modulator control OA

                d_mod_OA = abs(c_sr_high_mod_OA) - abs(c_sr_low_mod_OA); % coherence difference modulator
                d_ctrl_OA = abs(c_sr_high_ctrl_OA) - abs(c_sr_low_ctrl_OA); % coherence difference control
                
                % difference of the difference 
                d_mod_ctrl_OA = d_mod_OA - d_ctrl_OA;
                
                % store values 
                diff_mod_OA = [diff_mod_OA; d_mod_OA];
                diff_ctrl_OA = [diff_ctrl_OA; d_ctrl_OA];
                diff_mod_ctrl_OA = [diff_mod_ctrl_OA; d_mod_ctrl_OA];
                
                
                
                
            end
            cnt_m = cnt_m + 1;        
            
        end
        
    end
    
    
end

save(strcat(dir_high_low_theta,'/coh_diff_modulators_high_low_mod_SA_mav.mat'),'diff_mod','-v7.3')
save(strcat(dir_high_low_theta,'/coh_diff_ctlr_SA_high_low_mod_SA_mav.mat'),'diff_ctlr','-v7.3')
save(strcat(dir_high_low_theta,'/coh_diff_of_diff_mod_ctrl_SA_mav.mat'),'diff_mod_ctrl','-v7.3')


save(strcat(dir_high_low_theta,'/coh_diff_modulators_high_low_mod_OA_mav.mat'),'diff_mod_OA','-v7.3')
save(strcat(dir_high_low_theta,'/coh_diff_ctlr_SA_high_low_mod_OA_mav.mat'),'diff_ctlr_OA','-v7.3')
save(strcat(dir_high_low_theta,'/coh_diff_of_diff_mod_ctrl_OA_mav.mat'),'diff_mod_ctrl_OA','-v7.3')




