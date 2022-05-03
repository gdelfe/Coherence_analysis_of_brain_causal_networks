
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the theta-coherence as a function of the
% modulator score, where the modulator score is computed -in this case- as
% the modulator theta-power in the [4,8] Hz range
%
%    @ Gino Del Ferraro, May 2022, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

set(0,'DefaultFigureVisible','off')
% set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat 
recording = 'movie';

freq_band = 'theta_band';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONKEY MAVERICK 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

cnt_e = 1;
score_mat_mav = [];
for s = 1:size(sess_info{1},1)  % For each session with at least one modulator

    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_RS_Theta,sprintf('/Sess_%d',Sess));
    
%     dir_Sess_mod_rec_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_rec/Data',Sess));
    dir_Mod_recording = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/%s',Sess,recording));
    
    % %%%%%%%%% MODULATORS  %%%%%%
    load(strcat(dir_Mod_recording,sprintf('/sess_data_lfp_coherence_fk_200_W_5_%s.mat',recording))); % structure mod
  
        cnt_m = 1;
        for m = sess_data_lfp.mod_idx 
            
            if m ~= sess_data_lfp.receiver_idx % if modulator is different from receiver 
                  
            lfp_E = sq(sess_data_lfp.lfp_E(m,:,:));
            outliers_E = sess_data_lfp.outliers_E(cnt_m).idx;
            lfp_E(outliers_E,:) = [];

            % Compute the spectrum for each trial. Format: iTrial x times
            W = 3;
            [spec, f, err] = dmtspec(lfp_E,[1000/1e3,W],1e3,200,2,0.5,1);
            theta_score = mean(log(spec(9:18)));
            score_mat_mav = [score_mat_mav; double(theta_score), double(m), double(Sess), double(cnt_e)];
            cnt_e = cnt_e + 1;
            
            end 
            
        end    
        
end

score_sort_mav = sort(score_mat_mav,1,'descend')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONKEY ARCHIE  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

monkey = 'Archie';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

cnt_e = 1;
score_mat_arc = [];
for s = 1:size(sess_info{1},1)  % For each session with at least one modulator

    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_RS_Theta,sprintf('/Sess_%d',Sess));
    
%     dir_Sess_mod_rec_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_rec/Data',Sess));
    dir_Mod_recording = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/%s',Sess,recording));
    
    % %%%%%%%%% MODULATORS  %%%%%%
    load(strcat(dir_Mod_recording,sprintf('/sess_data_lfp_coherence_fk_200_W_5_%s.mat',recording))); % structure mod
  
        cnt_m = 1;
        for m = sess_data_lfp.mod_idx
                  
            if m ~= sess_data_lfp.receiver_idx % if modulator is different from receiver

                            
            lfp_E = sq(sess_data_lfp.lfp_E(m,:,:));
            outliers_E = sess_data_lfp.outliers_E(cnt_m).idx;
            lfp_E(outliers_E,:) = [];

            % Compute the spectrum for each trial. Format: iTrial x times
            W = 3;
            [spec, f, err] = dmtspec(lfp_E,[1000/1e3,W],1e3,200,2,0.5,1);
            theta_score = mean(log(spec(9:18)));
            score_mat_arc = [score_mat_arc; double(theta_score), double(m), double(Sess), double(cnt_e)];
            cnt_e = cnt_e + 1;
            
            end 
            
        end    
        
end

score_sort_arc = sort(score_mat_arc,1,'descend')

