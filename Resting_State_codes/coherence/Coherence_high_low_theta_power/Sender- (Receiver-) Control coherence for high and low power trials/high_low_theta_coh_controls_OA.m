
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code compute the theta coherence between the controls (Other Area) and
% sender/receiver for modulators having high/low theta power.
%
%    @ Gino Del Ferraro, March 2022, Pesaran lab, NYU
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
dir_high_low_theta = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Maverick/Resting_State/high_low_theta';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


% ---- parameters for the coherence-gram
nt = 150001;
fs = 1000;
fk = 200;
pad = 2;
N = 1;


coh_cs_high = [];
coh_cs_low = [];
coh_cr_high = [];
coh_cr_low = [];

for s = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    dir_Sess_mod_send_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_send/Data',Sess));
    dir_Ctrl =  strcat(dir_RS_Theta,sprintf('/Sess_%d/Controls_other_areas',Sess));
    
    load(strcat(dir_Modulators,name_struct_input)); % load sess_data_lfp, structure with session modulator info
    load(strcat(dir_Ctrl,'/sess_controls_other_areas_lfp_movie.mat')); % load controls lfp and data
    
    % outliers time series in sender and receiver
    outliers_S = sess_data_lfp.outliers_S;
    outliers_R = sess_data_lfp.outliers_R;
    
    lfp_S = sess_control_lfp.lfp_S;
    lfp_R = sess_control_lfp.lfp_R;
    
    
    mod_Ch = sess_data_lfp.mod_idx;
    RecordPairMRIlabels = sess_data_lfp.RecordPairMRIlabels; % -- MRI labels of the recorder pars
    MRIlabels = sess_data_lfp.MRIlabels; % -- all the available MRI labels
    receiver_idx = sess_data_lfp.receiver_idx; % -- receiver idx
    
    cnt_m = 1;
    for m = sess_data_lfp.mod_idx % for each modulator
        
        if m ~= receiver_idx % if modulator is not receiver
            display(['-- Modulator ',num2str(m),' -- : ',num2str(cnt_m),', out of tot  ',num2str(size(sess_data_lfp.mod_idx,2)),'  '])
            
            lfp_E = sq(sess_data_lfp.lfp_E(m,:,:));
            outliers_E = sess_data_lfp.outliers_E(cnt_m).idx; % trial with artifact for modulator
            
            % Compute the spectrum for each trial. Format: iTrial x times
            W = 3;
            [spec, f, err] = dmtspec(lfp_E,[1000/1e3,W],1e3,200);
            
            % Find low and high theta from the spectrum
            theta_pow = log(mean(spec(:,9:19),2)); % average the spectrum around theta frequencies (9:19) is the idx for theta range
            theta_pow_mean = mean(theta_pow); % get the average theta power
            theta_pow = theta_pow - theta_pow_mean; % rescale the theta power by removing the mean value
            
            [sort_theta, trial_idx] = sort(theta_pow); % sort theta power in ascending order
            cut = fix(length(theta_pow)/3); % find the index of 1/3 of the distribution
            
            % low and high theta power indexes
            low_theta_pow = sort_theta(1:cut);
            low_idx = trial_idx(1:cut);
            
            high_theta = sort_theta(end-cut+1:end);
            high_idx = trial_idx(end-cut+1:end);
            
            % get control other area index: for the case of other areas,
            % the controls are shared across modulators. This is because given a modulator, its region and the other modulators' region
            % are also escluded from the counting of the controls' region 
            ctrl_idx = sess_control_lfp.ctrl_idx;
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%
            % In relation to the sender
            
            for ctrl = ctrl_idx % for all the controls associated with the modulator
                
                cnt_c = find(sess_control_lfp.ctrl_idx == ctrl); % find the position of the control in the list
                lfp_C = sq(sess_control_lfp.lfp_E(ctrl,:,:)); % load control lfp
                
                
                outliers_C = sess_control_lfp.outliers_E(cnt_c).idx; % control outliers/artifact trials
                
                % SENDER - CONTROLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                outliers_CSE = [outliers_S, outliers_C, outliers_E];
                outliers_CSE = unique(outliers_CSE);  % -- remove repeated entries in outliers
                
                
                % COMPUTE HIGH/LOW THETA POWER FOR THE CONTROLS %%%%%%%
                % Compute the spectrum for each trial. Format: iTrial x times
                W = 3;
                [spec_C, f, err] = dmtspec(lfp_C,[1000/1e3,W],1e3,200);
                
                % Find low and high theta from the spectrum
                theta_pow_C = log(mean(spec_C(:,9:19),2)); % average the spectrum around theta frequencies (9:19) is the idx for theta range
                theta_pow_mean_C = mean(theta_pow_C); % get the average theta power
                theta_pow_C = theta_pow_C - theta_pow_mean_C; % rescale the theta power by removing the mean value
                
                [sort_theta_C, trial_idx_C] = sort(theta_pow_C); % sort theta power in ascending order
                cut_C = fix(length(theta_pow_C)/3); % find the index of 1/3 of the distribution
                
                % low and high theta power indexes
                low_theta_pow_C = sort_theta_C(1:cut_C);
                low_idx_C = trial_idx_C(1:cut_C);
                
                high_theta_C = sort_theta_C(end-cut_C+1:end);
                high_idx_C = trial_idx_C(end-cut_C+1:end);
                
                high_C = setdiff(high_idx_C,outliers_CSE); % remove trials with artifacts from high power trials
                low_C = setdiff(low_idx_C,outliers_CSE);    % remove trials with artifacts from low power trials
                
                % compute coherence between sender-control for high/low theta power trials
                W = 5;
                [c_cs_high,f] = coherency(lfp_S(high_C,:),lfp_C(high_C,:),[N W],fs,fk,pad,0.05,1,1);
                [c_cs_low,f] = coherency(lfp_S(low_C,:),lfp_C(low_C,:),[N W],fs,fk,pad,0.05,1,1);
                
                coh_cs_high = [coh_cs_high; c_cs_high];
                coh_cs_low = [coh_cs_low; c_cs_low];
                
                
                
                % RECEIVER - CONTROLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                outliers_CRE = [outliers_R, outliers_C, outliers_E];
                outliers_CRE = unique(outliers_CRE);  % -- remove repeated entries in outliers
                
                high_C = setdiff(high_idx_C,outliers_CRE); % remove trials with artifacts from high power trials
                low_C = setdiff(low_idx_C,outliers_CRE);    % remove trials with artifacts from low power trials
                
                % compute coherence between receiver-control for high/low theta power trials
                [c_cr_high,f] = coherency(lfp_R(high_C,:),lfp_C(high_C,:),[N W],fs,fk,pad,0.05,1,1);
                [c_cr_low,f] = coherency(lfp_R(low_C,:),lfp_C(low_C,:),[N W],fs,fk,pad,0.05,1,1);
                
                coh_cr_high = [coh_cr_high; c_cr_high];
                coh_cr_low = [coh_cr_low; c_cr_low];
                
            end
            
            cnt_m = cnt_m + 1;
        end
    end
end


ctrl_coh.cs_high = coh_cs_high;
ctrl_coh.cs_low = coh_cs_low;
ctrl_coh.cr_high = coh_cr_high;
ctrl_coh.cr_low = coh_cr_low;



keyboard;



save(strcat(dir_high_low_theta,'/coh_all_sess_controls_OA_high-low_pow_controls.mat'),'ctrl_coh')



% load(strcat(dir_high_low_theta,'/coh_all_sess_ms_high.mat'))
% load(strcat(dir_high_low_theta,'/coh_all_sess_ms_low.mat'))
% load(strcat(dir_high_low_theta,'/coh_all_sess_mr_high.mat'))
% load(strcat(dir_high_low_theta,'/coh_all_sess_mr_low.mat'))


mean_all_coh_ms_high = mean(abs(coh_cs_high),1);
mean_all_coh_ms_low = mean(abs(coh_cs_low),1);
mean_all_coh_mr_high = mean(abs(coh_cr_high),1);
mean_all_coh_mr_low = mean(abs(coh_cr_low),1);

err_all_coh_ms_high = std(abs(coh_cs_high),0,1)/sqrt(size(coh_cs_high,1));
err_all_coh_ms_low = std(abs(coh_cs_low),0,1)/sqrt(size(coh_cs_low,1));
err_all_coh_mr_high = std(abs(coh_cr_high),0,1)/sqrt(size(coh_cr_high,1));
err_all_coh_mr_low = std(abs(coh_cr_low),0,1)/sqrt(size(coh_cr_low,1));


% zscore_ms_high = (abs(coh_all_c_ms_high) - mean_all_coh_ms_high)./std(abs(coh_all_c_ms_high),0,1);
% zscore_ms_low = (abs(coh_all_c_ms_low) - mean_all_coh_ms_low)./std(abs(coh_all_c_ms_low),0,1);
% zscore_mr_high = (abs(coh_all_c_mr_high) - mean_all_coh_mr_high)./std(abs(coh_all_c_mr_high),0,1);
% zscore_mr_low = (abs(coh_all_c_mr_low) - mean_all_coh_mr_low)./std(abs(coh_all_c_mr_low),0,1);


% mean_zscore_ms_high = mean(zscore_ms_high,1);
% figure;
% plot(f,mean_all_coh_ms_high)
% hold on
% plot(f,abs(coh_all_c_ms_high))
%

%%%%%%%%%%%%%%%%%%%%%%
% FIGURES %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%


% %%%% RECEIVER %%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureVisible','on')
%     load(strcat(dir_high_low_theta,sprintf('/Sess_%d/coherence_all.mat')));

fig = figure;
shadedErrorBar(f,mean_all_coh_mr_high,err_all_coh_mr_high,'lineprops',{'color',[28 199 139]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_mr_low,err_all_coh_mr_low,'lineprops',{'color',[50 250 93]/255 },'patchSaturation',0.5);

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
title('OA: Coherence MR high vs low theta power trial of the controls','FontSize',12)
legend('high theta pow','low theta pow','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.25])
grid on


fname = strcat(dir_high_low_theta,sprintf('/MR_controls_OA_coherence_mean_high-low_pow_controls.jpg',cnt_m));
saveas(fig,fname);


% %%%% SENDER %%%%%%%%%%%%%%%%%%%%%%


fig = figure;
hold all

shadedErrorBar(f,mean_all_coh_ms_high,err_all_coh_ms_high,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_ms_low,err_all_coh_ms_low,'lineprops',{'color',[255, 51, 153]/255},'patchSaturation',0.4); hold on

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
title('OA: Coherence MS high vs low theta power trial of the controls','FontSize',12)
legend('high theta pow','low theta pow','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.36])
grid on


fname = strcat(dir_high_low_theta,sprintf('/MS_controls_OA_coherence_mean_high-low_pow_controls.jpg',cnt_m));
saveas(fig,fname);
