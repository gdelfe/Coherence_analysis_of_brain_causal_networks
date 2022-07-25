
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code compute the theta coherence between the sender-receiver
% for those trials having high/low modulator theta power for the control
% electrode in the Same Area as the modulators 
%
%    @ Gino Del Ferraro, June 2022, Pesaran lab, NYU
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

name_struct_input = '/session_controls_same_area_lfp.mat'; % -- name file to load
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Maverick';

dir_high_low_theta = strcat(dir_main,sprintf('/%s/Resting_State/high_low_theta',monkey));
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

% ---- parameters for the coherence-gram
tot_time = 150001;
nt = tot_time;
fs = 1000;
fk = 200;
pad = 2;
N = 1;


coh_all_c_sr_high = [];
coh_all_c_sr_low = [];
coh_OC_c_sr_high = [];
coh_OC_c_sr_low = [];


for s = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    close all
    clear send_rec sess_data_lfp
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess_Ctrl = strcat(dir_RS_Theta,sprintf('/Sess_%d/Controls_same_area',Sess));
    
    load(strcat(dir_Sess_Ctrl,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed
    
    
    % number of controls in that session
    ctrl_Ch = sess_control_lfp.ctrl_idx % control channel indexes
    cnt_m = 1; % counter for numb of modulators checked for the outlier trials remover
    n_ctrl = length(ctrl_Ch);
    
    send_rec = [];
    
    for Ch = ctrl_Ch % for all the controls
        
        
        if Ch ~= sess_control_lfp.receiver_idx % if the electrode is not the receiver itself
            
            
            display(['-- Control ',num2str(Ch),'  --- ',num2str(cnt_m),' out of ', num2str(n_ctrl)]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ONLY COMMON trials for Modulator, Sender, Receiver
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            lfp_E = sq(sess_control_lfp.lfp_E(Ch,:,:)); % control LFP
            lfp_S = sess_control_lfp.lfp_S;
            lfp_R = sess_control_lfp.lfp_R;
            
            
            % outliers time series in sender and receiver
            outliers_S = sess_control_lfp.outliers_S;
            outliers_R = sess_control_lfp.outliers_R;
            outliers_E = sess_control_lfp.outliers_E(cnt_m).idx;
            
            
            outliers_ESR = [outliers_S, outliers_E, outliers_R];
            outliers_ESR = unique(outliers_ESR);  % -- remove repeated entries in outliers
            
            lfp_S(outliers_ESR,:) = [];
            lfp_R(outliers_ESR,:) = [];
            lfp_E(outliers_ESR,:) = [];
            
            
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
            
            send_rec.mod(cnt_m).ch = Ch;
            send_rec.mod(cnt_m).all_theta_pow = theta_pow;
            send_rec.mod(cnt_m).high_theta_pow = high_theta;
            send_rec.mod(cnt_m).high_pow_idx = high_idx;
            send_rec.mod(cnt_m).low_theta_pow = low_theta_pow;
            send_rec.mod(cnt_m).low_pow_idx = low_idx;
            send_rec.mod(cnt_m).lfp_S_clean = lfp_S;
            send_rec.mod(cnt_m).lfp_R_clean = lfp_R;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Coherency Modulator-Sender        %%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            W = 5;
            
            % -- coherence calculation via coherency()
            [c_sr_high,f] = coherency(lfp_S(high_idx,:),lfp_R(high_idx,:),[N W],fs,fk,pad,0.05,1,1);
            [c_sr_low,f] = coherency(lfp_S(low_idx,:),lfp_R(low_idx,:),[N W],fs,fk,pad,0.05,1,1);
            
            
            send_rec.mod(cnt_m).c_sr_high = c_sr_high;
            send_rec.mod(cnt_m).c_sr_low = c_sr_low;
            send_rec.freq = f;
            
            
            % store the value for the c_sr for each modulator high/low
            coh_all_c_sr_high = [coh_all_c_sr_high; c_sr_high];
            coh_all_c_sr_low = [coh_all_c_sr_low; c_sr_low];
                     
        end
        cnt_m = cnt_m +1;
        
    end % end of for cycle for all the controls in a given session
    
    if ~isempty(send_rec)
        dir_Sess_send_rec_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/send_rec/Data',Sess));    
        save(strcat(dir_Sess_send_rec_data,'/send_rec_coh_for_session_controls_SA.mat'),'send_rec');
    end
    
end


save(strcat(dir_high_low_theta,'/coh_all_sess_sr_high_controls_SA.mat'),'coh_all_c_sr_high')
save(strcat(dir_high_low_theta,'/coh_all_sess_sr_low_controls_SA.mat'),'coh_all_c_sr_low')


% load(strcat(dir_high_low_theta,'/coh_all_sess_sr_high_controls_SA.mat'))


mean_all_coh_sr_high = mean(abs(coh_all_c_sr_high),1);
mean_all_coh_sr_low = mean(abs(coh_all_c_sr_low),1);


err_all_coh_sr_high = std(abs(coh_all_c_sr_high),0,1)/sqrt(size(coh_all_c_sr_high,1));
err_all_coh_sr_low = std(abs(coh_all_c_sr_low),0,1)/sqrt(size(coh_all_c_sr_low,1));



%%%%%%%%%%%%%%%%%%%%%%
% FIGURES %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%


% %%%% SENDER-RECEIVER %%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureVisible','on')
%     load(strcat(dir_high_low_theta,sprintf('/Sess_%d/coherence_all.mat')));

fig = figure;
shadedErrorBar(f,mean_all_coh_sr_high,err_all_coh_sr_high,'lineprops',{'color',[0, 102, 204]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_sr_low,err_all_coh_sr_low,'lineprops',{'color',[128, 191, 255]/255 },'patchSaturation',0.5);

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
title('Coherence SR high vs low theta power trial (control SA)','FontSize',12)
legend('high theta pow','low theta pow','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.25])
grid on


fname = strcat(dir_high_low_theta,sprintf('/SR_all_coherence_mean_controls_SA.jpg',cnt_m));
saveas(fig,fname);


