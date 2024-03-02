
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code compute the control theta coherence between the sender-receiver
% for those trials having high/low theta power. The theta power is computed
% wrt control electrodes in other brain regions (all other regions but the
% modulators').
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

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

% ---- parameters for the coherence-gram
tot_time = 150001;
nt = tot_time;
fs = 1000;
fk = 200;
pad = 2;
N = 1;


coh_cr_high = [];
coh_cr_low = [];

for s = 1:size(sess_info{1},1)  % For each session with at least one modulator
    
    
    close all
    clear send_rec sess_data_lfp
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    
    load(strcat(dir_Modulators,name_struct_input)); % load sess_data_lfp, structure with session modulator info
    
    
    
    % number of modulators in that session
    n_mod = size(sess_data_lfp.mod_idx,2);
    
    lfp_R = sess_data_lfp.lfp_R; % receiver lfp
    outliers_R = sess_data_lfp.outliers_R;
    
    send_rec.sess_idx = sess_data_lfp.sess_idx;
    send_rec.mod_idx = sess_data_lfp.mod_idx;
    
    mod_Ch = sess_data_lfp.mod_idx;
    s_area = sess_data_lfp.sender_area;
    RecordPairMRIlabels = sess_data_lfp.RecordPairMRIlabels; % -- MRI labels of the recorder pars
    MRIlabels = sess_data_lfp.MRIlabels; % -- all the available MRI labels
    receiver_idx = sess_data_lfp.receiver_idx; % -- receiver idx
    
    cnt_m = 1; % counter for numb of modulators checked for the outlier trials remover
    for m = sess_data_lfp.mod_idx % for each modulator
        
        if m ~= receiver_idx % if modulator is not receiver
            
            display(['-- Modulator ',num2str(m),'  --- ',num2str(cnt_m),' out of ', num2str(n_mod)]);
            
            
            lfp_E = sq(sess_data_lfp.lfp_E(m,:,:)); % modulator lfp
            outliers_E = sess_data_lfp.outliers_E(cnt_m).idx; % trial with artifact for modulator
            % Compute the spectrum for each trial. Format: iTrial x times
            W = 3;
            [spec, f, err] = dmtspec(lfp_E,[1000/1e3,W],1e3,200);
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % high and low theta trial for the modulator
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
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
            
            
            
            
            % get control for the sender other area index
            ctrl_idx = choose_ALL_control_sender_other_Region_one_mod(RecordPairMRIlabels,MRIlabels,receiver_idx,s_area,mod_Ch);
            
            for ctrl = ctrl_idx % for all the controls associated with the sender
                
                
                lfp_C = sq(sess_data_lfp.lfp_E(ctrl,:,:)); % load control lfp
                lfp_C_ns = reshape(lfp_C,1,size(lfp_C,1)*size(lfp_C,2)); % RS control lfp as a unique time series
                
                % -- outliers Controls
                th_C = 4*std(lfp_C_ns,[],2); % -- threshold for LFP Sender
                max_C_split = max(abs(lfp_C),[],2);  % -- max of LFP for each time window
                outliers_C = find(max_C_split > th_C)'; % outliers controls
                
                outliers_CR = [outliers_C, outliers_R, outliers_E];
                outliers_CR = unique(outliers_CR);  % -- remove repeated entries in outliers
                
                high = setdiff(high_idx,outliers_CR); % remove trials with artifacts from high power trials
                low = setdiff(low_idx,outliers_CR);    % remove trials with artifacts from low power trials
                
                % compute coherence between sender-control and the receiver for high/low theta power trials
                W = 5;
                [c_cr_high,f] = coherency(lfp_R(high,:),lfp_C(high,:),[N W],fs,fk,pad,0.05,1,1);
                [c_cr_low,f] = coherency(lfp_R(low,:),lfp_C(low,:),[N W],fs,fk,pad,0.05,1,1);
                
                coh_cr_high = [coh_cr_high; c_cr_high];
                coh_cr_low = [coh_cr_low; c_cr_low];
                
                
            end % end controls sender
        end
        cnt_m = cnt_m + 1;
    end % end modulators
end % end sessions

% control-sender receiver coherence
ctrl_send_OA_coh.cr_high = coh_cr_high;
ctrl_send_OA_coh.cr_low = coh_cr_low;


keyboard;

save(strcat(dir_high_low_theta,'/coh_all_sess_controls_sender_OA-receiver.mat'),'ctrl_send_SA_coh')


load(strcat(dir_high_low_theta,'/coh_all_sess_controls_sender_OA-receiver.mat'))
% control-sender receiver coherence
coh_cr_high = ctrl_send_OA_coh.cr_high;
coh_cr_low = ctrl_send_OA_coh.cr_low;
f = linspace(0,200,409);


mean_all_coh_cr_high = mean(abs(coh_cr_high),1);
mean_all_coh_cr_low = mean(abs(coh_cr_low),1);


err_all_coh_cr_high = std(abs(coh_cr_high),0,1)/sqrt(size(coh_cr_high,1));
err_all_coh_cr_low = std(abs(coh_cr_low),0,1)/sqrt(size(coh_cr_low,1));


%%%%%%%%%%%%%%%%%%%%%%
% FIGURES %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%


% %%%% RECEIVER %%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureVisible','on')
%     load(strcat(dir_high_low_theta,sprintf('/Sess_%d/coherence_all.mat')));

fig = figure;
shadedErrorBar(f,mean_all_coh_cr_high,err_all_coh_cr_high,'lineprops',{'color',[28 199 139]/255 },'patchSaturation',0.5); hold on
shadedErrorBar(f,mean_all_coh_cr_low,err_all_coh_cr_low,'lineprops',{'color',[50 250 93]/255 },'patchSaturation',0.5);

grid on
% title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontName','Arial','FontSize',15);
ylabel('Coherence','FontName','Arial','FontSize',15);
title('Coherence MR high vs low theta power trial','FontSize',12)
legend('high theta pow','low theta pow','FontSize',10,'FontName','Arial')
set(gcf, 'Position',  [100, 600, 898, 500])
xlim([1 95])
% ylim([0 0.25])
grid on


fname = strcat(dir_high_low_theta,'/MR_controls_SA_coherence_mean.jpg');
saveas(fig,fname);





