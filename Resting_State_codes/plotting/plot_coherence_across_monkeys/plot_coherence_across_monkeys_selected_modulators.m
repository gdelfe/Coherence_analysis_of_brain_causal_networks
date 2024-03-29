

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';


freq_band = 'theta_band';

N = 10 % --- max number of modulators 
N_mav = 101; 
N_arc = 44;
titleN = sprintf('%d perc top modulators - last_rec - rec001-002',N);

recording_mav = 'last_recording';
recording_arc = 'rec001_002_no_bad_sessions';

namef_mav = '';
namef_arc = '_rec001_002_no_bad_sess';

name = '_last_rec-rec001_002_no_bad_sess';
namefig = sprintf('%s.fig',name);
namepng = sprintf('%s.png',name);

recording_both = 'last_rec-rec001_002';
decoding = 'AUC';


dir_list_mav = strcat(dir_main,sprintf('Maverick/Resting_state/%s/Modulators_controls',freq_band));
dir_list_arc = strcat(dir_main,sprintf('Archie/Resting_state/%s/Modulators_controls',freq_band));

dir_Maverick = strcat(dir_main,sprintf('Maverick/Resting_State/%s/Modulators_Controls_avg_results/%s',freq_band));
dir_Archie = strcat(dir_main,sprintf('Archie/Resting_State/%s/Modulators_Controls_avg_results/%s',freq_band));

dir_Maverick_avg = strcat(dir_main,sprintf('Maverick/Resting_State/%s/Modulators_Controls_avg_results/%s',freq_band,recording_mav));
dir_Archie_avg = strcat(dir_main,sprintf('Archie/Resting_State/%s/Modulators_Controls_avg_results/%s',freq_band,recording_arc));

fk = 200; W = 5;
dir_both_monkeys = strcat(dir_main,sprintf('both_monkeys/%s/modulators_vs_controls/%s/%s',freq_band,recording_both,decoding));




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       MAVERICK              %
mod_list = importdata(strcat(dir_list_mav,'/modulators_sorted_decod_accuracy.txt'));
fid = fopen(strcat(dir_Maverick,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


% -- select the first N index
mod_idx = mod_list(1:N_mav,4);
sess_numb_M = unique(mod_list(1:N_mav,1));

% -- find the session index corresponding to the session with top modulators 
sess_idx_M = [];
for i=1:length(sess_numb_M)
    sess_idx_M = [sess_idx_M, find(sess_info{1}==sess_numb_M(i))];
end


% %%%%%%%%% MODULATORS Maverick %%%%%%
load(strcat(dir_Maverick_avg,sprintf('/coh_spec_m_fk_%d_W_%d%s.mat',fk,W,namef_mav))); % structure mod
load(strcat(dir_Maverick_avg,sprintf('/coh_spec_sr_fk_%d_W_%d%s.mat',fk,W,namef_mav))); % structure stim
struct_mod_mav = mod;
struct_stim_mav = stim;

% -- get coherences for the specific modulators and specific sessions 
struct_mod_mav = struct_mod_mav(mod_idx); 
struct_stim_mav = struct_stim_mav(sess_idx_M);




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ARCHIE                %
mod_list = importdata(strcat(dir_list_arc,'/modulators_sorted_decod_accuracy.txt'));
fid = fopen(strcat(dir_Archie,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


% -- select the first N index
mod_idx = mod_list(1:N_arc,4);
sess_numb_A = unique(mod_list(1:N_arc,1));

% -- find the session index corresponding to the session with top modulators 
sess_idx_A = [];
for i=1:length(sess_numb_A)
    sess_idx_A = [sess_idx_A, find(sess_info{1}==sess_numb_A(i))];
end


% %%%%%%%%% MODULATORS Archie %%%%%%
load(strcat(dir_Archie_avg,sprintf('/coh_spec_m_fk_%d_W_%d%s.mat',fk,W,namef_arc))); % structure mod
load(strcat(dir_Archie_avg,sprintf('/coh_spec_sr_fk_%d_W_%d%s.mat',fk,W,namef_arc))); % structure stim
struct_mod_arc = mod;
struct_stim_arc = stim;

struct_mod_arc = struct_mod_arc(mod_idx);
struct_stim_arc = struct_stim_arc(sess_idx_A);


% %%%%%%%%% Computes modulators for Maverick and Archie %%%%%%
modulators = mean_coh_and_spec_across_monkeys(struct_mod_mav,struct_stim_mav,struct_mod_arc,struct_stim_arc);






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLS SAME AREA MAVERICK              %

load(strcat(dir_Maverick_avg,sprintf('/coh_spec_m_Controls_same_area_fk_%d_W_%d%s',fk,W,namef_mav)));
load(strcat(dir_Maverick_avg,sprintf('/coh_spec_sr_Controls_same_area_fk_%d_W_%d%s',fk,W,namef_mav)));
mod_ctrl_SA_mav = mod;
stim_ctrl_SA_mav = stim;

ctrl_list = importdata(strcat(dir_list_mav,'/control_list_same_area.txt'));

% -- get the index of the controls associated with sessions where the top modulators are found 
ctrl_idx_M = [];
for i = sess_numb_M'
    ctrl_idx_M = [ctrl_idx_M; find(ctrl_list(:,1)== i)];
end

mod_ctrl_SA_mav = mod_ctrl_SA_mav(ctrl_idx_M);
stim_ctrl_SA_mav = stim_ctrl_SA_mav(sess_idx_M);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLS SAME AREA ARCHIE              %

load(strcat(dir_Archie_avg,sprintf('/coh_spec_m_Controls_same_area_fk_%d_W_%d%s',fk,W,namef_arc)));
load(strcat(dir_Archie_avg,sprintf('/coh_spec_sr_Controls_same_area_fk_%d_W_%d%s',fk,W,namef_arc)));
mod_ctrl_SA_arc = mod;
stim_ctrl_SA_arc = stim;

ctrl_list = importdata(strcat(dir_list_arc,'/control_list_same_area.txt'));

% -- get the index of the controls associated with sessions where the top modulators are found 
ctrl_idx_A = [];
for i = sess_numb_A'
    ctrl_idx_A = [ctrl_idx_A; find(ctrl_list(:,1)== i)];
end


mod_ctrl_SA_arc = mod_ctrl_SA_arc(ctrl_idx_A);
stim_ctrl_SA_arc = stim_ctrl_SA_arc(sess_idx_A);


ctrl_SA = mean_coh_and_spec_across_monkeys(mod_ctrl_SA_mav,stim_ctrl_SA_mav,mod_ctrl_SA_arc,stim_ctrl_SA_arc);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLS OTHER AREAS MAVERICK              %

load(strcat(dir_Maverick_avg,sprintf('/coh_spec_m_Controls_other_areas_fk_%d_W_%d%s',fk,W,namef_mav)));
load(strcat(dir_Maverick_avg,sprintf('/coh_spec_sr_Controls_other_areas_fk_%d_W_%d%s',fk,W,namef_mav)));
mod_ctrl_OA_mav = mod;
stim_ctrl_OA_mav = stim;

ctrl_list = importdata(strcat(dir_list_mav,'/control_list_other_areas.txt'));

% -- get the index of the controls associated with sessions where the top modulators are found 
ctrl_idx_M = [];
for i = sess_numb_M'
    ctrl_idx_M = [ctrl_idx_M; find(ctrl_list(:,1)== i)];
end

mod_ctrl_OA_mav = mod_ctrl_OA_mav(ctrl_idx_M);
stim_ctrl_OA_mav = stim_ctrl_OA_mav(sess_idx_M);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTROLS OTHER AREAS ARCHIE              %

load(strcat(dir_Archie_avg,sprintf('/coh_spec_m_Controls_other_areas_fk_%d_W_%d%s',fk,W,namef_arc)));
load(strcat(dir_Archie_avg,sprintf('/coh_spec_sr_Controls_other_areas_fk_%d_W_%d%s',fk,W,namef_arc)));
mod_ctrl_OA_arc = mod;
stim_ctrl_OA_arc = stim;

ctrl_list = importdata(strcat(dir_list_arc,'/control_list_other_areas.txt'));

% -- get the index of the controls associated with sessions where the top modulators are found 
ctrl_idx_A = [];
for i = sess_numb_A'
    ctrl_idx_A = [ctrl_idx_A; find(ctrl_list(:,1)== i)];
end


mod_ctrl_OA_arc = mod_ctrl_OA_arc(ctrl_idx_A);
stim_ctrl_OA_arc = stim_ctrl_OA_arc(sess_idx_A);


ctrl_OA = mean_coh_and_spec_across_monkeys(mod_ctrl_OA_mav,stim_ctrl_OA_mav,mod_ctrl_OA_arc,stim_ctrl_OA_arc);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = linspace(1,fk,size(modulators.mean_coh_ms,2)); % frequency values (range)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           FIGURES  COHERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% COHERENCES MODULATORS vs CONTROLS %%%%%%%%%%%%

 
set(0,'DefaultFigureVisible','on')
% -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
fig = figure;
hold all


shadedErrorBar(f,modulators.mean_coh_mr,modulators.err_mr,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA.mean_coh_mr,ctrl_SA.err_mr,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA.mean_coh_mr,ctrl_OA.err_mr,'lineprops',{'color',[102, 255, 217]/255},'patchSaturation',0.4); hold on

grid on
title(sprintf('Both animals: Abs MR coherence, %s - Resting State',titleN),'FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Receivers','Controls-Receivers  same area','Controls-Receiver  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_both_monkeys,sprintf('/coherency_MR_mod_vs_ctrl_both_monkeys_N_%d_W_%d_fk_%d%s',N,W,fk,namefig));
saveas(fig,fname)
fname = strcat(dir_both_monkeys,sprintf('/coherency_MR_mod_vs_ctrl_both_monkeys_N_%d_W_%d_fk_%d%s',N,W,fk,namepng));
saveas(fig,fname)

% --- ELECTRODE-SENDER coherence   -------%

fig = figure;
hold all


shadedErrorBar(f,modulators.mean_coh_ms,modulators.err_ms,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_SA.mean_coh_ms,ctrl_SA.err_ms,'lineprops',{'color',[255, 51, 153]/255},'patchSaturation',0.4); hold on
shadedErrorBar(f,ctrl_OA.mean_coh_ms,ctrl_OA.err_ms,'lineprops',{'color',[255, 128, 128]/255},'patchSaturation',0.4); hold on

grid on
title(sprintf('Both animals: Abs MS coherence, %s - Resting State',titleN),'FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('Modulators-Senders','Controls-Senders  same area','Controls-Senders  other areas','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
grid on

fname = strcat(dir_both_monkeys,sprintf('/coherency_MS_mod_vs_ctrl_both_monkeys_N_%d_W_%d_fk_%d%s',N,W,fk,namefig));
saveas(fig,fname)
fname = strcat(dir_both_monkeys,sprintf('/coherency_MS_mod_vs_ctrl_both_monkeys_N_%d_W_%d_fk_%d%s',N,W,fk,namepng));
saveas(fig,fname)

% keyboard 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      FIGURES  SPECTRUMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% SPECTRUMS MODULATORS vs RECEIVERS vs CONTROLS %%%%%%%%%%%%

 
% set(0,'DefaultFigureVisible','on')
% % -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
% fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
% hold all
% 
% 
% shadedErrorBar(f,modulators.mean_spec_m,modulators.err_S_m,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,modulators.mean_spec_r,modulators.err_S_r,'lineprops',{'color',[26 198 1]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,modulators.mean_spec_s,modulators.err_S_s,'lineprops',{'color',[102, 255, 217]/255},'patchSaturation',0.4); hold on
% 
% grid on
% title('RS: Modulators, Senders, Receivers spectrum','FontSize',11);
% xlabel('freq (Hz)');
% ylabel('spectrum');
% legend('Modulators','Receiver','Sender','FontSize',10)
% set(gcf, 'Position',  [100, 600, 1000, 600])
% grid on
% 
% fname = strcat(dir_Controls,sprintf('/spectrum_modulators_sender_receiver_W_%d_fk_%d.png',W,fk));
% saveas(fig,fname)
% fname = strcat(dir_Controls,sprintf('/spectrum_modulators_sender_receiver_W_%d_fk_%d.fig',W,fk));
% saveas(fig,fname)
% 
% 
% figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
% hold all
% for i=1:18
%     
% plot(f,abs(stim_mod(i).s_s))
% hold on
% end
% hold on 
% shadedErrorBar(f,modulators.mean_spec_s,modulators.err_S_s,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
% grid on
% 
% 
% set(0,'DefaultFigureVisible','on')
% % -- FIGURE: Plot average coherence across sessions for MR, CR same area, CR other areas
% fig = figure;
% hAx=axes;
% hAx.XScale='linear'
% hAx.YScale='log'
% hold all
% 
% 
% shadedErrorBar(f,modulators.mean_spec_m,modulators.err_S_m,'lineprops',{'color',[0, 153, 255]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,modulators.mean_spec_r,modulators.err_S_r,'lineprops',{'color',[0255, 102, 0]/255},'patchSaturation',0.4); hold on
% shadedErrorBar(f,modulators.mean_spec_s,modulators.err_S_s,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
% 
% grid on
% title('RS: Modulators, Senders, Receivers spectrum','FontSize',11);
% xlabel('freq (Hz)');
% ylabel('spectrum');
% legend('Modulators','Receiver','Sender','FontSize',10)
% set(gcf, 'Position',  [100, 600, 1000, 600])
% grid on
% 
% fname = strcat(dir_RS,sprintf('/spectrum_RS_modulators_sender_receiver_W_%d_fk_%d.png',W,fk));
% saveas(fig,fname)

