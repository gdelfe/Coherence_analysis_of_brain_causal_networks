
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the Resting State coherence between all the pairs in a
% given session. We call it "Grand Coherence"
%
%    @ Gino Del Ferraro, June 2022, Pesaran lab, NYU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

name_struct_input = '/session_data_lfp.mat';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

fs = 1000;
fk = 200;
pad = 2;
N = 1;
W = 5;
    
    
for s = 24 %4:length(sess_info{1})
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    
    load(strcat(dir_Modulators,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed
    
    lfp_E = sess_data_lfp.lfp_E;
    lfp_S = sess_data_lfp.lfp_S;
    S(1,:,:) = lfp_S;
    
    lfp_X = cat(1,lfp_E,S);
    
  
    coh = [];
    for i = 1:size(lfp_X,1)
        display(['Computing sender-receiver coherence...  -- ',num2str(i)])
        for j=(i+1):size(lfp_X,1)
            
            % -- coherence calculation via coherency()
            [c_el,f] = coherency(sq(lfp_X(i,:,:)),sq(lfp_X(j,:,:)),[N W],fs,fk,pad,0.05,1,1);
            coh = [coh; c_el];
        end
    end
    
    dir_Mod_results = strcat(dir_Modulators,'/grand_coherence');
    if ~exist(dir_Mod_results, 'dir')
        mkdir(dir_Mod_results)
    end
    
    save(strcat(dir_Mod_results,sprintf('/grand_coherence_sess_%d.mat',Sess)),'coh');
    
    %%%%%%%%%%%%%%%%%%%%%%
    % MEAN COHERENCE
    
    mean_coh = mean(abs(coh));
    % --- Error bars
    err = std(abs(coh))/sqrt(size(coh,1));
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % FIGURE %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    fig = figure;
    hold all
    
    shadedErrorBar(f,mean_coh,err,'lineprops',{'color',[0.4940, 0.1840, 0.5560]},'patchSaturation',0.4); hold on
    
    grid on
    title(sprintf('Pair-wise mean coherence for Sess %d - Maverick',Sess),'FontSize',10);
    xlabel('freq (Hz)');
    ylabel('coherence');
    
    legend('coherence','FontSize',10)
    set(gcf, 'Position',  [100, 600, 1000, 600])
    xlim([0 95])
    
    fname = strcat(dir_Mod_results,sprintf('/grand_coherence_sess_%d.jpg',Sess));
    saveas(fig,fname)
    fname = strcat(dir_Mod_results,sprintf('/grand_coherence_sess_%d.fig',Sess));
    saveas(fig,fname)
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SENDER-CONTROL COHERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RecordPairMRIlabels = sess_data_lfp.RecordPairMRIlabels;
MRIlabels = sess_data_lfp.MRIlabels;
sender_area = sess_data_lfp.sender_area;
receiver_idx = sess_data_lfp.receiver_idx;
[ctrl_el,ctrl_Reg] = receiver_control_other_region(RecordPairMRIlabels,MRIlabels,receiver_idx,sender_area);

[ctrl_el_NO_SR,ctrl_Reg_NO_SR] = receiver_control_no_send_no_rec(RecordPairMRIlabels,MRIlabels,receiver_idx,sender_area)


lfp_E = sess_data_lfp.lfp_E;

% Sender-Control Coherence with Sender and Receiver area excluded 
coh_sc = [];
for el = ctrl_el
    
    [c_sc,f] = coherency(lfp_S,sq(lfp_E(el,:,:)),[N W],fs,fk,pad,0.05,1,1);
    coh_sc = [coh_sc; c_sc];
    
end
mean_sc = mean(abs(coh_sc));
err_sc = std(abs(coh_sc))/sqrt(length(ctrl_el));


% Sender-Control Coherence with Sender and Receiver electrodes excluded 
coh_sc_NO_sr = [];
for el = ctrl_el_NO_SR
    
    [c_sc,f] = coherency(lfp_S,sq(lfp_E(el,:,:)),[N W],fs,fk,pad,0.05,1,1);
    coh_sc_NO_sr = [coh_sc_NO_sr; c_sc];
    
end
mean_sc_NO_sr = mean(abs(coh_sc_NO_sr));
err_sc_NO_sr = std(abs(coh_sc_NO_sr))/sqrt(length(ctrl_el_NO_SR));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SENDER-RECEIVER COHERENCE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lfp_S = sess_data_lfp.lfp_S;
lfp_R = sess_data_lfp.lfp_R;

% coherence computed using all the trials - errobar: Jackknife method 
[c_sr,f,S_s,S_r,coh_err] = coherency(lfp_S,lfp_R,[N W],fs,fk,pad,0.05,1,1);
% compute coherence for each trial 
[c_sr_trials,f] = coherency(lfp_S,lfp_R,[N W],fs,fk,pad,0.05,0,1);


load(strcat(dir_Mod_results,sprintf('/grand_coherence_sess_%d.mat',Sess)));

%%%%%%%%%%%%%%%%%%%%%%
% MEAN COHERENCE

% Grand coherence mean 
mean_coh = mean(abs(coh));
err = std(abs(coh))/sqrt(size(coh,1));

% Sender-Receiver coherence averaged across trials
mean_coh_sr = mean(abs(c_sr_trials));
err_sr = std(abs(c_sr_trials))/sqrt(size(c_sr_trials,1));

err_jackknife =  (coh_err(1,:)-coh_err(2,:))/2; % (top-bottom)/2

%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

fig = figure;
hold all

% sender-receiver coherence line 
% plot(f,abs(c_sr),'color',[0.4940, 0.1840, 0.5560]); hold on

% sender-receiver coherence with Jackknife err_bar
shadedErrorBar(f,abs(c_sr),err_jackknife,'lineprops',{'color',[220, 132, 37]/255},'patchSaturation',0.4); hold on

% grand coherence
shadedErrorBar(f,mean_coh,err,'lineprops',{'color',[0 105 143]/255},'patchSaturation',0.4); hold on

% sender-receiver coherence as average across all the trials
% shadedErrorBar(f,mean_coh_sr,err_sr,'lineprops',{'color',[51, 153, 255]/255},'patchSaturation',0.4); hold on

% sender-control coherence (ctrl in all other areas)
% shadedErrorBar(f,mean_sc,err_sc,'lineprops',{'color',[255, 153, 51]/255},'patchSaturation',0.4); hold on

% sender-control coherence with only sender-receiver excluded 
shadedErrorBar(f,mean_sc_NO_sr,err_sc_NO_sr,'lineprops',{'color',[60 219 216]/255},'patchSaturation',0.4); hold on
xlabel('freq (Hz)');
ylabel('coherence');

% legend('sender-receiver','grand coherence','SR Jackknife','SR average across trials','S-controls other Reg','FontSize',10)
legend('sender-receiver','grand coherence','S-Controls all electrodes','FontSize',10)

set(gcf, 'Position',  [100, 600, 800, 500])
xlim([0 95])
ylim([0 0.32])
grid on

fname = strcat(dir_Mod_results,sprintf('/SR_and_grand_coherence_and_Sender_Controls_sess_%d.pdf',Sess));
saveas(fig,fname)
fname = strcat(dir_Mod_results,sprintf('/SR_and_grand_coherence_and_Sender_Controls_sess_%d.fig',Sess));
saveas(fig,fname)





figure;
plot(f,abs(c_sr_trials))

grid on
title(sprintf('Pair-wise mean coherence for Sess %d - Maverick',Sess),'FontSize',10);
xlabel('freq (Hz)');
ylabel('coherence');

legend('coherence','FontSize',10)
set(gcf, 'Position',  [100, 600, 1000, 600])
xlim([0 95])

fname = strcat(dir_Mod_results,sprintf('/grand_coherence_sess_%d.jpg',Sess));
saveas(fig,fname)
fname = strcat(dir_Mod_results,sprintf('/grand_coherence_sess_%d.fig',Sess));
saveas(fig,fname)














