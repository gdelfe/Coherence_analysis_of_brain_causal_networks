% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots the average coherence for the modulators and controls for
% only a specific sender-receiver edge
%
% @ Gino Del Ferraro, NYU, June 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)


addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes/Resting_State_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

sess_list_file = '/Sessions_with_modulator_info_movie.txt';
freq_band = 'theta_band';
freq = 'theta band';
monkey = 'Archie';
sess_list_idx = [7];

filename_mod = '_rec002'; % -- loading file name for coherence averages ******************
filename_ctrl = '_rec002'; % -- loading file name for the list the coherences in sess_data_lfp_coherence
title_caption = 'S:CN - R:ACC'; % -- title caption 
SR_brain_areas = 'CN_ACC'; % -- name of SR brain area for the figures and coherence files 


recording = 'rec002'; % -- folder where to load coherency files  *************

dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));
dir_mod_network = strcat(dir_RS,sprintf('/Modulators_network/%s',SR_brain_areas));
if ~exist(dir_mod_network, 'dir')
    mkdir(dir_mod_network)
end

fid = fopen(strcat(dir_RS,sess_list_file)); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

fk = 200; W = 5;
mod = [];
ctrl_SA = [];
ctrl_OA = [];

for i = sess_list_idx
    
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess)]);
    
    dir_Modulators = strcat(dir_RS,sprintf('/Sess_%d/Modulators/%s',Sess,recording));
    load(strcat(dir_Modulators,sprintf('/sess_data_lfp_coherence_fk_200_W_5%s.mat',filename_mod)));
    
    mod = sess_data_lfp;
    
    dir_ctrl_SA = strcat(dir_RS,sprintf('/Sess_%d/Controls_same_area/%s',Sess,recording));
    load(strcat(dir_ctrl_SA,sprintf('/sess_data_lfp_coherence_fk_200_W_5%s.mat',filename_ctrl)));
    
    ctrl_SA = sess_control_lfp; % stack up control same area structure for coherence
    
    
    dir_ctrl_OA = strcat(dir_RS,sprintf('/Sess_%d/Controls_other_areas/%s',Sess,recording));
    load(strcat(dir_ctrl_OA,sprintf('/sess_data_lfp_coherence_fk_200_W_5%s.mat',filename_ctrl)));
    
    ctrl_OA = sess_control_lfp;  % stack up control other areas structure for coherence
    
end


fs = 1000;
fk = 200;
pad = 2;
N = 1;
W = 5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODULATOR - MODULATOR NETWORK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_mm = zeros(length(mod.mod_idx),length(mod.mod_idx));
for m1 = 1:length(mod.mod_idx)
    ch1 = mod.mod_idx(m1); % channel modulator 1
    outliers_m1 = mod.outliers_E(m1).idx;

    for m2 = (m1+1):length(mod.mod_idx)
        ch2 = mod.mod_idx(m2); % channel modulator 2
        lfp_m2 = sq(mod.lfp_E(ch2,:,:));
        lfp_m1 = sq(mod.lfp_E(ch1,:,:));

        outliers_m2 = mod.outliers_E(m2).idx;
        
        % -- outliers
        outliers_mm = [outliers_m1, outliers_m2];
        outliers_mm = unique(outliers_mm);
        % -- remove outliers from sender and modulator
        lfp_m1(outliers_mm,:) = [];
        lfp_m2(outliers_mm,:) = [];
        
        display(['Computing modulator-modulator coherence...'])
        [coh,f] = coherency(lfp_m1,lfp_m2,[N W],fs,fk,pad,0.05,1,1);
        
        c_mm(m1,m2) = mean(abs(coh(15:17)));
        
    end
    
end

% FIGURE: modulator-modulator network
fig = figure;
tvimage(c_mm,'CLim',[0,1]);
title(sprintf('%s: Theta Mod-Mod coherence, theta band, edge: %s',monkey,title_caption),'FontSize',10);
colorbar;

fname = strcat(dir_mod_network,'/mod_mod_coherence.png');
saveas(fig,fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODULATOR - CONTROL NETWORK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_mc = zeros(length(mod.mod_idx),length(mod.mod_idx));
for m1 = 1:length(mod.mod_idx)
    ch1 = mod.mod_idx(m1); % channel modulator 1
    outliers_m1 = mod.outliers_E(m1).idx;

    for m2 = 1:length(mod.mod_idx)
        ch2 = mod.mod_idx(m2); % channel modulator 2
        
        BrainRegM2 = ctrl_SA.RecordPairMRIlabels{ch2,1}; % get brain area of modulator 2
        brain_idx = ctrl_SA.MRIlabels.(BrainRegM2).ElecIndx; % get all the indexes of channels in the same brain area
        [ctrl_ch,pos]=intersect(brain_idx,ctrl_SA.ctrl_idx); % find common channels between controls and electrodes in this brain area
        
        for ch_c = ctrl_ch % for all the controls in the same brain areas as m2
            
            lfp_m1 = sq(mod.lfp_E(ch1,:,:));
            lfp_c2 = sq(ctrl_SA.lfp_E(ch_c,:,:)); % lfp control
            i = find(ctrl_SA.ctrl_idx == ch_c); % find index of the control in that specific brain area
            outliers_c2 = ctrl_SA.outliers_E(i).idx;
            
            % -- outliers
            outliers_mc = [outliers_m1, outliers_c2];
            outliers_mc = unique(outliers_mc);
            % -- remove outliers from sender and modulator
            lfp_m1(outliers_mc,:) = [];
            lfp_c2(outliers_mc,:) = [];
            
            display(['Computing modulator-control coherence...'])
            [coh,f] = coherency(lfp_m1,lfp_c2,[N W],fs,fk,pad,0.05,1,1);
            
            c_mc(m1,m2) = c_mc(m1,m2) + mean(abs(coh(15:17)));
        end 
        c_mc(m1,m2) = c_mc(m1,m2)/length(ctrl_ch); % normalize to get the average coherence
    end
    
end

% FIGURE: FULL MATRIX modulator-control network
fig = figure;
tvimage(c_mc,'CLim',[0,1]);
title(sprintf('%s: Theta Mod-Ctrl coherence full, theta band, edge: %s',monkey,title_caption),'FontSize',10);
colorbar;

fname = strcat(dir_mod_network,'/mod_ctrl_coherence_full.png');
saveas(fig,fname)


% Average over modulator 1 and modulator 2
c_mc_avg = zeros(length(mod.mod_idx),length(mod.mod_idx));
for m1=1:length(mod.mod_idx)
    for m2=m1:length(mod.mod_idx)
        c_mc_avg(m1,m2) = (c_mc(m1,m2) + c_mc(m2,m1))/2
    end
end


% FIGURE: modulator-control network

fig = figure;
tvimage(c_mc_avg,'CLim',[0,1]);
title(sprintf('%s: Theta Mod-Ctrl coherence, theta band, edge: %s',monkey,title_caption),'FontSize',10);
colorbar;

fname = strcat(dir_mod_network,'/mod_ctrl_coherence.png');
saveas(fig,fname)






