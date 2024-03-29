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
monkey = 'Maverick';
sess_list_idx = [16];
freq_range = 15:17;

filename_mod = ''; % -- loading file name for coherence averages ******************
filename_ctrl = ''; % -- loading file name for the list the coherences in sess_data_lfp_coherence
title_caption = 'S:CN - R:CN'; % -- title caption 
SR_brain_areas = 'CN_CN'; % -- name of SR brain area for the figures and coherence files 


recording = 'last_recording'; % -- folder where to load coherency files  *************

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


% -- Find indexes of common brain area
mod_areas = unique(mod.mod_areas);
area_idx = {}; % -- structure to store indexes 
cnt = 1;
for area = mod_areas
    area_idx{cnt} = find(ismember(mod.mod_areas,area));
    cnt =cnt +1;
end 

% -- sort in descending order the structure with indexes
[~,I] = sort(cellfun(@length,area_idx),'descend');
area_idx = area_idx(I);
idx_order = cell2mat(area_idx);
area_idx



fs = 1000;
fk = 200;
pad = 2;
N = 1;
W = 5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODULATOR - MODULATOR NETWORK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_mm = zeros(length(mod.mod_idx),length(mod.mod_idx));
for i = 1:length(mod.mod_idx)
    m1 = idx_order(i);
%     display(['----- m1 brain area: ',mod.mod_areas(m1)])
    ch1 = mod.mod_idx(m1); % channel modulator 1
    outliers_m1 = mod.outliers_E(m1).idx;

    for j = (i+1):length(mod.mod_idx)
        m2 = idx_order(j);
%         display(['m2 brain area: ',mod.mod_areas(m2)])
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
        
%         display(['Computing modulator-modulator coherence...'])
        [coh,f] = coherency(lfp_m1,lfp_m2,[N W],fs,fk,pad,0.05,1,1);
        
        c_mm(i,j) = mean(abs(coh(freq_range)));
        
    end
    
end


% FIGURE: modulator-modulator network
fig = figure;
imagesc(c_mm,[0,1]);
colorbar;
title(sprintf('%s: Theta Mod-Mod coherence, edge: %s',monkey,title_caption),'FontSize',10);

%%%%%%%%  BLOCKS LINES   %%%%%%%%
hold on;
x1 = 0.5 + length(area_idx{1});
x2 = x1;
y1 = 0;
y2 = x1;
line([x1,x2], [y1,y2], 'Color', 'w');

hx1 = x1; hx2 = x1; 
for i=2:length(area_idx)
    
    % horizontal lines 
    hx1 = hx2;
    hx2 = hx2 + length(area_idx{i});
    hy1 = hx1;
    hy2 = hx1;
    line([hx1,hx2], [hy1,hy2], 'Color', 'w');
    
    % vertical lines
    vx1 = hx2;
    vx2 = hx2;
    vy1 = hx1;
    vy2 = vy1 + length(area_idx{i});
    line([vx1,vx2], [vy1,vy2], 'Color', 'w');
    
end 


fname = strcat(dir_mod_network,'/mod_mod_coherence.png');
saveas(fig,fname)
writematrix(c_mm, strcat(dir_mod_network,'/c_mm.txt'),'delimiter',' ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODULATOR - CONTROL NETWORK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_mc = zeros(length(mod.mod_idx),length(mod.mod_idx));
for i = 1:length(mod.mod_idx)
    m1 = idx_order(i);
%     display(['----- m1 brain area: ',mod.mod_areas(m1)])

    ch1 = mod.mod_idx(m1); % channel modulator 1
    outliers_m1 = mod.outliers_E(m1).idx;

    display(['Computing modulator-control coherence...'])
    for j = 1:length(mod.mod_idx)
        m2 = idx_order(j);
%         display(['----- m2 brain area: ',mod.mod_areas(m2)])
        ch2 = mod.mod_idx(m2); % channel modulator 2
        
        BrainRegM2 = ctrl_SA.RecordPairMRIlabels{ch2,1}; % get brain area of modulator 2
        brain_idx = ctrl_SA.MRIlabels.(BrainRegM2).ElecIndx; % get all the indexes of recorded channels in the same brain area
        [ctrl_ch,pos]=intersect(brain_idx,ctrl_SA.ctrl_idx); % find common channels between controls and electrodes in this brain area
        
        for ch_c = ctrl_ch % for all the controls in the same brain areas as m2
            
            lfp_m1 = sq(mod.lfp_E(ch1,:,:));
            lfp_c2 = sq(ctrl_SA.lfp_E(ch_c,:,:)); % lfp control
            id = find(ctrl_SA.ctrl_idx == ch_c); % find index of the control in that specific brain area
            outliers_c2 = ctrl_SA.outliers_E(id).idx;
            
            % -- outliers
            outliers_mc = [outliers_m1, outliers_c2];
            outliers_mc = unique(outliers_mc);
            % -- remove outliers from sender and modulator
            lfp_m1(outliers_mc,:) = [];
            lfp_c2(outliers_mc,:) = [];
            
            [coh,f] = coherency(lfp_m1,lfp_c2,[N W],fs,fk,pad,0.05,1,1);
            
            c_mc(i,j) = c_mc(i,j) + mean(abs(coh(freq_range)));
        end 
        c_mc(i,j) = c_mc(i,j)/length(ctrl_ch); % normalize to get the average coherence
    end
    
end

% FIGURE: FULL MATRIX modulator-control network
fig = figure;
imagesc(c_mc,[0,1]);
title(sprintf('%s: Theta Mod-Ctrl coherence full, theta band, edge: %s',monkey,title_caption),'FontSize',10);
colorbar;

fname = strcat(dir_mod_network,'/mod_ctrl_coherence_full.png');
saveas(fig,fname)
writematrix(c_mc, strcat(dir_mod_network,'/c_mc_full.txt'),'delimiter',' ');



% Average over modulator 1 and modulator 2
c_mc_avg = zeros(length(mod.mod_idx),length(mod.mod_idx));
for m1=1:length(mod.mod_idx)
    for m2=(m1+1):length(mod.mod_idx)
        c_mc_avg(m1,m2) = (c_mc(m1,m2) + c_mc(m2,m1))/2
    end
end


% FIGURE: modulator-control network

fig = figure;
imagesc(c_mc_avg,[0,1]);
title(sprintf('%s: Theta Mod-Ctrl coherence, theta band, edge: %s',monkey,title_caption),'FontSize',10);
colorbar;


%%%%%%%%  BLOCKS LINES   %%%%%%%%
hold on;
x1 = 0.5 + length(area_idx{1});
x2 = x1;
y1 = 0;
y2 = x1;
line([x1,x2], [y1,y2], 'Color', 'w');

hx1 = x1; hx2 = x1; 
for i=2:length(area_idx)
    
    % horizontal lines 
    hx1 = hx2;
    hx2 = hx2 + length(area_idx{i});
    hy1 = hx1;
    hy2 = hx1;
    line([hx1,hx2], [hy1,hy2], 'Color', 'w');
    
    % vertical lines
    vx1 = hx2;
    vx2 = hx2;
    vy1 = hx1;
    vy2 = vy1 + length(area_idx{i});
    line([vx1,vx2], [vy1,vy2], 'Color', 'w');
    
end 

fname = strcat(dir_mod_network,'/mod_ctrl_coherence.png');
saveas(fig,fname)
writematrix(c_mc_avg, strcat(dir_mod_network,'/c_mc_avg.txt'),'delimiter',' ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HISTOGRAM OF THE NULL DISRTIBUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_mc_values = c_mc_avg(c_mc_avg ~= 0);
% -- FIGURE: histogram matrix c_mc_avg
fig = figure;
histogram(c_mc_values,20,'FaceAlpha',.6); grid on
legend('Mod-Ctrl coh distribution')
title(sprintf('%s - %s, Null model coh mod-ctrl',monkey,freq),'FontSize',12)
ylabel('mod-ctrl coherency')
xlabel('')

fname = strcat(dir_mod_network,'/null_distribution.png');
saveas(fig,fname)


% -- FIGURE: ATAN histogram matrix c_mc_avg
fig = figure;
histogram(atan(c_mc_values),20,'FaceAlpha',.6); grid on
legend('Mod-Ctrl coh distribution')
title(sprintf('%s - %s, Null model coh mod-ctrl',monkey,freq),'FontSize',12)
ylabel('mod-ctrl coherency')
xlabel('')

fname = strcat(dir_mod_network,'/null_distribution_atan.png');
saveas(fig,fname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PVALUE OF THE MOD-MOD NETWORK 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pval_mm = zeros(length(c_mm),length(c_mm));
for i=1:length(c_mm)
    for j=(i+1):length(c_mm)
        pval_mm(i,j) = nnz(c_mc_values > c_mm(i,j))/length(c_mc_values);   
    end    
end 

% Fill empty position with pval = 1
for j=1:length(c_mm)
    for i= j:length(c_mm)
        pval_mm(i,j) = 1;   
    end    
end 


% FIGURE: pvalue mod-mod network
fig = figure;
imagesc(pval_mm,[0,1]);
title(sprintf('%s: p-values, theta band, edge: %s',monkey,title_caption),'FontSize',10);
colormap(flipud(parula))
colorbar;

%%%%%%%%  BLOCKS LINES   %%%%%%%%
hold on;
x1 = 0.5 + length(area_idx{1});
x2 = x1;
y1 = 0;
y2 = x1;
line([x1,x2], [y1,y2], 'Color', 'w');

hx1 = x1; hx2 = x1; 
for i=2:length(area_idx)
    
    % horizontal lines 
    hx1 = hx2;
    hx2 = hx2 + length(area_idx{i});
    hy1 = hx1;
    hy2 = hx1;
    line([hx1,hx2], [hy1,hy2], 'Color', 'w');
    
    % vertical lines
    vx1 = hx2;
    vx2 = hx2;
    vy1 = hx1;
    vy2 = vy1 + length(area_idx{i});
    line([vx1,vx2], [vy1,vy2], 'Color', 'w');
    
end 


fname = strcat(dir_mod_network,'/pval_mod_mod_network.png');
saveas(fig,fname)
writematrix(c_mc_avg, strcat(dir_mod_network,'/pval_mod_mod_network.txt'),'delimiter',' ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Threshold p-values of the mod-mod network 
th = 0.1;
pval_mm_th = pval_mm;
pval_mm_th(pval_mm_th >= th) = 1;


% FIGURE: pvalue mod-mod network THRESHOLDED 
fig = figure;
imagesc(pval_mm_th,[0,1]);
title(sprintf('%s: th p-value = %.2f, theta band, edge: %s',monkey,th,title_caption),'FontSize',10);
colormap(flipud(parula))
colorbar;

%%%%%%%%  BLOCKS LINES   %%%%%%%%
hold on;
x1 = 0.5 + length(area_idx{1});
x2 = x1;
y1 = 0;
y2 = x1;
line([x1,x2], [y1,y2], 'Color', 'w');

hx1 = x1; hx2 = x1; 
for i=2:length(area_idx)
    
    % horizontal lines 
    hx1 = hx2;
    hx2 = hx2 + length(area_idx{i});
    hy1 = hx1;
    hy2 = hx1;
    line([hx1,hx2], [hy1,hy2], 'Color', 'w');
    
    % vertical lines
    vx1 = hx2;
    vx2 = hx2;
    vy1 = hx1;
    vy2 = vy1 + length(area_idx{i});
    line([vx1,vx2], [vy1,vy2], 'Color', 'w');
    
end 

fname = strcat(dir_mod_network,sprintf('/pval_mod_mod_network_p_th_%.2f.png',th));
saveas(fig,fname)










