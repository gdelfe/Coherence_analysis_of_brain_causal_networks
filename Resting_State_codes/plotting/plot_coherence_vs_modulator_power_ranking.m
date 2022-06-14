
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes the theta-coherence as a function of the
% modulator score, where the modulator score is computed -in this case- as
% the modulator theta-power in the [4,8] Hz range
%
%    @ Gino Del Ferraro, May 2022, Pesaran lab, NYU
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
filename = '.mat'; % -- filename for sess_data_info.mat
recording = 'movie';
freq_band = 'theta_band';
dir_both_monkeys = strcat(dir_main,sprintf('/both_monkeys/%s/modulators_vs_controls/%s',freq_band,recording));


N_list = [10,20,30,40,50,100];  % percentages of modulators 
N_mav_list = [10,20,30,40,50,101];  % modulators maverick
N_arc_list = [4,8,12,17,22,44];     % modulators archie 


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

score_sort_mav = sortrows(score_mat_mav,1,'descend') % theta score / mod / Sess / cnt_e


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONKEY ARCHIE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

monkey = 'Archie';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

% % -- exclude bad sessions
% sess_list = 1:size(sess_info{1},1);
% % excluded_sess = [8,22,30,31]; % archie
% excluded_idx = [2,5,8,9]; % archie
% sess_list(excluded_idx) = [];

cnt_e = 1;
score_mat_arc = [];
for s = sess_list  % For each session with at least one modulator
    
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

score_sort_arc = sortrows(score_mat_arc,1,'descend') % theta score / mod / Sess / cnt_e




N_mav = N_mav_list(end);
N_arc = N_arc_list(end);
c_mr_mav = zeros(N_mav,409);
c_mr_arc = zeros(N_arc,409);
c_ms_mav = zeros(N_mav,409);
c_ms_arc = zeros(N_arc,409);

%%%%%%%%%%%%%%%%%%%%%%%
% MAVERICK
%%%%%%%%%%%%%%%%%%%%%%

monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));

display(['Maverick  # mod = ',num2str(N_mav)]);
for m_idx = 1:N_mav
    
    Sess = score_sort_mav(m_idx,3);
    dir_Mod_recording = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/%s',Sess,recording));
    % %%%%%%%%% MODULATORS  %%%%%%
    load(strcat(dir_Mod_recording,sprintf('/sess_data_lfp_coherence_fk_200_W_5_%s.mat',recording))); % structure mod
    mav = sess_data_lfp;
    m = score_sort_mav(m_idx,2); % modulator
    cnt_m = find(mav.mod_idx == m); % modulator cnt
    
    mav_c_mr = mav.mod(cnt_m).c_mr; % modulator-receiver coherence
    c_mr_mav(m_idx,:) = mav_c_mr;
    mav_c_ms = mav.mod(cnt_m).c_ms; % modulator-sender coherence
    c_ms_mav(m_idx,:) = mav_c_ms;
    
end

%%%%%%%%%%%%%%%%%%%%%%%
% ARCHIE
%%%%%%%%%%%%%%%%%%%%%%

monkey = 'Archie';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));

display(['Archie  # mod = ',num2str(N_arc)]);
for m_idx = 1:N_arc
    
    Sess = score_sort_arc(m_idx,3);
    dir_Mod_recording = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/%s',Sess,recording));
    % %%%%%%%%% MODULATORS  %%%%%%
    load(strcat(dir_Mod_recording,sprintf('/sess_data_lfp_coherence_fk_200_W_5_%s.mat',recording))); % structure mod
    arc = sess_data_lfp;
    m = score_sort_arc(m_idx,2); % modulator
    cnt_m = find(arc.mod_idx == m); % modulator cnt
    
    arc_c_mr = arc.mod(cnt_m).c_mr; % modulator-receiver coherence
    c_mr_arc(m_idx,:) = arc_c_mr;
    arc_c_ms = arc.mod(cnt_m).c_ms; % modulator-sender coherence
    c_ms_arc(m_idx,:) = arc_c_ms;
    
end


c_mr_mean = [];
err_c_mr = [];
c_ms_mean = [];
err_c_ms = [];
for i=1:length(N_list)
    
    close all
    N = N_list(i); % --- max number of modulators
    N_mav = N_mav_list(i);
    N_arc = N_arc_list(i);
   
    
    % Modulator-Receiver 
    c_mr = [c_mr_mav(1:N_mav,:); c_mr_arc(1:N_arc,:)]; % coherence maverick + archie 
    c_mr_mean = [c_mr_mean, mean(mean(abs(c_mr(:,9:18)),1))];% average coherence across theta-band and across modulators 
    std_c_mr = std(abs(c_mr(:,9:18)),1)/sqrt(size(c_mr,1)); % SEM wrt the theta range
    err_c_mr = [err_c_mr, sqrt(sum(std_c_mr.^2))/length(std_c_mr)] % MSE
    
    % Modulator-Sender
    c_ms = [c_ms_mav(1:N_mav,:); c_ms_arc(1:N_arc,:)]; % coherence maverick + archie 
    c_ms_mean = [c_ms_mean, mean(mean(abs(c_ms(:,9:18)),1))];% average coherence across theta-band and across modulators 
    std_c_ms = std(abs(c_ms(:,9:18)),1)/sqrt(size(c_ms,1)); % SEM wrt the theta range
    err_c_ms = [err_c_ms, sqrt(sum(std_c_ms.^2))/length(std_c_ms)] % MSE
    
end




fig = figure;
errorbar(N_list,c_mr_mean,err_c_mr,'color',[28 199 139]/255,'LineWidth',1); hold on % modulator-receiver 
errorbar(N_list,c_ms_mean,err_c_ms,'color',[0.4940, 0.1840, 0.5560],'LineWidth',1); % modulator-sender 

xlim([0,110])
ylim([0.05 0.57])
xticks([0 10 20 30 40 50 60 80 100])
grid on
set(gca,'FontSize',14)
title('Theta-Coherence vs theta-power ranking ','FontSize',11')
xlabel('Modulator theta-power ranking','FontName','Arial','FontSize',15);
ylabel('Theta-Coherence','FontName','Arial','FontSize',15);
legend('Modulator - Receiver','Modulator-Sender','FontSize',10)
set(gcf, 'Position',  [100, 600, 650, 450])

fname = strcat(dir_both_monkeys,'/theta-coh_vs_theta-power-ranking.pdf');
saveas(fig,fname);
fname = strcat(dir_both_monkeys,'/theta-coh_vs_theta-power-ranking.png');
saveas(fig,fname);



