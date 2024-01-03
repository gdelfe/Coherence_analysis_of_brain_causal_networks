
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes');
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data';
dir_high_low_theta = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/Maverick/Resting_State/high_low_theta';

name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % -- filename for sess_data_info.mat

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

n_iter = 5000;

% ---- parameters for the coherence-gram
tot_time = 150001;
nt = tot_time;
fs = 1000;
fk = 200;
pad = 2;
N = 1;
W = 5;


Sess = sess_info{1}(1); % Session number
display(['-- Session ',num2str(1),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));

load(strcat(dir_Modulators,name_struct_input)); % load sess_data_lfp, structure with session modulator info

session(1).count = 1;
session(1).receivers = sess_data_lfp.receiver_idx;
session(1).sender_pair = sess_data_lfp.sender_pair;
session(1).labels = 1;

% %%%%%%%%%%%%%%%%%%%%
% Create structure to count how many session have the same sender but
% different receiver. Receivers must be taken out when doing the permutaion
% test of SR -coherence vs grand_coherence

i = 1;
sess_info{2}
for s = 2:size(sess_info{2},1)
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    load(strcat(dir_Modulators,name_struct_input)); % load sess_data_lfp, structure with session modulator info


    if isequal(sess_info{2}(s),sess_info{2}(s-1)) % if it is the same session (different receiver, but same sender)
        session(i).count = session(i).count + 1; % increase the count 
        session(i).receivers = [session(i).receivers, sess_data_lfp.receiver_idx];
        session(i).sender_pair = [session(i).sender_pair, sess_data_lfp.sender_pair];
        session(i).labels = [session(i).labels, s]; 
        
    elseif not(isequal(sess_info{2}(s),sess_info{2}(s-1))) % if it is a new session (different sender)
        i = i +1;
        session(i).count = 1;
        session(i).receivers = sess_data_lfp.receiver_idx;
        session(i).sender_pair = sess_data_lfp.sender_pair;
        session(i).labels = s;
        
    end
    session(i)
    
end



SR_coh = [];
grand_coh = [];
for s = 1:size(session,2)  % For each session (with one sender)
         
    Sess = sess_info{1}(session(s).labels(1)); % Session number 
    display(['Session ',num2str(Sess)])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    dir_Mod_results = strcat(dir_Modulators,'/grand_coherence');
    
    load(strcat(dir_Modulators,name_struct_input)); % load sess_data_lfp, structure with session modulator info 
    load(strcat(dir_Mod_results,sprintf('/grand_coherence_sess_%d.mat',Sess)));
    % mean grand coherence 
    mean_coh = mean(abs(coh));
    grand_coh = [grand_coh; mean_coh];
    
    for i=1:session(s).count % for all the sessions with the same sender 
        
        lfp_S = sess_data_lfp.lfp_S; % sender 
        r = session(s).receivers(i); % receiver index
        lfp_R = sq(sess_data_lfp.lfp_E(r,:,:));
        
        % coherence computed using all the trials - errobar: Jackknife method
        [c_sr,f] = coherency(lfp_S,lfp_R,[N W],fs,fk,pad,0.05,1,1);
        SR_coh = [SR_coh; c_sr];
        session(s).labels(i)
    end 
      
    
end 


coher.mean_SR = mean(abs(SR_coh));
coher.mean_grand = mean(grand_coh);
coher.diff = coher.mean_SR - coher.mean_grand;

save(strcat(dir_RS_Theta,'/SR_coh_and_grand_coherence.mat'),'coher');

    
fig = figure;
plot(f,coher.mean_SR)
hold on
plot(f,coher.mean_grand)
hold on
plot(f,coher.diff)
hold on
xlim([0, 95])gend('SR','grand','diff')
title(sprintf('%s SR mean coh vs grand coh across all sess',monkey),'FontSize',10)
grid on
fname = strcat(dir_RS_Theta,sprintf('/SR_and_grand_coherence.fig',Sess));
saveas(fig,fname)
fname = strcat(dir_RS_Theta,sprintf('/SR_and_grand_coherence.jpg',Sess));
saveas(fig,fname)


keyboard
%%%%%%%%%%%%%%%%%%%%%%%%
% BOTH MONKEYS 
% %%%%%%%%%%%%%%%%%%%%%%

dir_RS_Theta = strcat(dir_main,sprintf('/Maverick/Resting_state/%s',freq_band));
mav = load(strcat(dir_RS_Theta,'/SR_coh_and_grand_coherence.mat'));
dir_RS_Theta = strcat(dir_main,sprintf('/Archie/Resting_state/%s',freq_band));
arc = load(strcat(dir_RS_Theta,'/SR_coh_and_grand_coherence.mat'));


mean_SR = mean([mav.coher.mean_SR; arc.coher.mean_SR]);
grand = mean([mav.coher.mean_grand; arc.coher.mean_grand]);
diff = mean_SR - grand;


fig = figure;
plot(f,mean_SR)
hold on
plot(f,grand)
hold on
plot(f,diff)
hold on
xlim([0, 95])
legend('SR','grand','diff')
title(sprintf('Both Monkeys: SR mean coh vs grand coh across all sess',monkey),'FontSize',10)
grid on

dir_both = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/both_monkeys/theta_band';
fname = strcat(dir_both,sprintf('/SR_and_grand_coherence.fig',Sess));
saveas(fig,fname)
fname = strcat(dir_both,sprintf('/SR_and_grand_coherence.jpg',Sess));
saveas(fig,fname)
