
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots the "Grand Coherence" averaged across all sessions and
% averaged across monkeys
%
%    @ Gino Del Ferraro, Dec 2022, Pesaran lab, NYU
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
dir_out = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/both_monkeys/grand_coherence/'

freq_band = 'theta_band';


%%%%%%%%%%%%%%%%
% ARCHIE 
% %%%%%%%%%%%%%%

monkey = 'Archie';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

coher = [];
for s = 1:length(sess_info{1})
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Grand_Coh = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/grand_coherence/',Sess));
    load(strcat(dir_Grand_Coh,sprintf('grand_coherence_sess_%d.mat',Sess))); % RS LFP split into 1 sec window and artifacts removed
    coher = [coher; coh]; 
end


%%%%%%%%%%%%%%%%
% MAVERICK
% %%%%%%%%%%%%%%

monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);

for s = 1:length(sess_info{1})
    
    Sess = sess_info{1}(s); % Session number
    display(['-- Session ',num2str(s),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Grand_Coh = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/grand_coherence/',Sess));
    load(strcat(dir_Grand_Coh,sprintf('grand_coherence_sess_%d.mat',Sess))); % RS LFP split into 1 sec window and artifacts removed
    coher = [coher; coh];  
end

mean_coh = mean(abs(coher));
err_coh= std(abs(coher))/sqrt(size(coher,1));
f = linspace(1,200,size(mean_coh,2)); % frequency values (range)

% -- FIGURE: Plot grand coherence across monkeys and sessions
set(0,'DefaultFigureVisible','on')
fig = figure;
shadedErrorBar(f,mean_coh,err_coh,'lineprops',{'color',[0, 51, 0]/255},'patchSaturation',0.4); hold on

grid on
title('Grand Coherence across monkeys and sessions','FontSize',11);
xlabel('freq (Hz)');
ylabel('coherence');
legend('grand coherence','FontSize',10)
set(gcf, 'Position',  [100, 600, 600, 400])
grid on
xlim([0, 95])
ylim([0, 0.18])

fname = strcat(dir_out,sprintf('/grand_coherence_across_all.pdf'));
saveas(fig,fname)
fname = strcat(dir_out,sprintf('/grand_coherence_across_all.fig'));
saveas(fig,fname)



