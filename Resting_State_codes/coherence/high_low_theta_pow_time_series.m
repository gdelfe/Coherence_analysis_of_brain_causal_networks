
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code plots the time series (LFP) for low and high theta power of the
% modulator. High theta power time series should show oscillatory rhythms
% at the theta frequency. 
%
%    @ Gino Del Ferraro, November 2020, Pesaran lab, NYU
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
recording = 'last_recording';

freq_band = 'theta_band';
monkey = 'Maverick';
dir_RS_Theta = strcat(dir_main,sprintf('/%s/Resting_state/%s',monkey,freq_band));


fid = fopen(strcat(dir_RS_Theta,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


s = 16 % 1:size(sess_info{1},1)
Sess = sess_info{1}(s); % Session number
cnt_m = 1;

Sess = sess_info{1}(s); % Session number
dir_Sess_mod_send_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_send/Data',Sess));
dir_Sess_mod_rec_data = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_rec/Data',Sess));

dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));

load(strcat(dir_Modulators,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed


% for all the modulators in that session 
for cnt_m = 1:length(sess_data_lfp.mod_idx)
    
    load(strcat(dir_Sess_mod_send_data,sprintf('/modulator_%d_send.mat',cnt_m)));
    load(strcat(dir_Sess_mod_rec_data,sprintf('/modulator_%d_rec.mat',cnt_m)));
    
    % plot the time series lfp
%     fig_ts = figure;
%     for i = 0:4
%         
%         plot(mod_rec.lfp_E_high(end-i,:)+400-40*i,'color','b')
%         hold on
%         plot(mod_rec.lfp_E_low(i+1,:)+40*i,'color','r')
%         hold on
%         
%     end
    
    
%     
%     figure; 
%     title('Lfp high and low theta power trials','FontSize',12)
%     for i = 0:4
%         
%         ha(i+1) = subplot(10,1,i+1);
%         plot(mod_rec.lfp_E_high(end-i,:),'color','b')
%         set(gca,'xticklabel',[]);
%         grid on
%   
%         subplot(10,1,10-i);
%         plot(mod_rec.lfp_E_low(i+1,:),'color','r')
%         grid on 
%     end
%     
%     xlabel('Time (ms)','FontName','Arial','FontSize',12);
%     ylabel('Lfp','FontName','Arial','FontSize',12);
%     set(gcf, 'Position',  [100, 600, 1200, 800]);
    
    
    
    fig_ts = figure;
    [ha, pos] = tight_subplot(10,1,[.01 .03],[.2 .01],[.1 .1])
    for i = 0:4
        
        axes(ha(i+1))
        plot(mod_rec.lfp_E_high(end-i,:),'color','b')
        %set(gca,'xticklabel',[]);
        grid on
        
        axes(ha(10-i))
        plot(mod_rec.lfp_E_low(i+1,:),'color','r')
        grid on

    end
    % set(ha(1:5),'XTickLabel',''); set(ha,'YTickLabel','')
    
    ylabel('Lfp','FontName','Arial','FontSize',12);
    set(gcf, 'Position',  [100, 600, 1500, 1000]);
    
    


    
    % save the plots
    dir_ts_rec_fig = strcat(dir_high_low_theta,sprintf('/Sess_%d/mod_rec/Figures/Time_Series',Sess));
    if ~exist(dir_ts_rec_fig, 'dir')
        mkdir(dir_ts_rec_fig)
    end
    
    fname = strcat(dir_ts_rec_fig,sprintf('/modulator_%d_lfp_high_low_theta.jpg',cnt_m));
    saveas(fig_ts,fname);
    

    
end




