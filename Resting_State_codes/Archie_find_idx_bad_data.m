clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Gino_codes')
dir_main = '/mnt/pesaranlab/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';

freq_band = 'theta_band';
monkey = 'Archie';
dir_RS = strcat(dir_main,sprintf('%s/Resting_state/%s',monkey,freq_band));

fid = fopen(strcat(dir_RS,'/Sessions_with_modulator_info_movie.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


% -- define list of sessions
if strcmp(monkey,'Maverick')
    list_sess = 1:19;
    list_sess(17) = [];
else
    list_sess = 1:length(sess_info{3});
end


% -- exclude bad sessions
excluded_sess = [8,22,30,31];
excluded_idx = [2,5,8,9];
sess_list = 1:size(sess_info{1},1);
sess_list(excluded_idx) = [];

cnt = 0;
cnt_idx = [];
for i= list_sess
    
    Sess = sess_info{1}(i); % Session number
    %     display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Sess = strcat(dir_RS,sprintf('/Sess_%d/Modulators',Sess));
    
    load(strcat(dir_Sess,'/session_data_info.mat')); % --- dataG: all data info and LFP
    sess_data
    for Ch = sess_data.mod_idx
        
        if Ch ~= sess_data.receiver_idx
            cnt = cnt + 1;
            
            if any(excluded_sess == Sess)
                cnt_idx = [cnt_idx, cnt];
            end
            if Sess == 19 && Ch == 29
                cnt_idx = [cnt_idx, cnt];
            elseif Sess == 29 && Ch == 68
                cnt_idx = [cnt_idx, cnt];
            elseif Sess == 41 && Ch == 8
                cnt_idx = [cnt_idx, cnt];
            end
            
        end
    end
    
end

cnt

writematrix(cnt_idx,strcat(dir_RS,'/Modulators_Controls/Archie_modulator_idx_bad_data.txt'),'Delimiter',' ');






