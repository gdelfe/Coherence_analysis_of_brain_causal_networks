

clear all; close all;

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%
addpath('/mnt/pesaranlab/People/Gino/DL-modulators/Gino_codes')
dir_base = '/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State';

fid = fopen(strcat(dir_base,'/Sessions_with_modulator_info.txt')); % load session info with no repetition
sess_info = textscan(fid,'%d%s%s'); % sess label, date, RS label
fclose(fid);


% --- load structure with data about session, modulators, send-rec
temp_AM = load(strcat(dir_base,'/session_AM.mat'),'session_AM');
temp_MA = load(strcat(dir_base,'/session_MA.mat'),'session_MA');

session_AM = temp_AM.session_AM;
session_MA = temp_MA.session_MA;



for i=1 % :size(sess_info{1},1) % for all the sessions with modulator
    
    addpath(sprintf('/vol/sas8/Maverick_RecStim_vSUBNETS220/%s/%s/',sess_info{2}{i},sess_info{3}{i})) % add path of the specific RS session
    
    % file = 'rec004.Frontal.lfp.dat'
    file = sprintf('rec%s.Frontal.lfp.dat',sess_info{3}{i})
    fid = fopen(file);
    format = 'float=>single';
    
    CH = 220; % tot number of channels
    FS = 1000; % sampling
    
    
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',   out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    
    data = fread(fid,[CH,inf],format); % load the RS data
    % h = fread(fid,[CH,diff(bn)*FS./1e3],format);
    % ---- bipolar referencing, pairs of electrodes
    elect_dir = sprintf('/mnt/pesaranlab/People/Gino/DL-modulators/Shaoyu_data/Resting_State/Sess_%d',Sess);
    
    
    x = TS(:,1); % RMED
    y = TS(:,2); % RMEV
    
    alpha = 0.000000000001;
    lapsemax = 30;
    
    Two_causes_One = []; % RIML causes RMEV
    One_causes_Two = []; % RMEV cause RIML
    
    for lapse = 1:lapsemax
        [F1,c_v1] = granger_cause(x,y,alpha,lapse);
        Two_causes_One = [Two_causes_One; F1 c_v1];
        
        [F2,c_v2] = granger_cause(y,x,alpha,lapse);
        One_causes_Two = [One_causes_Two; F2 c_v2];
    end
    
    diff_Two_causes_One = zeros(length(Two_causes_One),1);
    diff_One_causes_Two = zeros(length(One_causes_Two),1);
    k = zeros(length(One_causes_Two),1);
    
    for i= 1:lapsemax
        diff_Two_causes_One(i) =  Two_causes_One(i,1) - Two_causes_One(i,2);
        diff_One_causes_Two(i) = One_causes_Two(i,1)- One_causes_Two(i,2);
        k(i) = i;
    end
    
    figure(1)
    plot(k,diff_Two_causes_One,k,diff_One_causes_Two);
    grid on
    xticks([1:2:30]);
    fig = title(name_title, 'FontSize', 18); % set title
    xlabel('time lapse', 'FontSize', 15);
    ylabel('difference', 'FontSize', 15);
    h = legend('2 causes 1','1 causes 2');
    set(h,'FontSize', 15);
    
    saveas(fig,name_fig);
    
end