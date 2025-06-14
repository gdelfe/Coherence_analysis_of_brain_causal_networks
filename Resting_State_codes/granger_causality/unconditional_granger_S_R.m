%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%  @ Gino Del Ferraro, July 2024, NYU, Pesaran Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% set(0,'DefaultFigureVisible','off')
set(0,'DefaultFigureVisible','on')
set(0,'DefaultLineLineWidth',2)

%%%%%%%%%%%%%%%%%%%
% - LOAD DATA --- %
%%%%%%%%%%%%%%%%%%%

addpath('T:/People/Gino/Coherence_modulator_analysis/Gino_codes');
% Main directory path
dir_main = '/vol/brains/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
% dir_main = 'T:\People\Gino\Coherence_modulator_analysis\Shaoyu_data';

% File names and parameters
name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % Filename for sess_data_info.mat
recording = 'last_recording';

% Frequency band and monkey information
freq_band = 'theta_band';
monkey = 'Maverick';
output_filename = "gc_SR_struct.mat"

% Construct the directory path for Resting State Theta
dir_RS_Theta = fullfile(dir_main, monkey, 'Resting_state', freq_band);

% Maximum lag for the computation of GC test
Lag = 50;

% Full path to the session info file
session_info_file = fullfile(dir_RS_Theta, 'Sessions_with_modulator_info_movie.txt');
session_info_file = strrep(session_info_file, '\', '/');  % Ensure all separators are consistent

% Attempt to open the session info file
fid = fopen(session_info_file);
if fid == -1
    error('Failed to open file. Check if the file exists and the path is correct: %s', session_info_file);
end

% Read the session information from the file
sess_info = textscan(fid, '%d%s%s'); % sess label, date, RS label

% Close the file after processing
fclose(fid);

% Display the contents of sess_info (optional)
disp(sess_info);


sess_cnt = 0; sr_cnt = 1; rs_cnt = 1; % Session count, SR granger count, RS granger count
gc = {};
for i = 1:size(sess_info{1},1)  % For each session with at least one modulator

    close all
    Sess = sess_info{1}(i); % Session number
    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),', out of tot  ',num2str(size(sess_info{1},1)),' sessions'])
    dir_Modulators = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators',Sess));
    load(strcat(dir_Modulators,name_struct_input)); % RS LFP split into 1 sec window and artifacts removed

    dir_Mod_recording = strcat(dir_RS_Theta,sprintf('/Sess_%d/Modulators/%s',Sess,recording));

    if ~exist(dir_Mod_recording, 'dir')
        mkdir(dir_Mod_recording)
    end

    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- true mod_Ch:  ',num2str(mod_Ch)])
    sess_data_lfp


    sess_cnt = sess_cnt + 1;

    close all
    lfp_S = sess_data_lfp.lfp_S; % Sender lfp
    lfp_R = sess_data_lfp.lfp_R; % Receiver lfp

    outliers_SR = [sess_data_lfp.outliers_S, sess_data_lfp.outliers_R]; % get outliers trials in Sender, Receiver, Modulator
    outliers_SR = unique(outliers_SR);  % -- remove repeated entries in outliers

    % -- remove outliers from sender, receiver, modulator
    lfp_S(outliers_SR,:) = [];
    lfp_R(outliers_SR,:) = [];


    Xtemp = cat(3, lfp_S, lfp_R);
    X = permute(Xtemp, [3, 2, 1]);  % vector for time series data

    %% Parameters

    ntrials   = size(X,3);     % number of trials
    nobs      = size(X,2);   % number of observations per trial

    regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
    icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

    morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
    momax     = 20;     % maximum model order for model order estimation

    acmaxlags = 2000;   % maximum autocovariance lags (empty for automatic calculation)

    tstat     = 'F';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
    alpha     = 0.05;   % significance level for significance test
    mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

    fs        = 1000;    % sample rate (Hz)
    fres      = 1000;     % frequency resolution (empty for automatic calculation)

    seed      = 0;      % random seed (0 for unseeded)


    %% Model order estimation (<mvgc_schema.html#3 |A2|>)

    % Calculate information criteria up to specified maximum model order.

    % ptic('\n*** tsdata_to_infocrit\n');
    [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
    % ptoc('*** tsdata_to_infocrit took ');

    dir_model = fullfile(dir_main, sprintf('%s/Resting_state/granger_gino/model_order/SR/', monkey));
    if ~exist(dir_model, 'dir')
        mkdir(dir_model);
    end

    % Plot information criteria.

    fig = figure(1); clf;
    plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
    title(sprintf('Model order estimation\nSess = %d',i));

    fprintf('\nbest model order (AIC) = %d\n',moAIC);
    fprintf('best model order (BIC) = %d\n',moBIC);
    saveas(fig,fullfile(dir_model,sprintf('Model_order_sess_%d.jpg',i)));

    % Select model order.

    if strcmpi(morder,'AIC')
        morder = moAIC;
        fprintf('\nusing AIC best model order = %d\n',morder);
    elseif strcmpi(morder,'BIC')
        morder = moBIC;
        fprintf('\nusing BIC best model order = %d\n',morder);
    else
        fprintf('\nusing specified model order = %d\n',morder);
    end


    %% VAR model estimation (<mvgc_schema.html#3 |A2|>)

    % Estimate VAR model of selected order from data.

    ptic('\n*** tsdata_to_var... ');
    [A,SIG] = tsdata_to_var(X,morder,regmode);
    ptoc;

    % Check for failed regression

    assert(~isbad(A),'VAR estimation failed');

    % NOTE: at this point we have a model and are finished with the data! - all
    % subsequent calculations work from the estimated VAR parameters A and SIG.

    %% Autocovariance calculation (<mvgc_schema.html#3 |A5|>)

    % The autocovariance sequence drives many Granger causality calculations (see
    % next section). Now we calculate the autocovariance sequence G according to the
    % VAR model, to as many lags as it takes to decay to below the numerical
    % tolerance level, or to acmaxlags lags if specified (i.e. non-empty).

    ptic('*** var_to_autocov... ');
    [G,info] = var_to_autocov(A,SIG,acmaxlags);
    ptoc;

    % The above routine does a LOT of error checking and issues useful diagnostics.
    % If there are problems with your data (e.g. non-stationarity, colinearity,
    % etc.) there's a good chance it'll show up at this point - and the diagnostics
    % may supply useful information as to what went wrong. It is thus essential to
    % report and check for errors here.

    var_acinfo(info,true); % report results (and bail out on error)

    %% Granger causality calculation: time domain  (<mvgc_schema.html#3 |A13|>)

    % Calculate time-domain pairwise-conditional causalities - this just requires
    % the autocovariance sequence.

    ptic('*** autocov_to_pwcgc... ');
    F = autocov_to_pwcgc(G);
    ptoc;

    % Check for failed GC calculation

    assert(~isbad(F,false),'GC calculation failed');

    % Significance test using theoretical null distribution, adjusting for multiple
    % hypotheses.
    nvars = size(X,1); % nvars = 3 in our case: Sender, Receiver, Modulator
    ntrials = size(X,3);

    pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
    sig  = significance(pval,alpha,mhtc);



    % Plot time-domain causal graph, p-values and significance.

    fig = figure(2); clf;
    set(gcf, 'Position', [100, 100, 1400, 400]); % [left, bottom, width, height]
    sgtitlex(sprintf('Pairwise Granger causality - time domain \nSession = %d',i));
    subplot(1,4,1);
    plot_pw_Gino(F);
    title('Pairwise-conditional GC');
    subplot(1,4,2);
    plot_pw_Gino(pval);
    title('p-values');
    subplot(1,4,3);
    plot_pw_Gino(sig);
    title(['Significant at p = ' num2str(alpha)])
    subplot(1,4,4)
    Granger_graph(sig,i)


    dir_granger = fullfile(dir_main, sprintf('%s/Resting_state/granger/granger_graphs/SR/', monkey));
    if ~exist(dir_granger, 'dir')
        mkdir(dir_granger);
    end

    saveas(fig,fullfile(dir_granger,sprintf('Sess_%d_granger_graph.jpg',i)))



    cd = mean(F(~isnan(F)));

    fprintf('\ncausal density = %f\n',cd);

    %% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

    % Calculate spectral pairwise-conditional causalities at given frequency
    % resolution - again, this only requires the autocovariance sequence.

    ptic('\n*** autocov_to_spwcgc... ');
    [GC, freq_axis] = autocov_to_spwcgc(G,fres, [], fs);
    ptoc;

    if sig(1,2)== 1 % Receveir to sender
        gc(rs_cnt).RS = GC;
        gc(rs_cnt).freq = freq_axis;
        gc(rs_cnt).pval = pval(1,2);
        gc(rs_cnt).sess = Sess;
        gc(rs_cnt).sess_i = i;
        rs_cnt = rs_cnt + 1;
    end
    if sig(2,1)== 1 % Sender to receiver
        gc(sr_cnt).SR = GC;
        gc(sr_cnt).freq = freq_axis;
        gc(sr_cnt).pval = pval(2,1);
        gc(sr_cnt).sess = Sess;
        gc(sr_cnt).sess_i = i;
        sr_cnt = sr_cnt + 1;
    end
    % Check for failed spectral GC calculation

    assert(~isbad(GC,false),'spectral GC calculation failed');

    % Plot spectral causal graph.

    fig = figure(3); clf;
    set(gcf, 'Position', [300, 300, 1000, 1000]); % [left, bottom, width, height]
    sgtitlex(sprintf('Pairwise Granger causality - frequency domain\nSess = %d',i));
    plot_spw_Gino(GC,fs,[0,50]);
    saveas(fig,fullfile(dir_granger,sprintf('Sess_%d_granger_frequency.jpg',i)))

    %% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

    % Check that spectral causalities average (integrate) to time-domain
    % causalities, as they should according to theory.

    fprintf('\nfrequency-domain GC integration check... ');
    Fint = smvgc_to_mvgc(GC); % integrate spectral MVGCs
    amax = maxabs(F+Fint)/2;
    if amax < 1e-5; amax = 1; end % in case all GCs very small
    mre = maxabs(F-Fint)/amax;
    if mre < 1e-5
        fprintf('OK (maximum relative error ~ %.0e)\n',mre);
    else
        fprintf(2,'WARNING: high maximum relative error ~ %.0e\n',mre);
    end

    %end
end

RS_rate = (rs_cnt-1)/sess_cnt;
SR_rate = (sr_cnt-1)/sess_cnt;

gc(1).RS_rate = RS_rate;
gc(1).SR_rate = SR_rate;
gc(1).rs_cnt_tot = rs_cnt-1;
gc(1).sr_cnt_tot = sr_cnt-1;
gc(1).sess_tot = sess_cnt;

dir_name = sprintf('%s/Resting_state/granger/granger_analysis/', monkey);
file_path = make_dir_get_file_path(dir_main, dir_name, output_filename)
save(file_path,'gc')