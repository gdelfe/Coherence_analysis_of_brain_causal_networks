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
dir_main = '/vol/bd5/People/Gino/Coherence_modulator_analysis/Shaoyu_data/';
% dir_main = 'T:\People\Gino\Coherence_modulator_analysis\Shaoyu_data';

% File names and parameters
name_struct_input = '/session_data_lfp.mat';
filename = '.mat'; % Filename for sess_data_info.mat
recording = 'last_recording';

% Frequency band and monkey information
freq_band = 'theta_band';
monkey = 'Maverick';

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
    %     % ---  time parameter
    tot_time = 150001;

    mod_Ch = sess_data_lfp.mod_idx; % -- modulators index

    % --- Remove channels with artifacts for Maverick
    if strcmp(monkey,"Maverick")
        if Sess == 19
            mod_Ch(mod_Ch == 60) = []; % remove channel
        end
        if Sess == 41
            mod_Ch(mod_Ch == 8) = []; % remove channel
        end
    end

    display(['-- Session ',num2str(i),' -- label: ',num2str(Sess),',  -- true mod_Ch:  ',num2str(mod_Ch)])

    lfp_E_all = sess_data_lfp.lfp_E; % All electrode lfps

    cnt_m = 1; % modulator count
    for Ch = mod_Ch % for all the modulators in the session

        if Ch ~= sess_data_lfp.receiver_idx % if modulator is not receiver

            close all
            lfp_S = sess_data_lfp.lfp_S; % Sender lfp
            lfp_R = sess_data_lfp.lfp_R; % Receiver lfp
            lfp_E = squeeze(lfp_E_all(Ch,:,:)); % Modulator lfp

            outliers_SRE = [sess_data_lfp.outliers_S, sess_data_lfp.outliers_R, sess_data_lfp.outliers_E(cnt_m).idx]; % get outliers trials in Sender, Receiver, Modulator
            outliers_SRE = unique(outliers_SRE);  % -- remove repeated entries in outliers

            % -- remove outliers from sender, receiver, modulator
            lfp_S(outliers_SRE,:) = [];
            lfp_R(outliers_SRE,:) = [];
            lfp_E(outliers_SRE,:) = [];

            Xtemp = cat(3, lfp_S, lfp_R, lfp_E);
            X = permute(Xtemp, [3, 2, 1]);  % vector for time series data

            %% Parameters

            ntrials   = size(X,3);     % number of trials
            nobs      = size(X,2);   % number of observations per trial

            regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
            icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

            morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
            momax     = 20;     % maximum model order for model order estimation

            acmaxlags = 2000;   % maximum autocovariance lags (empty for automatic calculation)

            tstat     = 'chi2';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
            alpha     = 0.05;   % significance level for significance test
            mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

            fs        = 1000;    % sample rate (Hz)
            fres      = [];     % frequency resolution (empty for automatic calculation)

            seed      = 0;      % random seed (0 for unseeded)


            %% Model order estimation (<mvgc_schema.html#3 |A2|>)

            % Calculate information criteria up to specified maximum model order.

            % ptic('\n*** tsdata_to_infocrit\n');
            [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
            % ptoc('*** tsdata_to_infocrit took ');

            dir_model = fullfile(dir_main, sprintf('%s\\Resting_state\\granger\\model_order\\', monkey));
            if ~exist(dir_model, 'dir')
                mkdir(dir_model);
            end

            % Plot information criteria.

            fig = figure(1); clf;
            plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
            title(sprintf('Model order estimation\nSess = %d, mod count = %d',i,cnt_m));

            fprintf('\nbest model order (AIC) = %d\n',moAIC);
            fprintf('best model order (BIC) = %d\n',moBIC);
            saveas(fig,fullfile(dir_model,sprintf('Model_order_sess_%d_cntm_%d.jpg',i,cnt_m)));

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
            sgtitlex(sprintf('Pairwise-conditional Granger causality - time domain \nSession = %d, mod cnt = %d',i,cnt_m));
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
            Granger_graph(sig,i,cnt_m)

            dir_granger = fullfile(dir_main, sprintf('%s\\Resting_state\\granger\\granger_graphs\\', monkey));
            if ~exist(dir_granger, 'dir')
                mkdir(dir_granger);
            end

            saveas(fig,fullfile(dir_granger,sprintf('Sess_%d_cntm_%d_granger_graph.jpg',i,cnt_m)))



            cd = mean(F(~isnan(F)));

            fprintf('\ncausal density = %f\n',cd);

            %% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

            % Calculate spectral pairwise-conditional causalities at given frequency
            % resolution - again, this only requires the autocovariance sequence.

            ptic('\n*** autocov_to_spwcgc... ');
            f = autocov_to_spwcgc(G,fres);
            ptoc;

            % Check for failed spectral GC calculation

            assert(~isbad(f,false),'spectral GC calculation failed');

            % Plot spectral causal graph.

            fig = figure(3); clf;
            set(gcf, 'Position', [300, 300, 1000, 1000]); % [left, bottom, width, height]
            sgtitlex(sprintf('Pairwise-conditional Granger causality - frequency domain\nSess = %d, mod count = %d',i,cnt_m));
            plot_spw_Gino(f,fs,[0,150]);
            saveas(fig,fullfile(dir_granger,sprintf('Sess_%d_cntm_%d_granger_frequency.jpg',i,cnt_m)))

            %% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)

            % Check that spectral causalities average (integrate) to time-domain
            % causalities, as they should according to theory.

            fprintf('\nfrequency-domain GC integration check... ');
            Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
            amax = maxabs(F+Fint)/2;
            if amax < 1e-5; amax = 1; end % in case all GCs very small
            mre = maxabs(F-Fint)/amax;
            if mre < 1e-5
                fprintf('OK (maximum relative error ~ %.0e)\n',mre);
            else
                fprintf(2,'WARNING: high maximum relative error ~ %.0e\n',mre);
            end


        end
    end
end