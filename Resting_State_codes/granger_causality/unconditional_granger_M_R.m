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
monkey = 'Archie';

U = 'S';
V = 'R ctrl';

output_filename = "gc_RM_struct.mat"

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


sess_cnt = 0; 
gc = {};
MR_all = [];
RM_all = [];
cnt_m_tot = 0; % global modulator count across sessions

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


    sess_data_lfp
    
    cnt_m = 0; % local modulator count per session
    for ch = mod_Ch

        cnt_m = cnt_m + 1;

        if ~ismember(sess_data_lfp.receiver_idx,ch) % if receiver is not modulator

            cnt_m_tot = cnt_m_tot + 1;

            close all
            lfp_M = sq(sess_data_lfp.lfp_E(ch,:,:)); % Modulator lfp
            lfp_R = sess_data_lfp.lfp_R; % Receiver lfp

            outliers_MR = [sess_data_lfp.outliers_E(cnt_m).idx, sess_data_lfp.outliers_R]; % get outliers trials in Receiver, Modulator
            outliers_MR = unique(outliers_MR);  % -- remove repeated entries in outliers

            % -- remove outliers from receiver, modulator
            lfp_M(outliers_MR,:) = [];
            lfp_R(outliers_MR,:) = [];


            Xtemp = cat(3, lfp_M, lfp_R);
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
            plot_pw_Gino(F,U,V);
            title('Pairwise-conditional GC');
            subplot(1,4,2);
            plot_pw_Gino(pval,U,V);
            title('p-values');
            subplot(1,4,3);
            plot_pw_Gino(sig,U,V);
            title(['Significant at p = ' num2str(alpha)])
            subplot(1,4,4)
            Granger_graph(sig,i)


            dir_granger = fullfile(dir_main, sprintf('%s/Resting_state/granger/granger_graphs/MR/', monkey));
            if ~exist(dir_granger, 'dir')
                mkdir(dir_granger);
            end

            saveas(fig,fullfile(dir_granger,sprintf('Sess_%d_granger_graph_M_%d_R_%d.jpg',i,ch,sess_data_lfp.receiver_idx)))



            cd = mean(F(~isnan(F)));

            fprintf('\ncausal density = %f\n',cd);

            %% Granger causality calculation: frequency domain  (<mvgc_schema.html#3 |A14|>)

            % Calculate spectral pairwise-conditional causalities at given frequency
            % resolution - again, this only requires the autocovariance sequence.

            ptic('\n*** autocov_to_spwcgc... ');
            [GC, freq_axis] = autocov_to_spwcgc(G,fres, [], fs);
            ptoc;

            if sig(1,2)== 1 % Receveir to Modulator
                RM_all = [RM_all; GC(1,2,:)];
            end
            if sig(2,1)== 1 % Sender to receiver
                MR_all = [MR_all; GC(2,1,:)];
            end
            % Check for failed spectral GC calculation

            assert(~isbad(GC,false),'spectral GC calculation failed');

            % Plot spectral causal graph.

            fig = figure(3); clf;
            set(gcf, 'Position', [300, 300, 1000, 1000]); % [left, bottom, width, height]
            sgtitlex(sprintf('Pairwise Granger causality - frequency domain\nSess = %d',i));
            plot_spw_Gino(GC,fs,[0,50]);
            saveas(fig,fullfile(dir_granger,sprintf('Sess_%d_granger_frequency_M_%d_R_%d.jpg',i,ch,sess_data_lfp.receiver_idx)))

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

        end
    end
end


RM_all = sq(RM_all);
MR_all = sq(MR_all);

gc.RM_all = RM_all;
gc.MR_all = MR_all;

gc.RM_rate = size(RM_all,1)/cnt_m_tot;
gc.MR_rate = size(MR_all,1)/cnt_m_tot;
gc.cnt_m_tot = cnt_m_tot;
gc.freq = freq_axis;

dir_name = sprintf('%s/Resting_state/granger/granger_analysis/', monkey);
file_path = make_dir_get_file_path(dir_main, dir_name, output_filename)
save(file_path,'gc')

