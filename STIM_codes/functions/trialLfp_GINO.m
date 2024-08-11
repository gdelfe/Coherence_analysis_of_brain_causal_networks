function [Lfp] = trialLfp_cand(Trials,sys,electrode,contact,field,bn, MonkeyDir,lfpType, Num)
%  TRIALLFP loads lfp data for a trial
%
%  [LFP] = TRIALLFP(TRIALS, SYS, ELECTRODE, CONTACT, FIELD, BN, MONKEYDIR)
%
%  Inputs:	TRIALS = Trials data structure
%               SYS =   Scalar/String.  Recording system.
%                   Defaults to 1.  Could also be 'LIP'
%               ELECTRODE      = Scalar/Vector.  Electrodes(s) to load data from.
%                   Defaults to all electrodes in recording
%		CONTACT - Scalar/vector.  Contacts to load data from.
%		    Defaults to all contacts in a recording
%            	FIELD   = Scalar.  Event to align data to.
%                   Defaults to 'StartOn'
%            	BN      = Vector.  Time to start and stop loading data.
%                   Defaults to [-500,1500].
%               
%               MONKEYDIR = String.  Directory of project data
%                   Defaults to the global MONKEYDIR variable
%
%               lfpType = string of which lfp file to load (e.g. clfp,
%                         bplfp, etc. Defaults to [], in which case code searches for
%                         lfp files that exist and picks whichever is 'highest'
%                         priority.
%               Num     = scalar. number of segment within 'field' to use. default: 1
%               Num - scalar. Input to note segment within 'field' to load.
%               
%
%  Outputs:	LFP	= [TRIAL,NCH,TIME] or [TRIAL,TIME]. Lfp data
%
global MONKEYDIR;


if nargin < 2 || isempty(sys); sys = Trials(1).MT(1); end
if nargin < 3 || isempty(electrode); electrode = []; end
if nargin < 4 || isempty(contact); contact = 1; end
if nargin < 5 || isempty(field); field = 'StartOn'; end
if nargin < 6 || isempty(bn); bn = [-500,1500]; end
if ~exist('Num', 'var') || isempty(Num)
    Num = 1;
end

if nargin < 7 || isempty(MonkeyDir); MonkeyDir = MONKEYDIR; end
if nargin < 8
    lfpType = [];
end

if ~exist('Num', 'var') || isempty(Num)
    Num = 1;
end

lfpTypeOptions = {'dsplfp', 'clfp', 'bplfp', 'lfp', 'mflp'};
if ~isempty(lfpType)
    if ~ismember(lfpType, lfpTypeOptions)
        error(['Specified lfp type unsupported. Options are ' lfpTypeOptions{:}])
    end
end

if iscell(contact); contact = contact{1}; end

%Is this what works for Nan drive?
if ischar(sys) || iscell(sys)
    sysnum = 1;%findSys(Trials,sys);
else
    sysnum = sys;
end
if(iscell(sys))
    sys = sys{1};
end

ntr = length(Trials);

if isfield(Trials(1), 'HardwareType')
     switch Trials(1).HardwareType
        case 'neuropixel'
            channel_ids = Trials(1).ChannelID(sysnum).Electrode;
            CH = sum(~isnan(channel_ids));
        otherwise
            CH = length(getChannelIndex(Trials(1), sys)); % if you do not input electrode, returns total number of channels in system
            CH;            
     end
else    
    CH = length(getChannelIndex(Trials(1), sys)); % if you do not input electrode, returns total number of channels in system
    CH;   
end


Days = {Trials.Day};
Recs = {Trials.Rec};
day = Trials(1).Day;
rec = Trials(1).Rec;
if exist([MonkeyDir '/' day '/' rec '/rec' rec '.experiment.mat'],'file')
    load([MonkeyDir '/' day '/' rec '/rec' rec '.experiment.mat'])
    if isfield(experiment.hardware.acquisition,'lfp')
        samplingrate = experiment.hardware.acquisition.lfp.samplingrate;
    else
        samplingrate=1e3;
    end
            
else
    samplingrate = 1e3;
end

DayList = unique(Days)
if length(DayList)>1
    warning('trialLfp not written to load trial lists that span multiple days.')
end
RecList = unique(Recs);
dN = round(diff(bn)*samplingrate/1e3);
Lfp = zeros(ntr,length(electrode),dN);
if strcmp(field, 'PulseStarts')
    Lfp = zeros(5e4,length(electrode),dN);
    i = 1;
end

day = Trials(1).Day;
for iRec = 1:length(RecList)
    clfp_flag = 0;
    dsplfp_flag = 0;
    mlfp_flag = 0;
    bplfp_flag=0;

    rec = RecList{iRec};
    EventsFile = [MonkeyDir '/' day '/' rec '/rec' rec '.Events.mat'];
    MocapEventsFile = [MonkeyDir '/' day '/' rec '/rec' rec '.MocapEvents.mat'];
    SpontEventsFile = [MonkeyDir '/' day '/' rec '/rec' rec '.SpontEvents.mat'];
    StimEventsFile  = [MonkeyDir '/' day '/' rec '/rec' rec '.StimEvents_trcompat.mat'];
    

    if exist(EventsFile,'file')
        load(EventsFile);
    elseif exist(SpontEventsFile, 'file')
        load(SpontEventsFile);
        Events = SpontEvents;
    elseif exist(MocapEventsFile,'file')
        load(MocapEventsFile);
        Events = MocapEvents;
    elseif exist(StimEventsFile, 'file')
        load(StimEventsFile);
        Events = StimEvents;
        disp('Loading stim events...')
    else
        EventsFile
        error('No events file')
    end
    
    %check if trials are stim--if so load stimEvents
    stimtrflag = isStimTrial(Trials);
    if any(stimtrflag)
        if exist(StimEventsFile, 'file')
            load(StimEventsFile);
            Events = StimEvents;
            disp('Loading stim events...')
        end
    end
    
    fileprefix = [MonkeyDir '/' day '/' rec '/rec' rec '.' sys];
    if isempty(lfpType) %if user didn't specify type to load, find it for them
        if exist([fileprefix '.dsplfp.dat'],'file');
            dsplfp_flag = 1; %disp('Loading spike denoised data');
        elseif exist([fileprefix '.bplfp.dat'],'file');
            bplfp_flag = 1; %disp('Loading bandpass data');
        elseif exist([fileprefix '.clfp.dat'],'file');
            clfp_flag = 1; %disp('Loading cleaned  data');
        elseif exist([fileprefix '.mlfp.dat'],'file');
            mlfp_flag = 1; %disp('Loading median filtered data');
        end
    else %otherwise, set flag for specified type to 1
        eval([lfpType '_flag = 1;'])
    end
    
    RecInd = find(strcmp(Recs,RecList{iRec}));
    
    for tr = RecInd
        subtrial = Trials(tr).Trial;
        %disp(['Trial ' num2str(tr) ' Subtrial ' num2str(subtrial)]);

        if length(electrode)
            mtch = getChannelIndex(Trials(tr),sys,electrode,contact);
        else
            mtch = getChannelIndex(Trials(tr),sys,electrode,contact);
        end
        %fileprefix

        if dsplfp_flag
%             Lfp_tmp = loaddsplfp(fileprefix,Events,subtrial,field,bn,CH,samplingrate);
            Lfp_tmp = loaddsplfp(fileprefix,Events,subtrial,field,bn,CH,samplingrate, Num);
        elseif clfp_flag
            %disp('Clfp_flag')
            Lfp_tmp = loadclfp(fileprefix,Events,subtrial,field,bn,CH,samplingrate, Num);
        elseif bplfp_flag
%             Lfp_tmp = loadclfp(fileprefix,Events,subtrial,field,bn,CH,samplingrate);
            %disp('Bplfp_flag')
            Lfp_tmp = loadbplfp(fileprefix,Events,subtrial,field,bn,CH,samplingrate, Num);
        elseif mlfp_flag
%             Lfp_tmp = loadmlfp(fileprefix,Events,subtrial,field,bn,CH,samplingrate);
            Lfp_tmp = loadmlfp(fileprefix,Events,subtrial,field,bn,CH, samplingrate, Num);
        else
%             Lfp_tmp = loadlfp(fileprefix,Events,subtrial,field,bn,CH, samplingrate);
            Lfp_tmp = loadlfp(fileprefix,Events,subtrial,field,bn,CH, samplingrate, Num);
        end
        if strcmp(field, 'PulseStarts')
            if(CH == 1)
                for iPulse = 1:size(Lfp_tmp,1)
                    Lfp(i,:) = Lfp_tmp(iPulse, :);
                    i = i + 1;
                end
            else
                for iPulse = 1:size(Lfp_tmp,1)
                    Lfp(i,:,:) = Lfp_tmp(iPulse,mtch,:);
                    i = i + 1;
                end
            end
        else
            if(CH == 1)
                Lfp(tr,:,:) = Lfp_tmp(:);
            else
                if size(Lfp_tmp(1,mtch,:),3) == dN
                    Lfp(tr,:,:) = Lfp_tmp(1,mtch,:);
                else
                    Lfp(tr,:,:) = nan(1,length(mtch),dN);
                end
            end            
        end
    end
end

if strcmp(field, 'PulseStarts')
    if CH==1
        Lfp = Lfp(find(sq(sum(sum(Lfp,2),3))~=0),:);
    else
        Lfp = Lfp(find(sq(sum(sum(Lfp,2),3))~=0),:,:);
    end
end

if ntr || length(mtch)==1
    Lfp = sq(Lfp);
end

