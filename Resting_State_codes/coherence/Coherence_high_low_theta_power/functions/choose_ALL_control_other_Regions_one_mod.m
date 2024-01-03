%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function picks a control electrode (index) in all the brain 
% regions which are neither the modulator(s)' nor the receiver's
%
% OUTPUT: indexes of the controls 
%         brain areas of the controls 
%
% @ Gino Del Ferraro, December 2020, Pesaran Lab, NYU
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mod_Ch_rand,ControlsReg] = choose_ALL_control_other_Regions(RecordPairMRIlabels,MRIlabels,receiver_idx,m,mod_Ch)
    
    mod_Ch_rand = []; % -- list to store the control electrodes
    mod_Ch(mod_Ch == receiver_idx) = []; % -- remove the receiver from the modulators in case it is one of them
    
    ModBrainReg = RecordPairMRIlabels(mod_Ch,1);  % -- get the brain region for the modulator
    ReceiverReg = RecordPairMRIlabels(receiver_idx,1); % -- receiver brain region 
    ModBrainReg = [ModBrainReg; ReceiverReg]; % -- add receiver region to modulator', in order to exclude it for the controls' regions
    ModBrainReg = ModBrainReg(~cellfun('isempty',ModBrainReg)); % -- remove empty cells (in Archie)
    ModBrainReg = unique(ModBrainReg);  % -- remove duplicates
    
    RecordPairMRIlabels = RecordPairMRIlabels(~cellfun('isempty',RecordPairMRIlabels)); % -- remove empty cells (in Archie)
    AllBrainReg = unique(RecordPairMRIlabels(:,1)); % -- all recorded regions for this session 
    ControlsReg = setdiff(AllBrainReg,ModBrainReg) % -- all other regions which are not modulators'
    
    L = length(ControlsReg);
    
    for reg = 1:L % -- for all the control(s) brain regions
        
        brainRegCtrl = ControlsReg{reg} % -- get the brain region 
        brain_idx = MRIlabels.(brainRegCtrl).ElecIndx;  % -- get the indexes of the electrodes in the same brain region
       
        brain_idx = setdiff(brain_idx,mod_Ch); % -- remove all the modulator indexes from the list of brain indexes for that region

        brain_idx(brain_idx == receiver_idx) = [] % -- esclude the receiver from the control modulators, in case it belongs to this brain area
        mod_Ch_rand = [mod_Ch_rand, brain_idx]; % -- return all the control electrodes in the same area as the modulator(s)'
        
    end
        
        
end 