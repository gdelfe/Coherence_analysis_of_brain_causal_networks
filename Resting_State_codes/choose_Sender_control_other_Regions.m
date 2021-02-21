

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function picks a control electrode (index) in the all the brain 
% regions which are neither the modulator(s)' nor the receiver's
%
% OUTPUT: indexes of the controls 
%         brain areas of the controls 
%
% @ Gino Del Ferraro, December 2020, Pesaran Lab, NYU
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mod_Ch_rand,ControlsReg] = choose_Sender_control_other_Regions(RecordPairMRIlabels,MRIlabels,send_area,receiver_idx,mod_Ch)
    
    mod_Ch_rand = []; % -- list to store the control electrodes
    ModBrainReg = {}; % -- list to store the control electrodes' brain areas 
    mod_Ch(mod_Ch == receiver_idx) = []; % -- remove the receiver from the modulators in case it is one of them
    
    ModBrainReg = RecordPairMRIlabels(mod_Ch,1);  % -- get all the brain regions for the modulator(s)
    ReceiverReg = RecordPairMRIlabels(receiver_idx,1); % -- receiver brain region 
    ModBrainReg = [ModBrainReg; ReceiverReg] % -- add receiver region to modulators', in order to exclude it for the controls' regions
    ModBrainReg = unique(ModBrainReg);  % -- remove duplicates
    
    AllBrainReg = unique(RecordPairMRIlabels(:,1)); % -- all recorded regions for this session 
    ControlsReg = setdiff(AllBrainReg,ModBrainReg); % -- all other regions which are not modulators'
    
    L = length(ControlsReg);
    
    for reg = 1:L % -- for all the modulator(s) brain regions
        
        brainRegCtrl = ControlsReg{reg}; % -- get the brain region 
        brain_idx = MRIlabels.(brainRegCtrl).ElecIndx;  % -- get the indexes of the electrodes in the same brain region
       
        
        brain_idx(brain_idx == receiver_idx) = []; % -- esclude the receiver from the control modulators, in case it belongs to this brain area
        mod_Ch_rand = [mod_Ch_rand, brain_idx]; % -- return all the control electrodes in the same area as the modulator(s)'
        
    end
        
        
end 