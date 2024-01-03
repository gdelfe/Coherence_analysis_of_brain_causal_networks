%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function picks a control electrode (index) all the brain 
% regions which are neither in the sender nor the receiver's region
%
% OUTPUT: indexes of the controls 
%         brain areas of the controls 
%
% @ Gino Del Ferraro, June 2022, Pesaran Lab, NYU
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ctrl_el,ControlsReg] = receiver_control_other_region(RecordPairMRIlabels,MRIlabels,receiver_idx,sender_area)
    
    ctrl_el = []; % -- list to store the control electrodes
    
    receiver_area = RecordPairMRIlabels(receiver_idx,1); % -- receiver brain region 
    Reg = [sender_area; receiver_area]; % -- add receiver region to modulator', in order to exclude it for the controls' regions
    Reg = Reg(~cellfun('isempty',Reg)); % -- remove empty cells (in Archie)
    Reg = unique(Reg);  % -- remove duplicates
    
    RecordPairMRIlabels = RecordPairMRIlabels(~cellfun('isempty',RecordPairMRIlabels)); % -- remove empty cells (in Archie)
    AllBrainReg = unique(RecordPairMRIlabels(:,1)); % -- all recorded regions for this session 
    ControlsReg = setdiff(AllBrainReg,Reg) % -- all other regions which are not modulators'
    
    L = length(ControlsReg);
    
    for reg = 1:L % -- for all the control(s) brain regions
        
        brainRegCtrl = ControlsReg{reg}; % -- get the brain region 
        brain_idx = MRIlabels.(brainRegCtrl).ElecIndx;  % -- get the indexes of the electrodes in the same brain region
       
%         brain_idx = setdiff(brain_idx,mod_Ch); % -- remove all the modulator indexes from the list of brain indexes for that region

%         brain_idx(brain_idx == receiver_idx) = [] % -- esclude the receiver from the control modulators, in case it belongs to this brain area
        ctrl_el = [ctrl_el, brain_idx]; % -- return all the control electrodes in the same area as the modulator(s)'
        
    end
        
        
end 