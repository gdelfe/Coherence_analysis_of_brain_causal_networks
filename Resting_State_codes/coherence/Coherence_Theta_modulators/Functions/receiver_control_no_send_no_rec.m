%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function picks a control electrode (index) in the all the brain 
% regions. It only exclude the sender and the receiver itself (actually the
% sender never appears anyway, it is not in the list of electrodes), so it
% only needs to exclude the receiver electrode.
%
% This operation can also be computed in a simpler way, without this
% function, just by excluding the receiver from all recorded electrodes. 
% OUTPUT: indexes of the controls 
%         brain areas of the controls 
%
% @ Gino Del Ferraro, June 2022, Pesaran Lab, NYU
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ctrl_el_NO_SR,ctrl_Reg_NO_SR] = receiver_control_no_send_no_rec(RecordPairMRIlabels,MRIlabels,receiver_idx,sender_area)
    
    ctrl_el = []; % -- list to store the control electrodes
    
    RecordPairMRIlabels = RecordPairMRIlabels(~cellfun('isempty',RecordPairMRIlabels)); % -- remove empty cells (in Archie)
    AllBrainReg = unique(RecordPairMRIlabels(:,1)); % -- all recorded regions for this session 
    
    L = length(AllBrainReg);
    
    for reg = 1:L % -- for all the control(s) brain regions
        
        brainRegCtrl = AllBrainReg{reg}; % -- get the brain region 
        brain_idx = MRIlabels.(brainRegCtrl).ElecIndx;  % -- get the indexes of the electrodes in the same brain region
       
%         brain_idx = setdiff(brain_idx,mod_Ch); % -- remove all the modulator indexes from the list of brain indexes for that region

        ctrl_el = [ctrl_el, brain_idx]; % -- return all the control electrodes in the same area as the modulator(s)'
        
    end
        
    ctrl_el(ctrl_el == receiver_idx) = [] % -- esclude the receiver from the controls
    ctrl_el_NO_SR = ctrl_el;
    ctrl_Reg_NO_SR = AllBrainReg;
end 


