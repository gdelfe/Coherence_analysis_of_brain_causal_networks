%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function picks control electrodes (indexes) in the same 
% brain area as the sender (all the available electrodes)
%
% OUTPUT: indexes of the controls
%         Brain region of the controls
%
% @ Gino Del Ferraro, March 2022, Pesaran Lab, NYU
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [brain_idx] = choose_ALL_control_same_Region(RecordPairMRIlabels,MRIlabels,receiver_idx,s_area,mod_Ch)
    
  
    brain_idx = MRIlabels.(s_area).ElecIndx;  % -- get the indexes of the electrodes in the same brain region as the sender
    
    brain_idx = setdiff(brain_idx,mod_Ch); % -- remove all the modulator indexes from the list of brain indexes for that region
    brain_idx(brain_idx == receiver_idx) = []; % -- esclude the receiver from the controls, in case it belongs to this brain area
    
        
end 