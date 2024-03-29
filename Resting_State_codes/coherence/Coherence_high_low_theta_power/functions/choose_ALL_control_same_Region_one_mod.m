%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function picks control electrodes (indexes) in the same 
% brain area as the modulator (all the available electrodes)
%
% OUTPUT: indexes of the controls
%         Brain region of the controls
%
% @ Gino Del Ferraro, March 2022, Pesaran Lab, NYU
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [brain_idx,BrainReg] = choose_ALL_control_same_Region(RecordPairMRIlabels,MRIlabels,receiver_idx,m,mod_Ch)
    
    mod_Ch(mod_Ch == receiver_idx) = []; % -- remove the receiver from the modulators in case it is one of them
    
    BrainReg = RecordPairMRIlabels(m,1);  % -- get all the brain regions for the modulator(s)
    BrainReg = BrainReg(~cellfun('isempty',BrainReg)); % -- remove empty cells (in Archie)
    
    brain_idx = MRIlabels.(BrainReg{1}).ElecIndx;  % -- get the indexes of the electrodes in the same brain region
    
    brain_idx = setdiff(brain_idx,mod_Ch); % -- remove all the modulator indexes from the list of brain indexes for that region
    brain_idx(brain_idx == receiver_idx) = []; % -- esclude the receiver from the control modulators, in case it belongs to this brain area
    
        
end 