%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function picks control electrodes for the receiver in the same
% brain area as the receiver itself.
%
% OUTPUT: indexes of the controls
%         Brain region of the controls
%
% @ Gino Del Ferraro, February 2021, Pesaran Lab, NYU
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rec_Ch,RecBrainReg] = choose_Receiver_control_same_Region(RecordPairMRIlabels,MRIlabels,receiver_idx,mod_Ch)
    
    RecBrainReg = RecordPairMRIlabels(receiver_idx,1);  % -- get the receiver's brain region 
    brain_idx = MRIlabels.(RecBrainReg{1}).ElecIndx;  % -- get the indexes of the electrodes in the same brain region
    brain_idx(brain_idx == receiver_idx) = []; % --- remove the index of the receiver itself 
    brain_idx = setdiff(brain_idx,mod_Ch); % --- remove the index of the modulators, in case there is any in the same brain area
    
    
    rec_Ch = brain_idx; % -- return all the control electrodes in the same area as the receiver
    
end 