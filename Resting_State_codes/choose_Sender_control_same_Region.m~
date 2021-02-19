%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function picks control electrodes for the receiver in the same
% brain area as the receiver itself.
%
% OUTPUT: indexes of the controls
%         Brain region of the controls
%
% @ Gino Del Ferraro, February 2021, Pesaran Lab, NYU
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [send_Ch,RecBrainReg] = choose_Sender_control_same_Region(RecordPairMRIlabels,MRIlabels,send_area)
        
    brain_idx = MRIlabels.(send_area).ElecIndx;  % -- get the indexes of the electrodes in the same brain region
    brain_idx(brain_idx == receiver_idx) = []; % --- remove the index of the receiver itself 
    
    send_Ch = brain_idx; % -- return all the control electrodes in the same area as the receiver
    
end 