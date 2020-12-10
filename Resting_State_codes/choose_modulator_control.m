%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function picks a control electrode (index) in the same 
% brain area as the modulator.
% OUTPUT: indexes of the control modulators - as many as the true
% modulators. Each one of them in the same brain area
%
% @ Gino Del Ferraro, December 2020, Pesaran Lab, NYU

function [mod_Ch_rand] = choose_modulator_control(RecordPairMRIlabels,MRIlabels,receiver_idx,Ch,mod_Ch)
    
    mod_Ch_rand = [];
    mod_Ch_check = mod_Ch;
    
    for mod = mod_Ch
        
        brainRegions = RecordPairMRIlabels(:,1); % -- get all the brain regions
        brainRegMod = brainRegions{mod};  % -- get the brain region of the modulator in channel Ch
        brain_idx = MRIlabels.(brainRegMod).ElecIndx;  % -- get the indexes of the electrodes in the same brain region
        
        % -- remove all the modulator indexes from the list of brain indexes
        % This is done to avoid resampling another modulator
        for elec = mod_Ch_check
            brain_idx(brain_idx == elec) = [];
        end
        brain_idx(brain_idx == receiver_idx) = []; % -- esclude the receiver from the control modulators, in case it belongs to this brain area
        
        
        L = length(brain_idx);
        brain_idx_rand = brain_idx(randperm(L)); % -- randomly permute the indexes of the electrodes in the same brain region

        mod_Ch_rand = [mod_Ch_rand, brain_idx_rand(1)]; % return as many modulator controls as the real modulators number
        mod_Ch_check = [mod_Ch_check, mod_Ch_rand]; % -- add the rand control to list of modulators in order not to pick it twice (This is a very important line)
        
    end
    
end 