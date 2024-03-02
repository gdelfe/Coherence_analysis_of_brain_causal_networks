
% select trials with amplitude lower than a given threshold

function [amp,idx] = amplitude_L(lfp,th)

    amp = [];
    idx = [];

    for i=1:size(lfp,1)
        if max(abs(lfp(i,:)))< th
            amp = [amp; lfp(i,:)];
            idx = [idx; i];
        end
    end

end