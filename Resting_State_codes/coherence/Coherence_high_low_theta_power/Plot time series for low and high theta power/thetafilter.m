function theta = thetafilter(data, sampling)
% theta band filter (4-10 Hz)

% theta = thetafilter(data, sampling)
%  sampling defaults to 2e4 (in Hz)
%

if nargin < 2; sampling = 2e4; end

theta = mtfilter(data,[1,3],sampling,7);
%lfp = lfp(:,1:sampling./fs:end);