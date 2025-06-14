%% plot_spw
%
% Plot spectral pairwise quantities on a grid
%
% <matlab:open('plot_spw.m') code>
%
%% Syntax
%
%     plot_spw(P,fs)
%
%% Arguments
%
% _input_
%
%     P          matrix of spectral pairwise quantities
%     fs         sample rate in Hz (default: normalised freq as per routine 'sfreqs')
%     frange     frequency range to plot: empty for all (default) else an ascending 2-vector
%
%% Description
%
% Plot pairwise spectral quantities in |P|, a 3-dim numerical matrix with
% first index representing target ("to"), second index source ("from")
% quantities and third index frequencies - typically spectral causalities (see
% e.g. <autocov_to_spwcgc.html |autocov_to_spwcgc|>).
%
%% See also
%
% <autocov_to_spwcgc.html |autocov_to_spwcgc|> |
% <mvgc_demo.html |mvgc_demo|> |
% <sfreqs.html |sfreqs|>
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%%

function plot_spw_Gino(P,fs,frange,J,K)

n = size(P,1);
assert(ndims(P) == 3 && size(P,2) == n,'must be a 3-dim matrix with the first two dims square');
h = size(P,3);

if nargin < 2, fs = []; end % default to normalised frequency as per 'sfreqs'

if nargin < 3, frange = []; end; % default to all
if ~isempty(frange)
    assert(isvector(frange) && length(frange) == 2 && frange(1) < frange(2),'frequency range must be an ascending 2-vector of frequencies');
end

fres = h-1;
lam = sfreqs(fres,fs)';
if ~isempty(frange)
    idx = lam >= frange(1) & lam <= frange(2);
    lam = lam(idx);
    P = P(:,:,idx,:);
end
xlims = [lam(1) lam(end)];
ylims = [min(P(:)) 1.1*max(P(:))];
if isempty(fs), xlab = 'normalised frequency'; else xlab = 'frequency (Hz)'; end

V = {J,K};
k = 0;
for i = 1:n
    for j = 1:n
        k = k+1;
        if i ~= j
            subplot(n,n,k);
            plot(lam,squeeze(P(i,j,:)));
            axis('square');
            xlim(xlims);
            ylim(ylims);
            xlabel(xlab);
            title(sprintf('\n%s -> %s',V{j},V{i}));
        end
    end
end
