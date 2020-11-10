function [tfproj, f] = tfsp_proj(X, tapers, sampling, dn, fk, pad)
%TFSP_PROJ  Spectral projection using the multitaper techniques
%
% [TFPROJ, F] = TFSP_PROJ(X, TAPERS, SAMPLING, DN, FK, PAD)
%
%  X: Input data matrix: Trial, Time for one channel. 
%		If multiple channels, assumes 1 trial: Ch, Time.

%  TAPERS 	=  Data tapers in [K,TIME], [N,W] form.
%			   	    [N,W] Form:  N = duration of analysis window in s.
%                                W = bandwidth of frequency smoothing in Hz.
%               Defaults to [N,3,5] where N is NT/10
%				and NT is duration of X.
%
%  TFPROJ: 4D Matrix.  Trial,Time,K,Freq
%

% Modification History: Written by Bijan Pesaran 3 July, 2000

ntr = size(X, 1);
nt = size(X, 2);

if nargin < 3
    sampling = 1;
end
nt = nt./sampling;
if nargin < 2
    tapers = [nt,3,5];
end
if length(tapers) == 2
    n = tapers(1);
    w = tapers(2);
    p = n*w;
    k = floor(2*p-1);
    tapers = [n,p,k];
end
if length(tapers) == 3
    tapers(1) = tapers(1).*sampling;
    tapers = dpsschk(tapers);
end
if nargin < 4 || isempty(dn)
    dn = n./10;
end
if nargin < 5 || isempty(fk)
    fk = [0,sampling./2];
end
if length(fk) == 1
    fk = [0,fk];
end
if nargin < 6 || isempty(pad)
    pad = 2;
end

N = length(tapers(:,1));
nt = nt.*sampling;
K = length(tapers(1,:));
display(['N = ',num2str(N),' Number of tapers K = ',num2str(K)])
dn = floor(dn.*sampling);
nf = max(256,pad*2.^(nextpow2(N+1)));
nfk = floor(fk./sampling.*nf);
nfk1 = nfk(1); nfk2 = nfk(2);
nwin = floor((nt-N)./dn);           % calculate the number of windows
tapers = single(tapers);
f = linspace(fk(1),fk(2),diff(nfk));

% Create input matrix (samples x timesteps x trials x tapers)

% Don't eliminate this for-loop. The time to create index variables
% dominates the savings.
X = bsxfun(@minus, X, sum(X)./ntr).'; % demean input
starts = 1 + (0:(nwin-1))*dn;
stops = starts + N - 1;
Xin = zeros(N,nwin,1,ntr,'single');
for iWin = 1:nwin
    Xin(:,iWin,1,:) = X(starts(iWin):stops(iWin),:);
end
%(Time,Win,1,Trial)
Xin = bsxfun(@times, Xin, permute(tapers,[1 3 2]) );
%(Time,
% perform fft and select relevant frequencies and return in proper order
tfproj = fft(Xin,nf);
%Freq,Win,K,Trial
tfproj = tfproj(nfk1+1:nfk2,:,:,:);
tfproj = permute(tfproj, [4 2 3 1]);
%Trial,Win,K,Freq

