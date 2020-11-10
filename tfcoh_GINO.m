function [coh,tf,f,S_X,S_Y] = ...
    tfcoh_GINO(X,Y,tapers,sampling,dn,fk,pad,pval,flag,contflag,detrendflag)
%  TFCOH Moving window time-frequency coherency between two time series.
%
% [COH, F, S_X, S_Y] = ...
%	TFCOH(X, Y, TAPERS, SAMPLING, DN, FK, PAD, PVAL, FLAG, CONTFLAG,DETRENDFLAG)
%
%  Inputs:  X		=  Time series array in [Space/Trials,Time] form.
%	    Y		=  Time series array in [Space/Trials,Time] form.
%	    TAPERS 	=  Data tapers in [K,TIME], [N,P,K] or [N, W] form.
%			   	Defaults to [N,5,9] where N is the duration
%				of X and Y.
%	    SAMPLING 	=  Sampling rate of time series X, in Hz.
%				Defaults to 1.
%	    DN		=  Overlap in time between neighbouring windows.
%			       	Defaults to N./10;
%	    FK 	 	=  Frequency range to return in Hz in
%                               either [F1,F2] or [F2] form.
%                               In [F2] form, F1 is set to 0.
%			   	Defaults to [0,SAMPLING/2]
%	    PAD		=  Padding factor for the FFT.
%			      	i.e. For N = 500, if PAD = 2, we pad the FFT
%			      	to 1024 points; if PAD = 4, we pad the FFT
%			      	to 2048 points.
%				Defaults to 2.
%
%	   FLAG = 0:	calculate COH seperately for each channel/trial.
%	   FLAG = 1:	calculate COH by pooling across channels/trials.
%      FLAG = 11:   calculate COH by pooling across channels/trials without
%      error bars
%                       Defaults to 11, for which the code is optimized.
%      CONTFLAG = 1; There is only a single continuous signal coming in.
%               Defaults to 0.
%
%      DETRENDFLAG = 1; Detrend each trial with a straight line for each
%       window.
%               Defaults to 0
%
%
%  Outputs: COH	        =  Coherency between X and Y in [Space/Trials,Freq].
%	    F		=  Units of Frequency axis for COH.
%	    S_X		=  Spectrum of X in [Space/Trials, Freq] form.
%	    S_Y		=  Spectrum of Y in [Space/Trials, Freq] form.
%

%  Written by:  Bijan Pesaran Caltech 1998
%                   Optimized when not computing error bars and pooling
%                   across trials
%


sX = size(X);
nt1 = sX(2);
nch1 = sX(1);

sY = size(Y);
nt2 = sY(2);
nch2 = sY(1);
 
if nt1 ~= nt2 error('Error: Time series are not the same length'); end
if nch1 ~= nch2 error('Error: Time series are incompatible'); end
nt = nt1;
nch = nch1;

if nargin < 4 sampling = 1; end
t = nt./sampling;
if nargin < 3 tapers = [t,5,9]; end
if length(tapers) == 2
    n = tapers(1);
    w = tapers(2);
    p = n*w;
    k = floor(2*p-1);
    tapers = [n,p,k];
    %disp(['Using ' num2str(k) ' tapers.']);
end
if length(tapers) == 3
    tapers(1) = floor(tapers(1).*sampling);
    tapers = single(dpsschk(tapers));
end
if nargin < 5 || isempty(dn); dn = n./10; end
if nargin < 6 || isempty(fk); fk = [0,sampling./2]; end
if length(fk) == 1
    fk = [0,fk];
end
if nargin < 7 || isempty(pad); pad = 2; end
if nargin < 8 || isempty(pval); pval = 0.05;  end
if nargin < 9 || isempty(flag); flag = 11; end
if nargin < 10 || isempty(contflag)
    contflag = 0;
end
if nargin < 11 || isempty(detrendflag)
    detrendflag = 0;
end
K = length(tapers(1,:));
N = length(tapers(:,1));
 

if N > nt error('Error: Tapers are longer than time series'); end

% Determine outputs
errorchk = 0;
if nargout > 3 errorchk = 1; end

dn = floor(dn.*sampling);
nf = max(256, pad*2^nextpow2(N+1));
nfk = floor(fk./sampling.*nf);
nwin = floor((nt-N)./dn);           % calculate the number of windows
f = linspace(fk(1),fk(2),diff(nfk));
K = single(K);
nch = single(nch);


if flag == 0
    S_X = zeros(nwin,nch,diff(nfk),'single');
    S_Y = zeros(nwin,nch,diff(nfk),'single');
    coh = zeros(nwin,nch,diff(nfk),'single') + i*zeros(nwin,nch,diff(nfk),'single');
    mX = sum(X)./nch; mY = sum(Y)./nch;
    X = (X - mX(ones(1,nch),:)).';
    Y = (Y - mY(ones(1,nch),:)).';
    
    for iWin=1:nwin
        tmp1 = X(dn*(iWin-1)+1:dn*(iWin-1)+N,:);
        tmp2 = Y(dn*(iWin-1)+1:dn*(iWin-1)+N,:);
        for ch = 1:nch
            if detrendflag
                tmp1(:,ch) = detrend(tmp1(:,ch),'linear');
                tmp2(:,ch) = detrend(tmp2(:,ch),'linear');
            end
            Xk = fft(tapers.*tmp1(:,ch*ones(1,K)),nf);
            Xk = Xk(nfk(1)+1:nfk(2),:);
            S_X(iWin,ch,:) = sum(Xk.*conj(Xk),2)./(K.*nch);
            Yk = fft(tapers.*tmp2(:,ch*ones(1,K)),nf);
            Yk = Yk(nfk(1)+1:nfk(2),:);
            S_Y(iWin,ch,:) = sum(Yk.*conj(Yk),2)./(K.*nch);
            coh(iWin,ch,:) = sum(Xk.*conj(Yk),2)./(K.*nch)./sqrt(sq(S_X(iWin,ch,:).*S_Y(iWin,ch,:)));
        end
    end
elseif flag == 11
    %  This is optimized
    %disp(['Flag is ' num2str(flag)]);
    
    S_X = zeros(nwin,diff(nfk),'single');
    S_Y = zeros(nwin,diff(nfk),'single');
    coh = zeros(nwin,diff(nfk),'single') + i*zeros(nwin,diff(nfk),'single');
    %     mX = sum(X)./nch; mY = sum(Y)./nch;
    %     X = (X - mX(ones(1,nch),:)).';
    %     Y = (Y - mY(ones(1,nch),:)).';
    for win = 1:nwin
        % Here the optimized spectral loop starts.
        if contflag
            tmpX = detrend(X(:,dn*win:dn*win+N-1))';
            tmpY = detrend(Y(:,dn*win:dn*win+N-1))';
        else
            mX = sum(X(:,dn*(win-1)+1:dn*(win-1)+N),1)./nch;
            tmpX = (X(:,dn*(win-1)+1:dn*(win-1)+N) - mX(ones(1,nch),:)).';
            mY = sum(Y(:,dn*(win-1)+1:dn*(win-1)+N),1)./nch;
            tmpY = (Y(:,dn*(win-1)+1:dn*(win-1)+N) - mY(ones(1,nch),:)).';
        end
        SX = zeros(diff(nfk),1,'single');
        SY = zeros(diff(nfk),1,'single');
        c = zeros(diff(nfk),1,'single');
        for ch = 1:nch
            if detrendflag
                tmpX(:,ch) = detrend(tmpX(:,ch),'linear');
                tmpY(:,ch) = detrend(tmpY(:,ch),'linear');
            end
            Xk = fft(tapers.*tmpX(:,ch*ones(1,K)),nf);
            Xk = Xk(nfk(1)+1:nfk(2),:);
            A = sum(Xk.*conj(Xk),2)./(K.*nch);
            SX = SX + A;
            Yk = fft(tapers.*tmpY(:,ch*ones(1,K)),nf);
            Yk = Yk(nfk(1)+1:nfk(2),:);
            A = sum(Yk.*conj(Yk),2)./(K.*nch);
            SY = SY + A;
            B = sum(Xk.*conj(Yk),2)./(K.*nch);
            c = c + B;
        end
        coh(win,:) = (c./(sqrt(SX.*SY))).';
        S_X(win,:) = SX.'; S_Y(win,:) = SY.';
    end
end

tf = linspace(N/2,nt-N/2,nwin);

if flag == 0
    coh = permute(coh,[2,1,3]);
end
