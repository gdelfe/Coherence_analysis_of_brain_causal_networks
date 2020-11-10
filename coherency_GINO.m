function [coh,f,S_X,S_Y,coh_err,SX_err,SY_err]=...
    coherency(X,Y,tapers,sampling,fk,pad,pval,flag,contflag)

% COHERENCY calculates the coherency between two time series, X and Y
%
% [COH, F, S_X, S_Y,, COH_ERR, SX_ERR, SY_ERR] = ...
%	COHERENCY(X, Y, TAPERS, SAMPLING, FK, PAD, PVAL, FLAG,CONTFLAG)
%
%  Inputs:  X		=  Time series array in [Space/Trials,Time] form.
%	    Y		=  Time series array in [Space/Trials,Time] form.
%	    TAPERS 	=  Data tapers in [K,TIME], [N,P,K] or [N, W] form.
%			   	Defaults to [N,5,9] where N is the duration
%				of X and Y.
%	    SAMPLING 	=  Sampling rate of time series X, in Hz.
%				Defaults to 1.
%	    FK 	 	=  Frequency range to return in Hz in
%                               either [F1,F2] or [F2] form.
%                               In [F2] form, F1 is set to 0.
%			   	Defaults to [0,SAMPLING/2]
%	    PAD		=  Padding factor for the FFT.
%			      	i.e. For N = 500, if PAD = 2, we pad the FFT
%			      	to 1024 points; if PAD = 4, we pad the FFT
%			      	to 2048 points.
%				Defaults to 2.
%	    PVAL	=  P-value to calculate error bars for.
%				Defaults to 0.05 i.e. 95% confidence.
%
%	    FLAG = 0:	calculate COH seperately for each channel/trial.
%	    FLAG = 1:	calculate COH by pooling across channels/trials.
%	    FLAG = 11 	calculation is done as for FLAG = 1
%		but the error bars cannot be calculated to save memory.
%	   	Defaults to FLAG = 11.
%     CONTFLAG = 1; There is only a single continuous signal coming in.
%               Defaults to 0.
%
%  Outputs: COH		=  Coherency between X and Y in [Space/Trials,Freq].
%                  F    =  Units of Frequency axis for COH
%	    S_X		=  Spectrum of X in [Space/Trials, Freq] form.
%	    S_Y		=  Spectrum of Y in [Space/Trials, Freq] form.
%	    COH_ERR 	=  Error bars for COH in [Hi/Lo, Space, Freq]
%			   form given by the Jacknife-t interval for PVAL.
% 	    SX_ERR 	=  Error bars for S_X.
% 	    SY_ERR 	=  Error bars for S_Y.
%

% Written by:  Bijan Pesaran Caltech 1998
%   Modified:  September 2003.
%

sX=size(X);
nt1=sX(2);
nch1=sX(1);

sY=size(Y);
nt2=sY(2);
nch2=sY(1);

if nt1 ~= nt2; error('Error: Time series are not the same length'); end
if nch1 ~= nch2; error('Error: Time series are incompatible'); end
nt = nt1;
nch = nch1;

if nargin < 4; sampling = 1; end
nt = nt./sampling;
if nargin < 3; tapers = [nt, 5, 9]; end
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
    tapers = dpsschk(tapers);
end
if nargin < 5 || isempty(fk); fk = [0,sampling./2]; end
if length(fk) == 1;  fk = [0,fk]; end
if nargin < 6 || isempty(pad); pad = 2; end
if nargin < 7 || isempty(pval); pval = 0.05;  end
if nargin < 8 || isempty(flag); flag = 11; end
if nargin < 9 || isempty(contflag); contflag = 0; end

N = length(tapers(:,1));
if N ~= nt*sampling
    error('Error: Tapers and time series are not the same length');
end

K = length(tapers(1,:));
nf = max(256,pad*2.^(nextpow2(N+1)));
nfk = floor(fk./sampling.*nf);
Nf = diff(nfk);
% Determine outputs
f = linspace(fk(1),fk(2),Nf);
errorchk = 0;
if nargout > 4; errorchk = 1; end

if flag == 0
    coh = zeros(nch, Nf);
    S_X = zeros(nch, Nf);
    S_Y = zeros(nch, Nf);
    
    if errorchk
        coh_err = zeros(2, nch, Nf);
        SX_err = zeros(2, nch, Nf);
        SY_err = zeros(2, nch, Nf);
    end
    
    if ~contflag
        mX = (sum(X,1)./nch)';
        mY = (sum(Y,1)./nch)';
    end
    for ch = 1:nch
        if contflag
            tmp1 = X(ch,:)'-sum(X(ch,:))./N;
            tmp2 = Y(ch,:)'-sum(Y(ch,:))./N;
        else
            tmp1 = X(ch,:)' - mX;
            tmp2 = Y(ch,:)' - mY;
        end
        Xk = fft(tapers(:,1:K).*tmp1(:,ones(1,K)),nf)';
        Xk = Xk(:,nfk(1)+1:nfk(2));
        Yk = fft(tapers(:,1:K).*tmp2(:,ones(1,K)),nf)';
        Yk = Yk(:,nfk(1)+1:nfk(2));
        SXk = Xk.*conj(Xk);
        S_X(ch,:) = sum(SXk,1)./K;
        SYk = Yk.*conj(Yk);
        S_Y(ch,:) = sum(SYk,1)./K;
        coh(ch,:) = sum(Xk.*conj(Yk),1)./K./sqrt(S_X(ch,:).*S_Y(ch,:));
        
        if errorchk	        %  Estimate error bars using Jacknife
            jcoh = zeros(K,Nf);
            jXlsp = zeros(K,Nf); jYlsp = zeros(K,Nf);
            for ik = 1:K
                indices = setdiff(1:K,ik);
                Xj = Xk(indices,:);
                Yj = Yk(indices,:);
                tmpx = sum(Xj.*conj(Xj),1)./(K-1);
                tmpy = sum(Yj.*conj(Yj),1)./(K-1);
                jcoh(ik,:)=atanh(abs(sum(Xj.*conj(Yj),1)./(K-1))./...
                    sqrt(tmpx.*tmpy));
                jXlsp(ik,:) = log(tmpx);
                jYlsp(ik,:) = log(tmpy);
            end
            lsigX = sqrt(K-1).*std(jXlsp,1);
            lsigY = sqrt(K-1).*std(jYlsp,1);
            lsigXY = sqrt(K-1).*std(jcoh,1);
            crit = tinv(1-pval./2,K-1);	%   Determine the scaling factor
            coh_err(1,ch,:) = tanh(atanh(abs(coh))+crit.*lsigXY);
            coh_err(2,ch,:) = tanh(atanh(abs(coh))-crit.*lsigXY);
            SX_err(1,ch,:) = exp(log(S_X)+crit.*lsigX);
            SX_err(2,ch,:) = exp(log(S_X)-crit.*lsigX);
            SY_err(1,ch,:) = exp(log(S_Y)+crit.*lsigY);
            SY_err(2,ch,:) = exp(log(S_Y)-crit.*lsigY);
        end
    end
end


if flag	== 1		% Pooling across trials
    Xk = zeros(nch*K, Nf);
    Yk = zeros(nch*K, Nf);
    if ~contflag
        mX = (sum(X,1)./nch)';
        mY = (sum(Y,1)./nch)';
    end
    for ch=1:nch
        if contflag
            tmp1 = X(ch,:)'-sum(X(ch,:))./N;
            tmp2 = Y(ch,:)'-sum(Y(ch,:))./N;
        else
            tmp1 = X(ch,:)' - mX;
            tmp2 = Y(ch,:)' - mY;
        end
        xk = fft(tapers(:,1:K).*tmp1(:,ones(1,K)),nf)';
        Xk((ch-1)*K+1:ch*K,:) = xk(:,nfk(1)+1:nfk(2));
        yk = fft(tapers(:,1:K).*tmp2(:,ones(1,K)),nf)';
        Yk((ch-1)*K+1:ch*K,:) = yk(:,nfk(1)+1:nfk(2));
    end
    S_X = sum(Xk.*conj(Xk),1)./K;
    S_Y = sum(Yk.*conj(Yk),1)./K;
    coh = sum(Xk.*conj(Yk),1)./K./sqrt(S_X.*S_Y);
    if errorchk			%  Estimate error bars using Jacknife
        jcoh = zeros(nch*K,Nf);  jXlsp = zeros(nch*K,Nf);
        jYlsp = zeros(nch*K,Nf);
        coh_err = zeros(2, Nf);
        SX_err = zeros(2, Nf);
        SY_err = zeros(2, Nf);
        for ik = 1:nch*K
            indices = setdiff(1:nch*K,ik);
            Xj = Xk(indices,:);  Yj = Yk(indices,:);
            tx = sum(Xj.*conj(Xj),1)./(nch*K-1);
            ty = sum(Yj.*conj(Yj),1)./(nch*K-1);
            %  Use atanh variance stabilizing transformation for coherence
            jcoh(ik,:)= atanh(abs(sum(Xj.*conj(Yj),1))./(nch*K-1)./...
                sqrt(tx.*ty));
            jXlsp(ik,:) = log(tx);
            jYlsp(ik,:) = log(ty);
        end
        lsigX = sqrt(nch*K-1).*std(jXlsp,1);
        lsigY = sqrt(nch*K-1).*std(jYlsp,1);
        lsigXY = sqrt(nch*K-1).*std(jcoh,1);
        crit = tinv(1-pval./2,nch*K-1);	%   Determine the scaling factor
        coh_err(1,:) = tanh(atanh(abs(coh))+crit.*lsigXY);
        coh_err(2,:) = tanh(atanh(abs(coh))-crit.*lsigXY);
        SX_err(1,:) = exp(log(S_X)+crit.*lsigX);
        SX_err(2,:) = exp(log(S_X)-crit.*lsigX);
        SY_err(1,:) = exp(log(S_Y)+crit.*lsigY);
        SY_err(2,:) = exp(log(S_Y)-crit.*lsigY);
    end
end


if flag == 11	%  Pooling across trials saving memory
    S_X = zeros(1,Nf);
    S_Y = zeros(1,Nf);
    coh = zeros(1,Nf);
    if ~contflag
        mX = (sum(X,1)./nch)';
        mY = (sum(Y,1)./nch)';
    end
    for ch = 1:nch
        if contflag
            tmp1 = X(ch,:)'-sum(X(ch,:))./N;
            tmp2 = Y(ch,:)'-sum(Y(ch,:))./N;
        else
            tmp1 = X(ch,:)' - mX;
            tmp2 = Y(ch,:)' - mY;
        end
        Xk = fft(tapers(:,1:K).*tmp1(:,ones(1,K)),nf)';
        S_X = S_X + sum(Xk(:,nfk(1)+1:nfk(2)).*conj(Xk(:,nfk(1)+1:nfk(2))),1)./K./nch;
        Yk = fft(tapers(:,1:K).*tmp2(:,ones(1,K)),nf)';
        S_Y = S_Y + sum(Yk(:,nfk(1)+1:nfk(2)).*conj(Yk(:,nfk(1)+1:nfk(2))),1)./K./nch;
        coh = coh + sum(Xk(:,nfk(1)+1:nfk(2)).*conj(Yk(:,nfk(1)+1:nfk(2))),1)./K./nch;
    end
    coh = coh./(sqrt(S_X.*S_Y));
end
end
