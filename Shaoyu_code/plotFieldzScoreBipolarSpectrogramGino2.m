function zSpec =  plotFieldzScoreBipolarSpectrogramGino(data,hitIndx,missIndx,bn,tapers,fs,dn,fk,pad,nPerm,alpha)

Spec = []; TestSpec = []; DiffSpec = []; AvDiffSpec = [];
keyboard
[Spec,f] = tfspec(data(hitIndx,:),tapers,fs,dn,fk,pad,0.05,0,1); % get the time-frequency spectrum of the hit trial
[spectrum, f2] = tfspec(Lfp,[0.5,2],1e3,0.05,200,2,0.05,1); % one single call


nWin = size(Spec,2);
TestSpec = sq(sum(log(Spec))./numel(hitIndx));

DiffSpec = tfspec(data(missIndx,:),tapers,fs,dn,fk,pad); % get the time-frequency spectrum of the miss trial 
AvDiffSpec =  sq(sum(log(DiffSpec)))./numel(missIndx);

figure; tvimage(log(Spec(10,:,:))); colorbar;
figure; tvimage(log(DiffSpec(10,:,:))); colorbar;

imagesc(log(Spec(1,:,:))); colorbar;

keyboard 

display('printing figure')
imagesc(flipud(zSpec'),clim)  
%colormap(gca,polarmap(winter,1))
colormap(gca,polarmap(jet,1.5))

figure; tvimage(log(spectrum)); colorbar;

display('permutation test')

keyboard
%% permutation test
D = []; Dperm = [];
D = TestSpec - AvDiffSpec;
numtrials = numel(hitIndx);
Diff_NumTr = numel(missIndx);
Tot_NumTr = numtrials + Diff_NumTr;
pSpec = []; zSpec = []; shuffleDperm = [];
for iWin = 1:nWin
    TotSpec = [sq(Spec(:,iWin,:));sq(DiffSpec(:,iWin,:))];
    Dperm = zeros(nPerm,size(DiffSpec,3));
    for iPerm = 1:nPerm
        Perm_ind = randperm(Tot_NumTr);
        Spec1 = sum(log(TotSpec(Perm_ind(1:numtrials),:)))./numtrials;
        Spec2 = sum(log(TotSpec(Perm_ind(numtrials+1:end),:)))./Diff_NumTr;
        Dperm(iPerm,:) = Spec1-Spec2;
    end
    
    for iF = 1:size(D,2)
        pSpec(iWin,iF) = length(find(abs(Dperm(:,iF))>abs(D(iWin,iF))))./nPerm;
    end
    pSpec(pSpec==0) = 1./nPerm;  pSpec(pSpec==1) = (nPerm-1)./nPerm;
    zSpec(iWin,:) = D(iWin,:)./std(Dperm,1);
    shuffleDperm(iWin,:,:) = Dperm;
end

display('end permutation test')

%% plot
%clim = [-ceil(max(max(abs(zSpec)))) ceil(max(max(abs(zSpec))))];
keyboard
clim = [-4 4];
%
NX = size(zSpec,1);
NY = size(zSpec,2);
%X = linspace(bn(1),bn(2),NX);
%Y = linspace(0,fk,NY);
%imagesc(X,Y,flipud(zSpec'),clim)

% Very important: be careful when you use zSpec, here it is flipup and transposed
figure;
display('printing figure')
imagesc(flipud(zSpec'),clim)  
%colormap(gca,polarmap(winter,1))
colormap(gca,polarmap(jet,1.5))
%colormap(gca,polarmap(polarmap,0.5))

%cmap = cbrewer('div','PuOr',9);
%colormap(cmap)

c = colorbar;
c.Label.String = 'Hit - Miss (Z-score)';

keyboard

FigDir = '/mnt/pesaranlab/People/Gino/code/DL-modulators/figures'
zSpecName = sprintf('%s/zSpec_Sess%03d_%s_Modulator_e%03d_e%03d_%02d_%02dHz_ROC.txt',FigDir,iSess,day,ModulatorPair(1),ModulatorPair(2),Fk(1),Fk(2))
dlmwrite(zSpecName,zSpec,'delimiter','\t')


print(sprintf('%sSpec_Sess%03d_%s_Modulator_e%03d_e%03d_%02d_%02dHz_ROC.svg',FigDir,iSess,day,ModulatorPair(1),ModulatorPair(2),Fk(1),Fk(2)),'-dsvg')



hold on

keyboard
% plot cluster boundary
% if ~isempty(sigClusterInds)
%     for iCluster = 1 : numel(sigClusterInds)
%         b = boundaries{sigClusterInds(iCluster)};
%         line(b(:,2),b(:,1),'LineWidth',2,'Color','k')
%     end
% end
% ax = gca;
% ax.XTick = [1 NX/2 NX];
% ax.XTickLabel = linspace(bn(1),bn(2),3);
% ax.YTick = 1:10:fk;
% ax.YTickLabel = fliplr(linspace(0,fk,7));
% ax.TickDir = 'out';
% xlim([1 NX])
% ylim([1 NY])
% xlabel('Time from pulse onset (ms)')
% ylabel('Frequency (Hz)')
% title(['Modulator (z = ' num2str(zThresh) ')'])
hold off