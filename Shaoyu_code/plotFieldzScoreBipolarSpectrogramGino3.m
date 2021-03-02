function plotFieldzScoreBipolarSpectrogramClusterCorrected(data,hitIndx,missIndx,bn,tapers,fs,dn,fk,pad,nPerm,alpha)

Spec = []; TestSpec = []; DiffSpec = []; AvDiffSpec = [];
[Spec,f] = tfspec(data(hitIndx,:),tapers,fs,dn,fk,pad);
nWin = size(Spec,2);
TestSpec = sq(sum(log(Spec))./numel(hitIndx));

DiffSpec = tfspec(data(missIndx,:),tapers,fs,dn,fk,pad);
AvDiffSpec =  sq(sum(log(DiffSpec)))./numel(missIndx);

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
% shuffleDperm = permute(shuffleDperm,[2 1 3]);
% %zSpec = sign(D).*norminv(1-pSpec,0,1);
% 
% % get distribution of largest signifciant cluster sizes
% cluster = zeros(1,1e3);
% for iPerm = 1 : nPerm
%     if(mod(iPerm,1000) == 0)
%         %feedback
%         disp(['Permutation ' num2str(iPerm) ' of ' num2str(nPerm)])
%     end
%     
%     %shuffleDperm  is the null dist
%     tmpD = sq(shuffleDperm(iPerm,:));
%     tmpD = tmpD';
%     
%     %convert tmpD current shuffle to a zscore
%     tmppSpec = zeros(nWin,numel(f));
%     for iWin = 1:nWin
%         for iF = 1:size(tmpD,2)
%             tmppSpec(iWin,iF) = length(find(shuffleDperm(:,iWin,iF)>tmpD(iWin,iF)))./nPerm;
%         end
%     end
%     tmppSpec(tmppSpec==0) = 1./nPerm;  tmppSpec(tmppSpec==1) = (nPerm-1)./nPerm;
%     
%     %convert zscore to sig/not sig(1/0)
%     for itmp =1:size(tmppSpec,1)
%         tmppSpec(itmp,tmppSpec(itmp,:) < (1-alpha)) = 0;
%         tmppSpec(itmp,tmppSpec(itmp,:) >= (1-alpha)) = 1;
%     end
%     
%     %cluster
%     [L,n] = bwlabel(tmppSpec);
%     bins = hist(reshape(L,size(L,1)*size(L,2),1),n+1);
%     if(length(bins) > 1)
%         cluster(iPerm) = max(bins(2:end));
%     else
%         cluster(iPerm) = 0;
%     end
% end
% 
% cluster = sort(cluster);
% minClustSize = cluster(round(0.99*numel(cluster))+1);
% %minClustSize = prctile(cluster,95);
% 
% % find significant pixels
% zSpec1 = flipud(zSpec');
% 
% if isequal(alpha,0.05)
%     zThresh = 1.96; % 5%
% elseif isequal(alpha,0.01)
%     zThresh = 2.58; % 1%
% end
% 
% L = [];
% bins = [];
% [L,n] = bwlabel(abs(zSpec1)>zThresh);
% bins = hist(reshape(L,size(L,1)*size(L,2),1),n+1);
% 
% % find the clusters that need to be kept
% keepClustID = [];
% keepClustID = find(bins(2:end)>minClustSize);
% 
% % convert to binary image
% LL = [];
% LL = ismember(L,keepClustID);
% 
% % find the cluster boundaries
% boundaries = [];
% %[boundaries,~] = bwboundaries(LL,'noholes');
% [boundaries,~] = bwboundaries(LL,'holes');
% 
% zSpec_cluster = [];
% zSpec_cluster_sum = [];
% sigClusterInds = [];
% if ~isempty(boundaries)
%     %%%%% calculate cluster-level z statistic %%%%%
%     for iCluster = 1 : numel(boundaries)
%         b = boundaries{iCluster};
%         for iB = 1 : size(b,1)
%             zSpec_cluster{iCluster}(iB) = zSpec1(b(iB,1),b(iB,2));
%         end
%         zSpec_cluster_sum(iCluster) = sum(zSpec_cluster{iCluster});
%     end
%     
%     clusterInd_pos = [];
%     clusterInd_pos = find(zSpec_cluster_sum > 0);
%     clusterInd_neg = [];
%     clusterInd_neg = find(zSpec_cluster_sum < 0);
%     
%     maxClusterInd = [];
%     [~,maxClusterInd] = max(abs(zSpec_cluster_sum));
%     
%     % calculate Monte Carlo p value
%     disp('Calculating Monte Carlo p value');
%     clusterInd = [];
%     clusterInd = setdiff(1 : numel(boundaries), maxClusterInd);
%     prob = [];
%     for iCluster = 1 : numel(clusterInd)
%         i = clusterInd(iCluster);
%         [prob(iCluster),~,~] = calcClusterCorrectPvalue(abs(zSpec_cluster{maxClusterInd}'),abs(zSpec_cluster{i}'),nPerm);
%     end
%     prob(prob==0) = 1./nPerm;
%     pcorrect = [];
%     
%     if ~isempty(prob)
%         [~,pcorrect,~] = false_discovery_rate(1-prob);
%         sigClusterInds = [clusterInd(pcorrect<=alpha/2) maxClusterInd];
%     else
%         sigClusterInds = maxClusterInd;
%     end
%     
%     
%     %     if ~isempty(clusterInd_pos)
%     %         maxClusterInd_pos = [];
%     %         [~,maxClusterInd_pos] = max(zSpec_cluster_sum(clusterInd_pos));
%     %         maxClusterInd_pos = clusterInd_pos(maxClusterInd_pos);
%     %
%     %         if numel(clusterInd_pos) > 1
%     %             % calculate Monte Carlo p value
%     %             disp('Calculating Monte Carlo p value');
%     %             clusters_pos = [];
%     %             clusters_pos = setdiff(clusterInd_pos, maxClusterInd);
%     %
%     %             p_pos = [];
%     %             for iCluster = 1 : numel(clusters_pos)
%     %                 i = clusters_pos(iCluster);
%     %                 [p_pos(iCluster),~,~] = calcClusterCorrectPvalue(abs(zSpec_cluster{i}'),abs(zSpec_cluster{maxClusterInd}'),nPerm);
%     %             end
%     %             p_pos(p_pos==0) = 1./nPerm;
%     %
%     %             pcorrect_pos = [];
%     %             [~,pcorrect_pos,~] = false_discovery_rate(p_pos);
%     %             sigClusterInds_pos = [];
%     %             sigClusterInds_pos = [clusters_pos(pcorrect_pos<pThresh) maxClusterInd];
%     %         else
%     %             sigClusterInds_pos = [];
%     %             sigClusterInds_pos = maxClusterInd_pos;
%     %         end
%     %     else
%     %         disp(['No positive clusters'])
%     %         sigClusterInds_pos = [];
%     %     end
%     %
%     %     if ~isempty(clusterInd_neg)
%     %         maxClusterInd_neg = [];
%     %         [~,maxClusterInd_neg] = max(abs(zSpec_cluster_sum(clusterInd_neg)));
%     %         maxClusterInd_neg = clusterInd_neg(maxClusterInd_neg);
%     %
%     %         if numel(clusterInd_neg) > 1 ||
%     %             % calculate Monte Carlo p value
%     %             disp('Calculating Monte Carlo p value');
%     %             clusters_neg = [];
%     %             clusters_neg = setdiff(clusterInd_neg, maxClusterInd);
%     %
%     %             p_neg = [];
%     %             for iCluster = 1 : numel(clusters_neg)
%     %                 i = clusters_neg(iCluster);
%     %                 [p_neg(iCluster),~,~] = calcClusterCorrectPvalue(abs(zSpec_cluster{i}'),abs(zSpec_cluster{maxClusterInd_neg}'),nPerm);
%     %             end
%     %             p_neg(p_neg==0) = 1./nPerm;
%     %
%     %             pcorrect_neg = [];
%     %             [~,pcorrect_neg,~] = false_discovery_rate(p_neg);
%     %             sigClusterInds_neg = [];
%     %             sigClusterInds_neg = [clusters_neg(pcorrect_neg<pThresh) maxClusterInd_neg];
%     %         else
%     %             sigClusterInds_neg = [];
%     %             sigClusterInds_neg = maxClusterInd_neg;
%     %         end
%     %     else
%     %         disp(['No negative clusters'])
%     %         sigClusterInds_neg = [];
%     %     end
%     %     sigClusterInds = [sigClusterInds_pos sigClusterInds_neg];
% end

%% plot
%clim = [-ceil(max(max(abs(zSpec)))) ceil(max(max(abs(zSpec))))];
clim = [-4 4];
%
NX = size(zSpec,1);
NY = size(zSpec,2);
%X = linspace(bn(1),bn(2),NX);
%Y = linspace(0,fk,NY);
%imagesc(X,Y,flipud(zSpec'),clim)
imagesc(flipud(zSpec'),clim)
%colormap(gca,polarmap(winter,1))
colormap(gca,polarmap(jet,1.5))
%colormap(gca,polarmap(polarmap,0.5))

%cmap = cbrewer('div','PuOr',9);
%colormap(cmap)

c = colorbar;
c.Label.String = 'Hit - Miss (Z-score)';
hold on

% plot cluster boundary
if ~isempty(sigClusterInds)
    for iCluster = 1 : numel(sigClusterInds)
        b = boundaries{sigClusterInds(iCluster)};
        line(b(:,2),b(:,1),'LineWidth',2,'Color','k')
    end
end
ax = gca;
ax.XTick = [1 NX/2 NX];
ax.XTickLabel = linspace(bn(1),bn(2),3);
ax.YTick = 1:10:fk;
ax.YTickLabel = fliplr(linspace(0,fk,7));
ax.TickDir = 'out';
xlim([1 NX])
ylim([1 NY])
xlabel('Time from pulse onset (ms)')
ylabel('Frequency (Hz)')
title(['Modulator (z = ' num2str(zThresh) ')'])
hold off