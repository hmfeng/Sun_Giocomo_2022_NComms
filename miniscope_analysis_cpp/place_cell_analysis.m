%% batch analysis of place cells
%% specify the location of datapath for all the mice need to be analyzed
load('E:\Sun_Giocomo_NComm_DataSharing\datapath.mat','datapath');
%% specify days of baseline sessions and conditioning sessions
for n = 1:length(datapath)
    cd(datapath{n})
    bsl = 1:2; %baseline
    cdn = 3:8; %conditioning
    bt = [bsl,max(cdn)+1,max(cdn)+2]; %bt: baseline and test sessions
    save('property.mat', 'bsl', 'cdn', 'bt');
end
%% Calculate spatial maps
for n = 1:length(datapath)
    cd(datapath{n})
    load ('neuronIndividualsf.mat');
    load ('behavIndividualsf.mat');
    load ('thresh.mat');
    meanFRindividuals = cell(1,length(neuronIndividualsf));
    firingrateAll = cell(1,length(neuronIndividualsf));
    countAll = cell(1,length(neuronIndividualsf));
    countTime = cell(1,length(neuronIndividualsf));
    % Calculate spatial map and mean firing rate for each individual sessions
    for k = 1:length(neuronIndividualsf)
        [firingrateAll{k},countAll{k},countTime{k}] = calculate_firing_ratemap(neuronIndividualsf{k},behavIndividualsf{k},thresh,2);
        for m = 1:size(firingrateAll{1,k},2)
            meanFRindividuals{1,k}(m,1)= sum(sum(countAll{1,k}{1,m}))/sum(sum(countTime{1,k}));
        end
    end
    save firingrateAll.mat firingrateAll countAll countTime meanFRindividuals;
end
%% Visually inspect individual spatial maps for each mouse (optional)
cell2plot = 1:10;
plot_rate_map_longitudinal(neuronIndividualsf,behavIndividualsf,firingrateAll,cell2plot,thresh,'S',2);
% Print figure into vector eps files
filepath = pwd;
neuronSelected = 142;
fig = openfig(fullfile(filepath, 'RatemapFigures', ['Cell',num2str(neuronSelected), '_ratemap.fig']));
print (fig, '-painters', '-depsc', fullfile(filepath,'RatemapFigures',['Cell',num2str(neuronSelected),'_ratemap.eps']));
% open a group of plots;
neuronSelected = [159, 199];
for ii = 1:length(neuronSelected)
    if exist(fullfile(filepath, 'RatemapFigures', ['Cell',num2str(neuronSelected(ii)), '_ratemap.fig']), 'file')
        openfig(fullfile(filepath, 'RatemapFigures', ['Cell',num2str(neuronSelected(ii)), '_ratemap.fig']));
    end
end
%% split neuron and behav acorrding to animals position at left or right chamber
%% split neuron and behav data to get NBindivLR
clearvars -except datapath
for n = 1:length(datapath)
    cd(datapath{n});
    load('neuronIndividualsf.mat');
    load('behavIndividualsf.mat');
    load('thresh.mat');
    load('property.mat');
    neuronIndivLR = cell(2,length(bt));
    behavIndivLR = cell(2,length(bt));
    midline = NaN(1, length(bt));
    for ii = 1:length(bt)
        jj = bt(ii);
    %     if ii==2 || ii==3
    %         centerdetect = 'mannual';
    %     else
    %         centerdetect = 'auto';
    %     end
        [neuronL, neuronR, behavL, behavR, autocenter] = split_neuron_behav_LR(neuronIndividualsf{jj}, behavIndividualsf{jj});
        neuronIndivLR{1,ii} = neuronL; neuronIndivLR{2,ii} = neuronR;
        behavIndivLR{1,ii} = behavL; behavIndivLR{2,ii} = behavR;
        midline(1,ii) = autocenter;
    end
    close all
    %align the edge of behavIndivLR
%     edgevalue = cellfun(@(x) min(x.position(:,1)), behavIndivLR);
%     edgeleft = max(edgevalue(1,:));
%     edgeright = max(edgevalue(2,:));
%     for ii = 1:length(behavIndivLR)
%         edgeleftdiff = min(behavIndivLR{1,ii}.position(:,1)) - edgeleft;
%         behavIndivLR{1,ii}.position(:,1) = behavIndivLR{1,ii}.position(:,1) - edgeleftdiff;
%         neuronIndivLR{1,ii}.pos(:,1) = neuronIndivLR{1,ii}.pos(:,1) - edgeleftdiff;
%         idxb = find(behavIndivLR{1,ii}.position(:,1) > 28);
%         behavIndivLR{1,ii}.position(idxb,1) = 28;
%         idxn = find(neuronIndivLR{1,ii}.pos(:,1) > 28);
%         neuronIndivLR{1,ii}.pos(idxn,1) = 28;
%         edgerightdiff = min(behavIndivLR{2,ii}.position(:,1)) - edgeleft;
%         behavIndivLR{2,ii}.position(:,1) = behavIndivLR{2,ii}.position(:,1) - edgerightdiff;
%         neuronIndivLR{2,ii}.pos(:,1) = neuronIndivLR{2,ii}.pos(:,1) - edgerightdiff;
%         idxb = find(behavIndivLR{2,ii}.position(:,1) > 28);
%         behavIndivLR{2,ii}.position(idxb,1) = 28;
%         idxn = find(neuronIndivLR{2,ii}.pos(:,1) > 28);
%         neuronIndivLR{2,ii}.pos(idxn,1) = 28;
%     end
    %get the raw ratemap
    meanfrLR = cell(2,length(bt));
    for k = 1:size(neuronIndivLR,2)
        [frLR{1,k},countLR{1,k},timeLR{1,k}] = calculate_subset_ratemap(neuronIndivLR{1,k},behavIndivLR{1,k},thresh,2);
        [frLR{2,k},countLR{2,k},timeLR{2,k}] = calculate_subset_ratemap(neuronIndivLR{2,k},behavIndivLR{2,k},thresh,2);
        for m = 1:size(frLR{1,k},2)
            meanfrLR{1,k}(m,1)= sum(sum(countLR{1,k}{1,m}))/sum(sum(timeLR{1,k}));
            meanfrLR{2,k}(m,1)= sum(sum(countLR{2,k}{1,m}))/sum(sum(timeLR{2,k}));
        end
    end
    save('NBindivLR.mat', 'neuronIndivLR', 'behavIndivLR', 'frLR', 'countLR',...
        'timeLR', 'meanfrLR','midline', 'v-7.3');
    clearvars -except datapath n
end
% check some plot if you want
load('NBindivLR.mat');
load('thresh.mat');
cell2plot = 95;
plot_rate_map_longitudinal(neuronIndivLR(1,:),behavIndivLR(1,:),frLR(1,:),cell2plot,thresh,'S',2,false);
plot_rate_map_longitudinal(neuronIndivLR(2,:),behavIndivLR(2,:),frLR(2,:),cell2plot,thresh,'S',2,false);
%% Calculate all the metrics and shuffle place cells by using left and right matched data
clearvars -except datapath
for k = 1:length(datapath)
    cd(datapath{k});
    load('property.mat');
    load('NBindivLR.mat', 'neuronIndivLR', 'behavIndivLR')
    load('thresh.mat');
    niter = 20;
    metrics = cell(2, length(bt), niter);
    metrics_avg = cell(2, length(bt));
    for n = 1:niter
        [neuronIndivLRm, behavIndivLRm] = match_neuron_behav(neuronIndivLR, behavIndivLR, 'LR');
        parfor ii = 1:length(bt)
            neuronL = neuronIndivLRm{1,ii};
            behavL = behavIndivLRm{1,ii};
            [~,TinfoL,~] = permutingSpike_LR(neuronL,behavL,thresh);
            metrics{1,ii,n} = table2array(TinfoL); %defined by bit/sec
        end
        parfor jj = 1:length(bt)
            neuronR = neuronIndivLRm{2,jj};
            behavR = behavIndivLRm{2,jj};
            [~,TinfoR,~] = permutingSpike_LR(neuronR,behavR,thresh);
            metrics{2,jj,n} = table2array(TinfoR); %defined by bit/sec
        end
    end
    %average the value for each iteration
    for ii = 1:size(metrics,1)
        for jj = 1:size(metrics,2)
            metrics_each = squeeze(metrics(ii,jj,:));
            metrics_each = cat(3, metrics_each{:});
            metrics_avg{ii,jj} = array2table(nanmean(metrics_each,3),...
                'VariableNames',{'neuron','bitpsec','bitpsec_pcidx','bitpspike','bitpspike_pcidx','meanfr','peakfr'});
        end
    end
    save('Metrics_match.mat', 'metrics', 'metrics_avg', '-v7.3');
    clear metrics metrics_avg;
end
for k = 1:length(datapath)
    cd(datapath{k});
    load('Metrics_match.mat', 'metrics_avg')
    placecellsLRm = cellfun(@(x) find(x.meanfr >= 0.1 & x.bitpspike_pcidx >= 0.3), metrics_avg,...
        'Uniformoutput', false);
    save('PlaceCellsLR.mat', 'placecellsLRm')
end
%% Calculate all the metrics and shuffle place cells by using longitudinally matched data
% for k = 1:length(datapath)
%     cd(datapath{k});
%     load('property.mat');
%     load('NBindivLR.mat', 'neuronIndivLR', 'behavIndivLR')
%     load('thresh.mat');
%     niter = 20;
%     metrics = cell(2, length(bt), niter);
%     metrics_avg = cell(2, length(bt));
%     for n = 1:niter
%         [neuronIndivLRm, behavIndivLRm] = match_neuron_behav(neuronIndivLR, behavIndivLR, 'longitudinal_full');
%         parfor ii = 1:length(bt)
%             neuronL = neuronIndivLRm{1,ii};
%             behavL = behavIndivLRm{1,ii};
%             [~,TinfoL,~] = permutingSpike_LR(neuronL,behavL,thresh);
%             metrics{1,ii,n} = table2array(TinfoL); %defined by bit/sec
%         end
%         parfor jj = 1:length(bt)
%             neuronR = neuronIndivLRm{2,jj};
%             behavR = behavIndivLRm{2,jj};
%             [~,TinfoR,~] = permutingSpike_LR(neuronR,behavR,thresh);
%             metrics{2,jj,n} = table2array(TinfoR); %defined by bit/sec
%         end
%     end
%     %average the value for each iteration
%     for ii = 1:size(metrics,1)
%         for jj = 1:size(metrics,2)
%             metrics_each = squeeze(metrics(ii,jj,:));
%             metrics_each = cat(3, metrics_each{:});
%             metrics_avg{ii,jj} = array2table(nanmean(metrics_each,3),...
%                 'VariableNames',{'neuron','bitpsec','bitpsec_pcidx','bitpspike','bitpspike_pcidx','meanfr','peakfr'});
%         end
%     end
%     save('Metrics_matchL.mat', 'metrics', 'metrics_avg', '-v7.3');
%     clear metrics metrics_avg;
% end
% load(path_all)
% for k = 1:length(datapath)
%     cd(datapath{k});
%     load('Metrics_matchL.mat', 'metrics_avg')
%     placecellsLRml = cellfun(@(x) find(x.meanfr >= 0.1 & x.bitpspike_pcidx >= 0.3), metrics_avg,...
%         'Uniformoutput', false);
%     save('PlaceCellsLR.mat', 'placecellsLRml','-append')
% end
%% shuffle to get place cells using the same simulated trajectories across multiple different sessions
% % use the same trajectory to match behavior across multiple sessions
% for n = 1:length(datapath)
%     cd(datapath{n});
%     load('NBindivLR.mat','neuronIndivLR','behavIndivLR');
%     load('behavior.mat','behavcover')
% %     datalength = cellfun(@(x) length(x.position),behavIndivLR);
% %     datalength(:,1) = inf; 
% %     [~,idxL] = min(datalength(1,:),[],2);
% %     [~,idxR] = min(datalength(2,:),[],2);
%     [~,idxL] = min(behavcover(1,:),[],2);
%     [~,idxR] = min(behavcover(2,:),[],2);
%     neuronIndivLRmt = cell(size(neuronIndivLR,1),size(neuronIndivLR,2));
%     behavIndivLRmt = cell(size(behavIndivLR,1),size(behavIndivLR,2));
%     % match the left compartment
%     neuron2 = neuronIndivLR{1,idxL};
%     behav2 = behavIndivLR{1,idxL};
%     for ii = 1:length(neuronIndivLR)
%         if ii == idxL
%             neuronIndivLRmt{1,ii} = neuron2;
%             behavIndivLRmt{1,ii} = behav2;
%         else
%             neuron1 = neuronIndivLR{1,ii};
%             behav1 = behavIndivLR{1,ii};
%             [neuronc, behavc, ~, ~] = downsample_match_trajectory(neuron1,behav1,neuron2,behav2);
%             neuronIndivLRmt{1,ii} = neuronc;
%             behavIndivLRmt{1,ii} = behavc;
%         end
%     end
%     % match the right compartment
%     neuron2 = neuronIndivLR{2,idxR};
%     behav2 = behavIndivLR{2,idxR};
%     for ii = 1:length(neuronIndivLR)
%         if ii == idxR
%             neuronIndivLRmt{2,ii} = neuron2;
%             behavIndivLRmt{2,ii} = behav2;
%         else
%             neuron1 = neuronIndivLR{2,ii};
%             behav1 = behavIndivLR{2,ii};
%             [neuronc, behavc, ~, ~] = downsample_match_trajectory(neuron1,behav1,neuron2,behav2);
%             neuronIndivLRmt{2,ii} = neuronc;
%             behavIndivLRmt{2,ii} = behavc;
%         end
%     end
%     load('thresh.mat');
%     meanfrLR = cell(2,length(neuronIndivLRmt));
%     for k = 1:size(neuronIndivLR,2)
%         [frLR{1,k},countLR{1,k},timeLR{1,k}] = calculate_subset_ratemap(neuronIndivLRmt{1,k},behavIndivLRmt{1,k},thresh,2);
%         [frLR{2,k},countLR{2,k},timeLR{2,k}] = calculate_subset_ratemap(neuronIndivLRmt{2,k},behavIndivLRmt{2,k},thresh,2);
%         for m = 1:size(frLR{1,k},2)
%             meanfrLR{1,k}(m,1)= sum(sum(countLR{1,k}{1,m}))/sum(sum(timeLR{1,k}));
%             meanfrLR{2,k}(m,1)= sum(sum(countLR{2,k}{1,m}))/sum(sum(timeLR{2,k}));
%         end
%     end
%     save('NBindivLRmt.mat', 'neuronIndivLRmt', 'behavIndivLRmt', 'frLR', 'countLR',...
%         'timeLR', 'meanfrLR', '-v7.3');
%     clearvars -except datapath n
% end
% %check for some rate map results
% cell2plot = 12;
% figure;
% ratemap = filter2DMatrices(frLR{1,4}{cell2plot},1);
% pcolor(ratemap);
% shading flat
% axis image
% colormap jet
% caxis([0 2])
% %shuffle NBindivLRmt to get place cells and metrics
% for n = 1:length(datapath)
%     cd(datapath{n});
%     load('NBindivLRmt.mat','neuronIndivLRmt','behavIndivLRmt');
%     load('thresh.mat');
%     metrics_mt = cell(2, length(neuronIndivLRmt));
%     parfor ii = 1:length(neuronIndivLRmt)
%         neuronL = neuronIndivLRmt{1,ii};
%         behavL = behavIndivLRmt{1,ii};
%         [~,TinfoL,~] = permutingSpike_LR(neuronL,behavL,thresh);
%         metrics_mt{1,ii} = TinfoL;
%     end
%     parfor jj = 1:length(neuronIndivLRmt)
%         neuronR = neuronIndivLRmt{2,jj};
%         behavR = behavIndivLRmt{2,jj};
%         [~,TinfoR,~] = permutingSpike_LR(neuronR,behavR,thresh);
%         metrics_mt{2,jj} = TinfoR;
%     end
%     placecellsLRmt = cellfun(@(x) find(x.meanfr >= 0.1 & x.bitpspike_pcidx == 1),...
%         metrics_mt, 'Uni', false);
%     save('Metrics_match.mat','metrics_mt','-append');
%     save('PlaceCellsLR.mat', 'placecellsLRmt', '-append');
% end
% for n = 1:length(datapath)
%     cd(datapath{n});
%     load('Metrics_match.mat','metrics_mt');
%     load('NBindivLR.mat','meanfrLR');
%     A = cellfun(@(x) find(x.bitpspike_pcidx == 1),metrics_mt, 'Uni', false);
%     B = cellfun(@(x) find(x >= 0.1),meanfrLR, 'Uni', false);
%     placecellsLRmt = cellfun(@(x,y) intersect(x,y), A,B, 'Uni', false);
%     save('PlaceCellsLR.mat', 'placecellsLRmt', '-append');
% end
%% Plot the proportion of place cells over time
pcnum_meth = struct;
pcnum_meth.p = []; pcnum_meth.np = [];
for ii = 1:length(datapath)
    cd(datapath{ii,1});
    load('CPPscore.mat');
    load('PlaceCellsLR.mat', 'placecellsLRm')
    load('property.mat');
    load('thresh.mat')
    ti = find(bt > max(cdn));
    pcnumber = cell2mat(cellfun(@(x) length(x)/length(thresh),placecellsLRm,'uni',0));
    if CPPtimen(2, 1) > CPPtimen(2, 2)
        pcnum_meth.p = padconcatenation(pcnum_meth.p, [pcnumber(1,bsl(end-1:end)),pcnumber(1,ti)],1);
        pcnum_meth.np = padconcatenation(pcnum_meth.np, [pcnumber(2,bsl(end-1:end)),pcnumber(2,ti)],1);
    else
        pcnum_meth.p = padconcatenation(pcnum_meth.p, [pcnumber(2,bsl(end-1:end)),pcnumber(2,ti)],1);
        pcnum_meth.np = padconcatenation(pcnum_meth.np, [pcnumber(1,bsl(end-1:end)),pcnumber(1,ti)],1);
    end
end
% save('pcnum.mat','pcnum_meth','-append')
[h,p1,~,t1] = ttest(pcnum_meth.p(:,2),pcnum_meth.p(:,3));
[h,p2,~,t2] = ttest(pcnum_meth.p(:,2),pcnum_meth.p(:,4));
[h,p3,~,t3] = ttest(pcnum_meth.np(:,2),pcnum_meth.np(:,3));
[h,p4,~,t4] = ttest(pcnum_meth.np(:,2),pcnum_meth.np(:,4));
figure;
subplot(2,1,1)
plot(pcnum_meth.p', '-o', 'Color', 'k'); 
xlim([0.5 4.5]);
ylim([0 0.8])
xticks([1 2 3 4])
xticklabels({'pre-bsl','bsl','test1','test2'})
ylabel('Proportion of place cells')
xlabel('Saline-paired side')
text(2.5,0.8,['P=',num2str(p1)])
text(2.5,0.7,['t=',num2str(t1.tstat)])
text(4,0.8,['P=',num2str(p2)])
text(4,0.7,['t=',num2str(t2.tstat)])
subplot(2,1,2)
plot(pcnum_meth.np', '-o', 'Color','k');
xlim([0.5 4.5]);
ylim([0 0.8])
xticks([1 2 3 4])
xticklabels({'pre-bsl','bsl','test1','test2'})
ylabel('Proportion of place cells')
xlabel('MA-paired side')
text(2.5,0.8,['P=',num2str(p3)])
text(2.5,0.7,['t=',num2str(t3.tstat)])
text(4,0.8,['P=',num2str(p4)])
text(4,0.7,['t=',num2str(t4.tstat)])

%% Identify functional cell types, disPCp, disPCnp, aPCp, aPCnp
clearvars -except datapath
for k = 1:length(datapath)
    cd(datapath{k});
    load('PlaceCellsLR.mat', 'placecellsLRm')
    load('CPPscore.mat','CPPtimen')
    load('property.mat')
    dispcm = struct;
    if CPPtimen(2, 1) > CPPtimen(2, 2)
        lia1 = ~ismember(placecellsLRm{1,max(bsl)}, union(placecellsLRm{1,max(bsl)+1},placecellsLRm{1,max(bsl)+2}));
        dispcm.p = placecellsLRm{1,max(bsl)}(lia1);
        lia2 = ~ismember(placecellsLRm{2,max(bsl)}, union(placecellsLRm{2,max(bsl)+1},placecellsLRm{2,max(bsl)+2}));
        dispcm.np = placecellsLRm{2,max(bsl)}(lia2);
    else
        lia1 = ~ismember(placecellsLRm{2,max(bsl)}, union(placecellsLRm{2,max(bsl)+1},placecellsLRm{2,max(bsl)+2}));
        dispcm.p = placecellsLRm{2,max(bsl)}(lia1);
        lia2 = ~ismember(placecellsLRm{1,max(bsl)}, union(placecellsLRm{1,max(bsl)+1},placecellsLRm{1,max(bsl)+2}));
        dispcm.np = placecellsLRm{1,max(bsl)}(lia2);
    end
    save('PlaceCellsLR.mat', 'dispcm', '-append')
end
% neurons that become place cell after cdn
for k = 1:length(datapath)
    cd(datapath{k});
    load('PlaceCellsLR.mat', 'placecellsLRm')
    load('CPPscore.mat','CPPtimen')
    load('property.mat')
    appcm = struct;
    if CPPtimen(2, 1) > CPPtimen(2, 2)
        ppc = intersect(placecellsLRm{1,max(bsl)+1},placecellsLRm{1,max(bsl)+2});
        lia1 = ~ismember(ppc, placecellsLRm{1,max(bsl)});
        appcm.p = ppc(lia1);
        nppc = intersect(placecellsLRm{2,max(bsl)+1},placecellsLRm{2,max(bsl)+2});
        lia2 = ~ismember(nppc, placecellsLRm{2,max(bsl)});
        appcm.np = nppc(lia2);
    else
        ppc = intersect(placecellsLRm{2,max(bsl)+1},placecellsLRm{2,max(bsl)+2});
        lia1 = ~ismember(ppc, placecellsLRm{2,max(bsl)});
        appcm.p = ppc(lia1);
        nppc = intersect(placecellsLRm{1,max(bsl)+1},placecellsLRm{1,max(bsl)+2});
        lia2 = ~ismember(nppc, placecellsLRm{1,max(bsl)});
        appcm.np = nppc(lia2);
    end
    save('PlaceCellsLR.mat', 'appcm', '-append')
end
%% Plot the percentage of functional cell types
dp=[]; dnp=[]; ap= []; anp = [];
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('PlaceCellsLR.mat','dispcm','appcm')
    load('thresh.mat')
    dp = [dp;length(dispcm.p)/length(thresh)];
    dnp = [dnp;length(dispcm.np)/length(thresh)];
    ap = [ap;length(appcm.p)/length(thresh)];
    anp = [anp;length(appcm.np)/length(thresh)];
end
cellp_all = padcat(dp,dnp,ap,anp);
[h,p1,~,t1] = ttest(cellp_all(:,1),cellp_all(:,2));
[h,p2,~,t2] = ttest(cellp_all(:,1),cellp_all(:,3));
[h,p3,~,t3] = ttest(cellp_all(:,1),cellp_all(:,4));
[h,p4,~,t4] = ttest(cellp_all(:,3),cellp_all(:,4));

data = cellp_all;
figure;
hold on;
colors = [128 130 133;0 0 0;255 100 100;255 0 0]/255; % ctrl MA
h = boxplot(data,'colors',colors,'symbol', '.');
xlim([0.5 size(data,2)+0.5])
xlabels = {'disPCp','disPCnp','aPCp','aPCnp'};
set(gca,'Xtick',1:size(data,2))
set(gca,'XtickLabel',xlabels,'FontSize',12,'FontName','Arial')
ylabel('Proportion of total neurons','FontSize',12,'FontName','Arial')
% Alter linestyle
h2 = findobj(gca,'Tag','Box');
set(h2,{'linew'},{1.5})
h1 = findobj(gca,'tag','Median');
set(h1,{'linew'},{2.5})
h1 = findobj(gca,'tag','Upper Whisker');
set(h1,'LineStyle','-')
h1 = findobj(gca,'tag','Lower Whisker');
set(h1,'LineStyle','-')
box off
plot([1,2],data(:,1:2)','-o','Color','k')
plot([3,4],data(:,3:4)','-o','Color','k')
title('MA CPP')
ylim([0.02,0.26])
text(2,0.2, ['P=',num2str(p1)]);
text(2,0.18, ['t=',num2str(t1.tstat)]);
text(3,0.2, ['P=',num2str(p2)]);
text(3,0.18, ['t=',num2str(t2.tstat)]);
text(4,0.2, ['P=',num2str(p3)]);
text(4,0.18, ['t=',num2str(t3.tstat)]);
text(4,0.25, ['P=',num2str(p4)]);
text(4,0.23, ['t=',num2str(t4.tstat)]);

%% Compute the correlation of spatial maps between the left and right CPP compartments
clearvars -except datapath
for n = 1:length(datapath)
    cd(datapath{n});
    load('NBindivLR.mat','neuronIndivLR','behavIndivLR');
    load('thresh.mat');
    load('property.mat');
    corr_LRm_all = [];
    for ii = 1:50
        [ratemapt_LRm, frLRm, countLRm, timeLRm] = match_smooth_trim_ratemap(neuronIndivLR, behavIndivLR, thresh, 'LR');
        [corr_LRm, avg_corr_LRm] = get_corr_matrix([], ratemapt_LRm(1,:), ratemapt_LRm(2,:), 'LR');
        corr_LRm_all = cat(3, corr_LRm_all, avg_corr_LRm);
    end
    close all
    save ('Corr_match.mat', 'corr_LRm_all');
end
%% Linear regression between CorrDiff and CPP score
clearvars -except datapath
% get the CPP score from all mice
CPPmeth = [];
for ii = 1:length(datapath)
    cd(datapath{ii,1});
    load('CPPscore.mat');
    temp = [CPPscore1, CPPscore2];
    CPPmeth = [CPPmeth; temp];
end
% perform the regression
dcorr = NaN(length(datapath),2);
for ii = 1:length(datapath)
    cd(datapath{ii})
    load('Corr_match.mat','corr_LRm_all')
    load('PlaceCellsLR.mat','dispcm')
    corrLRm = atanh(mean(corr_LRm_all,3));
    corrLRm_pc = mean(corrLRm(dispcm.p,:));
    corrLRm_pc = corrLRm_pc(:, end-2:end);
    dcorr(ii,1) = corrLRm_pc(1) - corrLRm_pc(2);
    dcorr(ii,2) = corrLRm_pc(1) - corrLRm_pc(3);
end
cfplot_cpp(dcorr(:), CPPmeth(:), 'correlation', 'CPP_corr_meth');
