function [ratemapt_LRm, frLRm, countLRm, timeLRm] = match_smooth_trim_ratemap(neuronIndivLR, behavIndivLR, thresh, condition)
% This function performs multiple steps for ratemap analysis
% 1. downsample data to match bahavior coverages across left and right cpp
% experment chambers or across different sessions;
% 2. calculate raw firing ratemap by using matched neuron and behav data;
% 3. smooth the raw firing ratemap;
% 4. trim the smoothed raw firing ratemap in order to prepare for
% correlation computation.
% Input args:
%  -neuronIndivLR, splitted neuron data from left and right cpp chambers
%  -behavIndivLR, splitted behav data from left and right cpp chamebrs
% Output:
%  -ratemapt_LRm: trimed, smoothed, and matched ratemaps
%  -frLRm: raw firingrate
%  -countLRm: spike counts for matched datasets
%  -timeLRm: occupancy time for matched datasets
% note: this function contains customized functions naming as follows
%  -downsample_match_position
%  -calculate_subset_ratemap
%  -smooth_trim_ratemap
% Yanjun Sun, Stanford University, 9/16/2019

if ~exist('condition', 'var') || isempty(condition)
    condition = 'xsession';
end

switch condition
    case {'LR', 'left right', 'Left Right', 'left and right', 'Left and Right'}
        % preallocate
        neuronIndivLRm = cell(size(neuronIndivLR,1), size(neuronIndivLR,2));
        behavIndivLRm = cell(size(behavIndivLR,1), size(behavIndivLR,2));
        frLRm = cell(size(neuronIndivLR,1), size(neuronIndivLR,2));
        countLRm = cell(size(neuronIndivLR,1), size(neuronIndivLR,2));
        timeLRm = cell(size(neuronIndivLR,1), size(neuronIndivLR,2));
        % downsample the data points to match behaviors of left and right chambers
        for ii = 1:length(neuronIndivLR);
            [neuronLm, behavLm, neuronRm, behavRm] = downsample_match_position(neuronIndivLR{1,ii},behavIndivLR{1,ii},neuronIndivLR{2,ii},behavIndivLR{2,ii});
            neuronIndivLRm{1,ii} = neuronLm; behavIndivLRm{1,ii} = behavLm;
            neuronIndivLRm{2,ii} = neuronRm; behavIndivLRm{2,ii} = behavRm;
        end
        for k = 1:size(neuronIndivLR,2);
            [frLRm{1,k},countLRm{1,k},timeLRm{1,k}] = calculate_subset_ratemap(neuronIndivLRm{1,k},behavIndivLRm{1,k},thresh,2);
            [frLRm{2,k},countLRm{2,k},timeLRm{2,k}] = calculate_subset_ratemap(neuronIndivLRm{2,k},behavIndivLRm{2,k},thresh,2);
%             for m = 1:size(frLR{1,k},2);
%                 meanfrLRm{1,k}(m,1)= sum(sum(countLR{1,k}{1,m}))/sum(sum(timeLR{1,k}));
%                 meanfrLRm{2,k}(m,1)= sum(sum(countLR{2,k}{1,m}))/sum(sum(timeLR{2,k}));
%             end;
        end;
        % smooth and trim frLRm
        ratemapt_LRm = cell(2,length(neuronIndivLR));
        for ii = 1:length(neuronIndivLR)
            [~, ratemapt] = smooth_trim_ratemap(frLRm(:,ii)', 'top right');
            ratemapt_LRm(:,ii) = ratemapt';
        end
%         % peak firing rate based on smoothed rate maps
%         peakfrLRm = cell(2,length(neuronIndivLR));
%         for k = 1:length(neuronIndivLR);
%             for m = 1:size(ratemapt_LRm{1,k},2);
%                 peakfrLRm{1,k}(m,1)= max(max(ratemapt_LRm{1,k}{1,m}));
%                 peakfrLRm{2,k}(m,1)= max(max(ratemapt_LRm{2,k}{1,m}));
%             end;
%         end;

    case{'xsession', 'cross-session', 'x-session', 'longitudinal'}
        % downsample the data points to match behaviors of longitudinal recordings
        % data will be stored in a cell array, with (1,2) stores matched data of
        % session 1 from 1 and 2 matching, (2,1) stores matched data of session 2
        % from 1 and 2 matching, so on and so forth. dimension will be how
        % many subsegments of the data are required
        dimension = size(neuronIndivLR,1);
        if dimension > 1
            neuronIndivLRm = cell(length(neuronIndivLR),length(neuronIndivLR), dimension);
            behavIndivLRm = cell(length(behavIndivLR),length(behavIndivLR), dimension);
            frLRm = cell(length(behavIndivLR),length(behavIndivLR), dimension);
            countLRm = cell(length(behavIndivLR),length(behavIndivLR), dimension);
            timeLRm = cell(length(behavIndivLR),length(behavIndivLR), dimension);
            for d = 1:dimension
                for ii = 1:length(neuronIndivLR);
                    for jj = ii:length(neuronIndivLR);
                    [neuronL1, behavL1, neuronL2, behavL2] = downsample_match_position(neuronIndivLR{d,ii},behavIndivLR{d,ii},neuronIndivLR{d,jj},behavIndivLR{d,jj});
                    neuronIndivLRm{ii,jj,d} = neuronL1; behavIndivLRm{ii,jj,d} = behavL1;
                    neuronIndivLRm{jj,ii,d} = neuronL2; behavIndivLRm{jj,ii,d} = behavL2;
                    end
                end
            end
            % smooth and trim frLRm
            for d = 1:dimension
                for ii = 1:size(neuronIndivLRm,1)
                    for jj = 1:size(neuronIndivLRm,2)
                    [frLRm{ii,jj,d},countLRm{ii,jj,d},timeLRm{ii,jj,d}] = calculate_subset_ratemap(neuronIndivLRm{ii,jj,d},behavIndivLRm{ii,jj,d},thresh,2);
                %     for m = 1:size(frLRm{ii,jj,d},2)
                %         meanfrLRm{ii,jj}(m,1)= sum(sum(countLm{ii,jj}{1,m}))/sum(sum(timeLm{ii,jj}));
                %     end
                    end
                end
            end
            % smooth and trim ratemap
            ratemapt_LRm = cell(length(neuronIndivLRm),length(neuronIndivLRm),dimension);
            for d = 1:dimension
                for ii = 1:length(neuronIndivLRm);
                    for jj = 1:length(neuronIndivLRm);
                        frMatch = {frLRm{ii,jj,d},frLRm{jj,ii,d}};
                        [~, ratemapt] = smooth_trim_ratemap(frMatch, 'top right');
                        ratemapt_LRm{ii,jj,d} = ratemapt{:,1};
                        ratemapt_LRm{jj,ii,d} = ratemapt{:,2};
                    end
                end
            end
        else
            neuronIndivLRm = cell(length(neuronIndivLR),length(neuronIndivLR));
            behavIndivLRm = cell(length(behavIndivLR),length(behavIndivLR));
            frLRm = cell(length(behavIndivLR),length(behavIndivLR));
            countLRm = cell(length(behavIndivLR),length(behavIndivLR));
            timeLRm = cell(length(behavIndivLR),length(behavIndivLR));
            for ii = 1:length(neuronIndivLR);
                for jj = ii:length(neuronIndivLR);
                [neuronL1, behavL1, neuronL2, behavL2] = downsample_match_position(neuronIndivLR{1,ii},behavIndivLR{1,ii},neuronIndivLR{1,jj},behavIndivLR{1,jj});
                neuronIndivLRm{ii,jj} = neuronL1; behavIndivLRm{ii,jj} = behavL1;
                neuronIndivLRm{jj,ii} = neuronL2; behavIndivLRm{jj,ii} = behavL2;
                end
            end

            for ii = 1:size(neuronIndivLRm,1)
                for jj = 1:size(neuronIndivLRm,2)
                [frLRm{ii,jj},countLRm{ii,jj},timeLRm{ii,jj}] = calculate_subset_ratemap(neuronIndivLRm{ii,jj},behavIndivLRm{ii,jj},thresh,2);
                [frLRm{ii,jj},countLRm{ii,jj},timeLRm{ii,jj}] = calculate_subset_ratemap(neuronIndivLRm{ii,jj},behavIndivLRm{ii,jj},thresh,2);
                end
            end
            % smooth and trim ratemap
            ratemapt_LRm = cell(length(neuronIndivLRm),length(neuronIndivLRm));
            for ii = 1:length(neuronIndivLRm);
                for jj = 1:length(neuronIndivLRm);
                    frMatch = {frLRm{ii,jj},frLRm{jj,ii}};
                    [~, ratemapt] = smooth_trim_ratemap(frMatch, 'top right');
                    ratemapt_LRm{ii,jj} = ratemapt{:,1};
                    ratemapt_LRm{jj,ii} = ratemapt{:,2};
                end
            end
        end

end
end