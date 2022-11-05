function [ratemap, ratemap_trim] = smooth_trim_ratemap(firingrateAll, trim_direct, countTime, edgeFilter, edgeFilterThresh)
% this function aims to trim nans outside the actuall ratemap and apply
% smoothing to the raw ratemap firingrateAll. In addition, this function
% also further trim the unequal size ratemap across different days to equal
% size to facilitate further analysis, e.g. 2D correlation
% developed by Yanjun Sun 10/22/2018, 
% updated 1/31/2019 for nature neuroscience paper revision
if ~exist('edgeFilterThresh','var') || isempty(edgeFilterThresh)
    edgeFilterThresh = 0.2; %sec
end
if ~exist('edgeFilter','var') || isempty(edgeFilter)
    edgeFilter = false;
end
if ~exist('trim_direct','var') || isempty(trim_direct)
    trim_direct = 'top left';
end

% to trim the ratemap matrix and smooth them
ratemap = cell(1,length(firingrateAll));
ratemap0 = cell(1,length(firingrateAll));
for jj = 1:length(firingrateAll)
    singleRatemap = firingrateAll{jj};
    for ii = 1:length(singleRatemap)
        singleRatemap{ii} = singleRatemap{ii}(any(~isnan(singleRatemap{ii}),2),:); % remove nan - rows
        singleRatemap{ii} = singleRatemap{ii}(:,any(~isnan(singleRatemap{ii})));   % remove nan - columns
        singleRatemap{ii} = filter2DMatrices(singleRatemap{ii}, 1);
    end
    ratemap{jj} = singleRatemap;
    if edgeFilter
        countTime{jj} = countTime{jj}(any((countTime{jj}),2),:); % remove zeros - rows
        countTime{jj} = countTime{jj}(:,any((countTime{jj}))); % remove zeros - columns
    end
end
% save Ratemap.mat ratemap; %ratemap can be used for plotting

% to replace nans to 0 in the ratemap
for jj = 1:length(firingrateAll);
    singleratemap = ratemap{jj};
    for ii = 1:length(singleratemap);
        singleratemap{ii}(isnan(singleratemap{ii})) = 0;
    end
    ratemap0{jj} = singleratemap;
end
% to filter out edge cols with sporadic bins that contains very low
% occupancy time
if edgeFilter
    for ii = 1:length(countTime);
        if sum(countTime{1,ii}(:,1)) < edgeFilterThresh %set sporadic bins on the edge to zero
            countTime{1,ii}(:,1) = 0;
        end
        if sum(countTime{1,ii}(:,end)) < edgeFilterThresh
            countTime{1,ii}(:,end) = 0;
        end
        cols = find(~all(countTime{1,ii}==0,1));
        rows = find(~all(countTime{1,ii}==0,2));
        for jj = 1:length(ratemap0{1,ii});
            singlecountAll = ratemap0{1,ii}{1,jj};
            if isempty(singlecountAll);
                singlecountAll = zeros(length(rows), length(cols));
            else
                singlecountAll = singlecountAll(:,cols); % remove nan - rows
                singlecountAll = singlecountAll(rows,:);   % remove nan - columns
            end
            ratemap_trim{1,ii}{1,jj} = singlecountAll;
        end
    end
else
    ratemap_trim = ratemap0;
end

% to further trim circle ratemap to a universal size across different days
for ii = 1:length(ratemap_trim);
    for jj = 1:length(ratemap_trim{ii});
    [row,col] = size(ratemap_trim{ii}{jj});
    if row >0 && col >0;
        break
    end
    end
    ratemapsize_circle(ii,1) = row;
    ratemapsize_circle(ii,2) = col;
end
sizelimit_circle = min(ratemapsize_circle);

for ii = 1:length(ratemap_trim);
    singleratemap = ratemap_trim{ii};
        for jj = 1:length(singleratemap);
            if isempty(singleratemap{jj});
                singleratemap{jj} = zeros(sizelimit_circle(1,1),sizelimit_circle(1,2));
            else
                %current trim takes out additional rows/cols from left and top of the matrix. change the trim direction if needed
                switch trim_direct;
                    case {'top left','Top Left','Top left'};
                        singleratemap{jj} = singleratemap{jj}(size(singleratemap{jj},1)-sizelimit_circle(1,1)+1:end, size(singleratemap{jj},2)-sizelimit_circle(1,2)+1:end);
                    case {'top right','Top Right','Top right'};
                        singleratemap{jj} = singleratemap{jj}(size(singleratemap{jj},1)-sizelimit_circle(1,1)+1:end, 1:sizelimit_circle(1,2));  
                    case {'bottom left','Bottom Left','Bottom left'};
                        singleratemap{jj} = singleratemap{jj}(1:sizelimit_circle(1,1), size(singleratemap{jj},2)-sizelimit_circle(1,2)+1:end);  
                    case {'bottom right','Bottom Right','Bottom right'};
                        singleratemap{jj} = singleratemap{jj}(1:sizelimit_circle(1,1), 1:sizelimit_circle(1,2));  
                end
            end
        end
    ratemap_trim{ii} = singleratemap;
end
% save Ratemap_trim.mat ratemap_trim;
end
