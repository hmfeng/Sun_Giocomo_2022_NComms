function [Corr_condition, Corr_condition2] = get_corr_matrix(PlaceCells_batch, ratemap_condition, ratemap_condition2, condition);
% This function aims to get ratemaps 2D correlation across different days
% for each place cell, which defined by the particular day and condition.
% input argument
% PlaceCells_batch: shuffled results of place cells for all the sessions
% ratemap_condition: trimed ratemap for each condition e.g. ratemap_square
% condition: 'Circle', 'Square', 'Any', 'LR'for now
% output:
% Corr_condition: a cell array contains all the correlation matrices for
% each neuron
% Corr_condition2: transform Corr_condition into a 3d array, but note for
% condition LR is a little different
% developed by Yanjun Sun, 10/23/2018
% updated by Yanjun Sun, 3/21/2019, with 'LR' condition

if ~exist('condition','var') || isempty(condition)
    condition = 'any';
end
% it will calculate all the cells if there is no place cell number feed in
if ~exist('PlaceCells_batch','var') || isempty(PlaceCells_batch)
    PlaceCells_batch = {};
    PlaceCells_batch{1,1} = (1:length(ratemap_condition{1,1}))';
    PlaceCells_batch = repmat(PlaceCells_batch,1,length(ratemap_condition)); 
end
% convert PlaceCells_batch to a cell array if input is not
if ~iscell(PlaceCells_batch)
    PlaceCells_batch2 = {};
    PlaceCells_batch2{1,1} = PlaceCells_batch;
    PlaceCells_batch = repmat(PlaceCells_batch2,1,length(ratemap_condition)); 
end

%% start to calculate correlations
switch condition;
case {'Circle','circle'};
Corr_condition = {};
Avg_Corr_condition = {};
ct = 0; %count loop time
for session_select = 1:2:length(PlaceCells_batch);
    ct = ct+1;
for k = 1:length(PlaceCells_batch{1,session_select});
    kk = PlaceCells_batch{1,session_select}(k,1);
    for ii = 1:length(ratemap_condition);
        for jj = 1:length(ratemap_condition);
            Fcorr(ii,jj) = corr2(ratemap_condition{ii}{kk},ratemap_condition{jj}{kk});
        end
    end
    Corr_condition{ct,k} = Fcorr;
end
    Corr_condition2 = cat(3,Corr_condition{ct,:});
    Avg_Corr_condition{ct} = nanmean(Corr_condition2,3);
    % plot averaged correlation matrix
    figure;
    imagesc(Avg_Corr_condition{ct});
    axis square;
end

case {'Square','square'};
Corr_condition = {};
Avg_Corr_condition = {};
ct = 0; %count loop time
for session_select = 2:2:length(PlaceCells_batch);
    ct = ct+1;
for k = 1:length(PlaceCells_batch{1,session_select});
    kk = PlaceCells_batch{1,session_select}(k,1);
    for ii = 1:length(ratemap_condition);
        for jj = 1:length(ratemap_condition);
            Fcorr(ii,jj) = corr2(ratemap_condition{ii}{kk},ratemap_condition{jj}{kk});
        end
    end
    Corr_condition{ct,k} = Fcorr;
end
    Corr_condition2 = cat(3,Corr_condition{ct,:});
    Avg_Corr_condition{ct} = nanmean(Corr_condition2,3);
    % plot averaged correlation matrix
    figure;
    imagesc(Avg_Corr_condition{ct});
    axis square;
end

case 'any';
Corr_condition = {};
Avg_Corr_condition = {};
ct = 0; %count loop time
for session_select = 1:length(PlaceCells_batch);
    ct = ct+1;
for k = 1:length(PlaceCells_batch{1,session_select});
    kk = PlaceCells_batch{1,session_select}(k,1);
    for ii = 1:length(ratemap_condition);
        for jj = 1:length(ratemap_condition);
            Fcorr(ii,jj) = corr2(ratemap_condition{ii}{kk},ratemap_condition{jj}{kk});
        end
    end
    Corr_condition{ct,k} = Fcorr;
end
    Corr_condition2 = cat(3,Corr_condition{ct,:});
    Avg_Corr_condition{ct} = nanmean(Corr_condition2,3);
    % plot averaged correlation matrix
    figure;
    imagesc(Avg_Corr_condition{ct});
    axis square;
end

case {'LR', 'L&R','Left&Right', 'left&right', 'left right', 'Left Right'}
Corr_condition = {};
Avg_Corr_condition = [];
for session_select = 1:length(PlaceCells_batch);
    for ii = 1:length(ratemap_condition{1,1});
        if and(sum(sum(ratemap_condition{1,session_select}{1,ii}))==0, sum(sum(ratemap_condition2{1,session_select}{1,ii}))==0);
            Fcorr(1,ii) = 1;
        elseif or(sum(sum(ratemap_condition{1,session_select}{1,ii}))==0, sum(sum(ratemap_condition2{1,session_select}{1,ii}))==0);
            Fcorr(1,ii) = 0;
        else
            Fcorr(1,ii) = corr2(ratemap_condition{1,session_select}{1,ii},ratemap_condition2{1,session_select}{1,ii});
        end
    end
    Corr_condition{1,session_select}=Fcorr';
end
    figure;
    Corr_condition2 = cell2mat(Corr_condition);
    imagesc(Corr_condition2);

case {'match', 'Match', 'matched', 'Matched'} %no place cell feed in
Corr_condition = cell(length(ratemap_condition), length(ratemap_condition));
for ii = 1:length(ratemap_condition)
    for jj = 1:length(ratemap_condition)
        for n = 1:length(ratemap_condition{ii,jj});
            if and(sum(sum(ratemap_condition{ii,jj}{1,n}))==0, sum(sum(ratemap_condition{jj,ii}{1,n}))==0);
                Fcorr(1,n) = 1;
            elseif or(sum(sum(ratemap_condition{ii,jj}{1,n}))==0, sum(sum(ratemap_condition{jj,ii}{1,n}))==0);
                Fcorr(1,n) = 0;
            else
                Fcorr(1,n) = corr2(ratemap_condition{ii,jj}{1,n},ratemap_condition{jj,ii}{1,n});
            end
        end
        Corr_condition{ii,jj} = Fcorr';
        Corr_condition{jj,ii} = Fcorr';
    end
end
Corr_condition2 = cellfun(@mean, Corr_condition);

    case{'xsessions'}
    %Compute ratemap correlations across different sessions for each
    %neuron. Using cellfun instead of loop. Oct.2021
    N = length(ratemap_condition);
    Corr_condition = cell(N,N);
    for ii = 1:N
        for jj = 1:N
            A = ratemap_condition{ii};
            B = ratemap_condition{jj};
            %replace NaN values to zero
            idxA = cellfun(@(x) isnan(x),A,'uni',0);
            idxB = cellfun(@(x) isnan(x),B,'uni',0);
            for n = 1:length(A)
                A{n}(idxA{n})=0;
                B{n}(idxB{n})=0;
            end
            %compute correlation for each cell
            Corr_condition{ii,jj} = cellfun(@(x,y) corr2(x,y), A,B);
        end
    end
    Corr_condition2 = cellfun(@mean, Corr_condition);

end
end