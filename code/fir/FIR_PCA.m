function [coeff,scores,explained,TS_reconstruct,betaFIR] = FIR_PCA(X,concTS)
    % INPUTS:
    % X: FIR design matrix, (nsubjs*nTRperscan)x(nsubjs*n_timepointsmodeled). plus the first
    % column is ones for intercept
    % concTS: (nsubjs*nTRperscan)xnregions time series
    %
    % OUTPUTS:
    % coeff,scores, explained from PCA on variance explained by X in concTS
    % if response data is missing leading to NaNs in X, these are ignored
    % in regression and PA
    % coeff: nregions x nregions
    % scores: (nsubjs*nTRperscan)x nregions    
    % TS_reconstruct: (nsubjs*nTRperscan)x nregions, variance in concTS captured by X
    % betaFIR: (nsubjs*ntimepointsmodeled)x nregions, betas from
    % regressing concTS on X
        
    % remove regressors for subjects EITHER with no response data at all OR
    % no responses of a particular type
    nanmask_y = ~isnan(sum(X,1)); 
    % remove subjects whose data is not captured by any regressors
    % (accounting for intercept that makes everybody sum to at least 1)
    nanmask_x = sum(X(:,nanmask_y),2)>1; 
    
    betaFIR = nan(size(X,2),size(concTS,2));
    tic
    betaFIR(nanmask_y,:) = [X(nanmask_x,nanmask_y)\concTS(nanmask_x,:)];
    toc
    TS_reconstruct = nan(size(concTS)); % reconstruct time series containing only variance explained by task
    TS_reconstruct(nanmask_x,:) = X(nanmask_x,nanmask_y)*betaFIR(nanmask_y,:);
    tic
    [coeff,scores,~,~,explained] = pca(TS_reconstruct);
    toc