function betaFIR_PCA = CPCA_SCORES_FIR_REGRESS(X,scores)

% INPUTS:
% X: FIR design matrix, (nTRsperscan*nobs by nobs*ntimepointsmodeled) 
%with NaNs for subjects who don't have a particular response, or who have 
% too many non responses
% scores: temporal loadings of PCA on some part of task related variance
% carried out by FIR_PCA(). (nTRsperscan*nobs by ncomponents)
%
% OUTPUTS:
% betaFIR_PCA: (nobs*ntimepointsmodeled)*ncomponents) matrix of betas
% that measure time course of PCA scores wrt stimuli specified in FIR

% which subjects have missing response data either in regression1 or regression2 - exclude them
ncomponents = size(scores,2);
% remove regressors for subjects EITHER with no response data at all OR
% no responses of a particular type
nanmask_y = ~isnan(sum(X,1));% & ~isnan(sum(X_reg1,1));    
% even after that, can still have some rank deficiency where two isolated
% responses overlap with no additional data points to resolve. remove those.
dup_mask_all = FIR_REGRESSOR_RM_DUPS(X);
nanmask_y(dup_mask_all) = false;
% remove subjects whose data is not captured by any regressors (accounting
% for presence of intercept which gives everybody 1 regressor
% also possible that some regressors are excluded from larger
% concatenated matrix due to duplicates across response types, so exclude
% those responses here as well (judging off of scores which is obtained
% from larger regressor in situation of including all responses in PCA)
nanmask_x = sum(X(:,nanmask_y),2)>1 & ~isnan(sum(scores,2));
% in case this made some regressors all 0's, exclude those too
nanmask_y(sum(X(nanmask_x,:),1)==0) = false;
betaFIR_PCA = nan(size(X,2),ncomponents); % only fit model on subjects with data for responses, but store in same size matrix to preserve subject order
betaFIR_PCA(nanmask_y,:) = X(nanmask_x,nanmask_y)\scores(nanmask_x,:); % subject level betas for principal component loadings of ROI-level betas
