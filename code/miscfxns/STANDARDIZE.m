function X_std = STANDARDIZE(X,zdim)
	% demean or z-score time-by-region time series X depending on indicator zdim
	% if zdim = 0, just demean each region
	% if zdim = 1 or 2, then z-score each region (column) or each time point (row), respectively
	% if zdim = 3, the n z-zscore each region (column)
	if zdim == 0
		X_std = X - mean(X,1); % demean each region's time series
	elseif zdim < 3
		X_std = zscore(X,[],zdim); % z-score each region's time series (zdim=1) or each time point (zdim=2)
	elseif zdim == 3
		X_std = zscore(X,[],1); % z-score each region's time series (zdim=1)
	end