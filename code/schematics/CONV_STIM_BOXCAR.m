function [conv_stim,boxcar,boxcar_upsample] = CONV_STIM_BOXCAR(stim_onset,stim_duration,TR,T)

	% INPUTS:
	% stim_onset: vector indicating stimulus onset times in seconds post scan start
	% stim_duration: vector indicating duration of each stimulus
	% TR: repetition time (1/sampling rate) in seconds
	% T: number of time points in scan (convolution will add extra time points to finish out kernel)
	%
	% OUTPUTS:
	% conv_stim: stim convolved with hemodynamic response function, truncated at T time points
    
    if length(stim_duration) == 1
        stim_duration = repmat(stim_duration,[length(stim_onset) 1]);
    end
	upsampleFactor = 100;
	T = T*upsampleFactor;
	TR = TR/upsampleFactor;	
	boxcar_upsample = zeros(T, 1); % initialize stimulus vector		
	stim_indices = floor(stim_onset/TR + 1); % shift so that you index the start of each TR -- i.e. if stimulus onsets at 3s, it's starting during TR 2
	stim_duration = ceil(stim_duration/TR); % number of TRs for each stimulus
    for iStim = 1:length(stim_indices)
		boxcar_upsample(stim_indices(iStim):stim_indices(iStim) + stim_duration(iStim)) = 1; % construct upsampled boxcar vector of stimulus
	end
    %% convolve binary regressor with HRF
    hrf = spm_hrf(TR); %get the HRF in the correct time resolution from SPM
    conv_stim = conv(boxcar_upsample, hrf); %stim is a binary vector indicating which TRs contained a stimuli during the scan
    conv_stim = conv_stim(1:T);
    % downsample convolved stimulus by taking a sample in the center of each TR
    conv_stim = conv_stim(floor(linspace(upsampleFactor/2,T-upsampleFactor/2,T/upsampleFactor)));
    boxcar = boxcar_upsample(floor(linspace(upsampleFactor/2,T-upsampleFactor/2,T/upsampleFactor)));
