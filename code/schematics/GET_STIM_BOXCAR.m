function [boxcar] = GET_STIM_BOXCAR(stim_onset,stim_duration,TR,T)

	% INPUTS:
	% stim_onset: vector indicating stimulus onset times in seconds post scan start
	% stim_duration: vector indicating duration of each stimulus
	% TR: repetition time (1/sampling rate) in seconds
	% T: number of time points in scan (convolution will add extra time points to finish out kernel)
	%
	% OUTPUTS:
	% boxcar: stim convolved with hemodynamic response function, truncated at T time points
    
	% this function is a hacked version of CONV_STIM_BOXCAR which uses spm hrf
	% this function takes out the spm part and just generates a boxcar version of the stimulus

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
    boxcar = boxcar_upsample(floor(linspace(upsampleFactor/2,T-upsampleFactor/2,T/upsampleFactor)));
