function Xout = HORZ_WHITEPAD(Xin,sz)
	% INPUTS:
	% Xin: image
	% sz: target size
	%
	% OUTPUTS:
	% Xout: Xin with symmetric white padding on left and right to make it size sz 
	pad = 0.5*(sz - size(Xin,2));
	Xout = [255*ones(size(Xin,1),floor(pad),3),Xin,255*ones(size(Xin,1),ceil(pad),3)];