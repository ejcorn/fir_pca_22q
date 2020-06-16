function [allScanTRs,scanlab,BOLDpaths] = SET_SCANS(scan,idemoNumTRs,idemofnames,restNumTRs,restfnames)

if strcmp(scan,'C')
	allScanTRs = [idemoNumTRs,restNumTRs];
	scanlab = {'IDEmo','Rest'};
	BOLDpaths = {idemofnames,restfnames};
elseif strcmp(scan,'ID')
	allScanTRs = [idemoNumTRs];
	scanlab = {'IDEmo'};
	BOLDpaths = {idemofnames};
elseif strcmp(scan,'R')
	allScanTRs = [restNumTRs];
	scanlab = {'Rest'};
	BOLDpaths = {restfnames};
end