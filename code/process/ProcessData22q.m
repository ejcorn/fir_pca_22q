%basedir = '/cbica/home/cornblae/ecornblath/fir_pca_22q/'; scan = 'ID'; zdim = 1; atlasName = 'HarvardOxford'; atlasScale = 112; extralabel = ''; XCP_folder = 'xcp_36p_despike';
cd(basedir);

name_root = ['CPCA_',scan,atlasName,num2str(atlasScale),'Z',num2str(zdim),XCP_folder,extralabel];
datadir = fullfile(basedir,'data');
dtiPathRoot = fullfile(datadir,'qsirecon/'); % this needs to be specified in advance
masterdir = fullfile(basedir,'results',name_root);
mkdir(masterdir);
addpath(genpath(fullfile(basedir,'code')))

%% Load sample of age-sex matched 22q and PNC subjects
grouplab = {'PNC','22q'};
demo_fname = fullfile(datadir,'PNC_PreQC_sample.csv');
demoMatch = readtable(demo_fname,'TreatAsEmpty','NA');
demoMatch = demoMatch(~isnan(demoMatch.idemo_meanrelrms),:);    % exclude people with missing motion for idemo... corresponds with missing files
idemoExcludePairs = demoMatch.match_pair(demoMatch.idemo_exclude == 1);     % apply idemo exclusion
demoMatch = demoMatch(~ismember(demoMatch.match_pair,idemoExcludePairs),:); % remove pairs for which idemo data is poor quality in 22q subject
nobs = length(demoMatch.bblid);
nobs_22q = sum(strcmp(demoMatch.study,'22q'));
nobs_PNC = nobs - nobs_22q;

disp('demo loaded');

%% Declare BOLD file names and choose scan composition

if strcmp(atlasName,'Laus')
	if atlasScale == 250
		nparc = 462;
	elseif atlasScale == 125
		nparc = 233;
	end
	prefix = [lower(atlasName),num2str(atlasScale)];
elseif strcmp(atlasName,'Brainnetome')
	nparc = 246;
	prefix = [lower(atlasName),num2str(atlasScale)];
elseif strcmp(atlasName,'Schaefer')
	if atlasScale == 200
		nparc = 200;
	elseif atlasScale == 400
		nparc = 400;
	end
	XCP_atlasName = [lower(atlasName),num2str(atlasScale),'x7']; % specify 7 networks
elseif strcmp(atlasName,'HarvardOxford')    
	nparc = 112;
	XCP_atlasName = atlasName;
end
%
atlasScale_str = num2str(atlasScale); % for use in filepaths

% use xcp data that comes parcellated
% ~~PATH~~ beginning stem for XCP data from Arun .. middle asterisk is separate 22q and pnc folders which doesn't matter bc scanids and bblids don't overlap    
idemoFileRoot_PNC = ['/cbica/home/mahadeva/22q_project/pnc_sample2/derivatives_uniqueSession/',XCP_folder,'/sub-']; 
idemoFileRoot_22q = ['/cbica/home/mahadeva/22q_project/22q_deleted/derivatives/',XCP_folder,'/sub-'];
bblids_22q = num2str(demoMatch.bblid(strcmp(demoMatch.study,'22q'))); % no 0 padding
bblids_PNC = cellstr(string(demoMatch.bblid(strcmp(demoMatch.study,'pnc_sample'))));% no 0 padding -- pnc has different length ids
scanIDs_22q = num2str(demoMatch.scanid(strcmp(demoMatch.study,'22q'))); % no 0 padding
scanIDs_PNC = cellstr(string(demoMatch.scanid(strcmp(demoMatch.study,'pnc_sample')))); % no 0 padding -- pnc has different length ids
session_PNC = cellstr(string(demoMatch.pnc_session(strcmp(demoMatch.study,'pnc_sample')))); % get session number for the scan for each pnc subject that has been matched to the single scan from a 22q subject

idemofnames_22q = cellfun(@(bblid,scanid) [idemoFileRoot_22q,bblid,'/*',scanid,'/idemo2/fcon/',XCP_atlasName,'/sub-',bblid,'*',scanid,'_idemo2_',XCP_atlasName,'_ts.1D'],cellstr(bblids_22q),cellstr(scanIDs_22q),'UniformOutput',false);
%idemofnames_22q = cellstr([repmat(idemoFileRoot_22q,[nobs_22q 1]),bblids_22q,repmat('/*',[nobs_22q 1]),scanIDs_22q,repmat(['/idemo2/fcon/schaefer',atlasScale_str,'x7/sub-'],[nobs_22q 1]),bblids_22q,repmat('*',[nobs_22q 1]),scanIDs_22q,repmat(['_idemo2_schaefer',atlasScale_str,'x7_ts.1D'],[nobs_22q 1])]);
% for the only 22q subject with multiple sessions, use this separate xcp outfolder
idemofnames_22q(find(strcmp(bblids_22q,'16095'))) = {['/cbica/home/mahadeva/22q_project/22q_deleted/derivatives_sub-16095/',XCP_folder,'/sub-16095/ses-20140715x8971/idemo2/fcon/',XCP_atlasName,'/sub-16095*_idemo_',XCP_atlasName,'_ts.1D']};
cellfun(@disp,idemofnames_22q,'UniformOutput',false)
idemofnames_22q = cellfun(@dir,idemofnames_22q,'UniformOutput',false); % i don't know all session numbers so just infer them assuming 1 session from each subject with scanid and bblid
idemofnames_22q = cellfun(@(x) fullfile(x(1).folder,x(1).name),idemofnames_22q,'UniformOutput',false); % take first session with (end)

idemofnames_PNC = cellfun(@(bblid,ses) [idemoFileRoot_PNC,bblid,'PNC',ses,'/ses-PNC',ses,'/idemo/fcon/',XCP_atlasName,'/sub-',bblid,'*_idemo_',XCP_atlasName,'_ts.1D'],bblids_PNC,session_PNC,'UniformOutput',false);
% for now just include subject 96030 because they have only one session
idemofnames_PNC(find(strcmp(bblids_PNC,'96030'))) = {['/cbica/home/mahadeva/22q_project/pnc_sample2/derivatives/',XCP_folder,'/sub-96030/ses-PNC1/idemo/fcon/',XCP_atlasName,'/sub-96030*_idemo_',XCP_atlasName,'_ts.1D']};
idemofnames_PNC = cellfun(@dir,idemofnames_PNC,'UniformOutput',false); % i don't know all session numbers so just infer them assuming 1 session from each subject with scanid and bblid
idemofnames_PNC = cellfun(@(x) fullfile(x(end).folder,x(end).name),idemofnames_PNC,'UniformOutput',false); % take most recent session with (end) .. if PNc has multiple sessions ... need
%}

% concatenate group label and scan id next to path

idemofnames = [idemofnames_22q repmat({'22q'},[nobs_22q 1]) regexprep(cellstr(scanIDs_22q),'^0*',''); idemofnames_PNC repmat({'pnc_sample'},[nobs_PNC 1]) regexprep(cellstr(scanIDs_PNC),'^0*','')];

restFileRoot = [datadir,'/',atlasName,num2str(atlasScale),'rest/rest_',atlasName,num2str(atlasScale),'_']; % beginning stem to all parcellated rest BOLD files
restfnames_22q = cellstr([repmat(restFileRoot,[nobs_22q 1]),num2str(demoMatch.bblid(strcmp(demoMatch.study,'22q'))),repmat('.txt',[nobs_22q 1])]);
restfnames_PNC = cellstr([repmat(restFileRoot,[nobs_PNC 1]),num2str(demoMatch.bblid(strcmp(demoMatch.study,'pnc_sample'))),repmat('.txt',[nobs_PNC 1])]);
% concatenate group label next to path
restfnames = [restfnames_22q repmat({'22q'},[nobs_22q 1]) regexprep(cellstr(scanIDs_22q),'^0*',''); restfnames_PNC repmat({'pnc_sample'},[nobs_PNC 1]) regexprep(cellstr(scanIDs_PNC),'^0*','')];

restNumTRs = 120; idemoNumTRs = 204;
% either analyze rest scans only, idemo scans only, or both together
[allScanTRs,scanlab,BOLDpaths] = SET_SCANS(scan,idemoNumTRs,idemofnames,restNumTRs,restfnames);

%% Concatenate BOLD data
totalNumTRs = sum(allScanTRs);
concTS = nan(totalNumTRs*nobs,nparc);
subjInd = nan(nobs*totalNumTRs,1); % indexes subjects
scanInd = nan(nobs*totalNumTRs,1); % indexes scans
SubjectIs22q = nan(nobs*totalNumTRs,1); % indexes which TRs are 22q
subjInd_scanID = cell(length(BOLDpaths),1); % holds scan ids for each subjInd entry. subjInd 1:nobs will be iterated through in other scripts to compute metrics. so in order to match with demographics need scanid
for j = 1:length(BOLDpaths)
	subjInd_scanID{j} = cell(nobs,1); % loop through subjects and store scanid in same loop as bold concatenation to be sure data is linked
end
for N = 1:nobs
	disp(['Subject ', num2str(N)]); 
	% iterate through subjects, then through scans, and add BOLD data to concTS matrix
	for S = 1:length(BOLDpaths) 
		disp(scanlab{S})
		st = 1 + allScanTRs(S)*(N-1); 
		nd = allScanTRs(S)*N;
		if S > 1
			st = 1 + sum(allScanTRs(1:(S-1)))*(nobs) + allScanTRs(S)*(N-1);    % all previous scans offset the index for scan i
			nd = sum(allScanTRs(1:(S-1)))*nobs + allScanTRs(S)*(N);
		end
		tmp = dlmread(BOLDpaths{S}{N,1});         
		tmp = tmp(7:end,1:nparc); % remove brainstem and first 6 volumes
		concTS(st:nd,:) = STANDARDIZE(tmp,zdim); % demean or zscore data for each subject, depending on zdim (see code/miscfxns/STANDARDIZE.m)
		scanInd(st:nd) = S;
		subjInd(st:nd) = N;
		% for a given subject label whether they are 22q
		SubjectIs22q(st:nd) = strcmp(BOLDpaths{S}{N,2},'22q');
		subjInd_scanID{S}{N} = BOLDpaths{S}{N,3}; % store scan id
	end
end

% scanids are in the same order as demo but this link is just to be positive
% isequal(subjInd_scanID{1},regexprep(cellstr(num2str(demoMatch.scanid,'%06d')),'^0*',''))

disp('Concatenation done');

%% load structural connectivity estimated from diffusion-weighted imaging
% this data is found in qsirecon folder
% dwi is not present for all subjects and some subjects had multiple runs
SC_IDs_load = [demoMatch.bblid(strcmp(demoMatch.study,'22q'));demoMatch.scanid(strcmp(demoMatch.study,'pnc_sample'))]; % used bblid to label 22q and scanid to label PNC for loading
SC_IDs = [demoMatch.bblid(strcmp(demoMatch.study,'22q'));demoMatch.bblid(strcmp(demoMatch.study,'pnc_sample'))]; % used bblid to label 22q and bblid to label PNC for future analyses
dti_missing = ones(nobs,1); % flag subjects with DTI missing -- default is all missing, mark 0 as you look through files

nparc_sc = 214; % SC only comes in schaefer + HO subcort L then R
% reorder the subcortex so it alternates L node, R node, L node, R node rather than 7 left, 7 right, following data/nifti/HarvardOxford/HarvardOxfordNodeNames.txt
sc_reorder_idx = [1:200,201,208,202,209,203,210,204,211,206,213,207,214,205,212]; 
QA_pass = nan(nparc_sc,nparc_sc,nobs); % matrix to store fractional anisotropy
GFA_Pass = nan(nparc_sc,nparc_sc,nobs); % matrix to store quantitative anisotropy
SIFT_radius2 = nan(nparc_sc,nparc_sc,nobs); % matrix to store prob track weighted streamline count
for N = 1:nobs    
	subj = num2str(SC_IDs_load(N));
	disp(['Subject ',num2str(N),' of ',num2str(nobs),', ID: ',num2str(subj)]);
	subjPath_det = dir(fullfile(dtiPathRoot,['sub-',subj],'dwi',['sub-',subj,'_acq-DTI_*space-T1w_desc-preproc_space-T1w_gqinetwork.mat'])); % deterministic tractography path
	subjPath_prob = dir(fullfile(dtiPathRoot,['sub-',subj],'dwi',['sub-',subj,'_acq-DTI_*space-T1w_desc-preproc_space-T1w_dhollanderconnectome.mat'])); % probabilistic tractography path
	if length(subjPath_det)>0 % everybody has either no dwi or both deterministic and probabilistic tractography
		dti_missing(N) = 0; % add 0 to indicate DTI is NOT missing
		subjQSIRecon_det = load(fullfile(subjPath_det.folder,subjPath_det.name),'SchaeferHOsub214x7_qa_pass_connectivity','SchaeferHOsub214x7_gfa_pass_connectivity');
		QA_pass(:,:,N) = subjQSIRecon_det.SchaeferHOsub214x7_qa_pass_connectivity(sc_reorder_idx,sc_reorder_idx);		
		GFA_Pass(:,:,N) = subjQSIRecon_det.SchaeferHOsub214x7_gfa_pass_connectivity(sc_reorder_idx,sc_reorder_idx);
		subjQSIRecon_prob = load(fullfile(subjPath_prob.folder,subjPath_prob.name),'SchaeferHOsub214x7_sift_radius2_count_connectivity');
		SIFT_radius2(:,:,N) = subjQSIRecon_prob.SchaeferHOsub214x7_sift_radius2_count_connectivity(sc_reorder_idx,sc_reorder_idx);
	end
end

savedir = fullfile(datadir,'sc'); mkdir(savedir); % make directory to hold SC data
sc_data_names = {'QA_pass','GFA_Pass','SIFT_radius2'}; % save each SC metric in a separate file as a variable called SC
sc_data = {QA_pass,GFA_Pass,SIFT_radius2};
for iSC = 1:numel(sc_data_names)
	SC = sc_data{iSC};
	save(fullfile(savedir,['StructuralConnectivity',sc_data_names{iSC},name_root,'.mat']),'SC','dti_missing','SC_IDs');
end

disp('SC loaded');
%}
%% Save data

cd(datadir);
save(['TimeSeriesIndicators',name_root,'.mat'],'scanInd','subjInd','SubjectIs22q','subjInd_scanID');
%save(['VolNormSC',num2str(atlasScale),'.mat'],'SCvolnorm');
save(['Demographics',name_root,'.mat'],'demoMatch','nparc','nobs','nobs_PNC','nobs_22q','atlasScale','atlasName','XCP_folder','scan','zdim','scanlab','grouplab','allScanTRs','basedir','datadir','masterdir');
writetable(demoMatch,['Demographics',name_root,'.csv'],'Delimiter',',')
save(['ConcTimeSeries',name_root,'.mat'],'concTS');
disp('.mat data saved');