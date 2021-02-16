% make schematic demonstrating FIR PCA method
addpaths;
cd(basedir);
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile('data',['TimeSeriesIndicators',name_root,'.mat']));
load(fullfile('data',['ConcTimeSeries',name_root,'.mat']));
masterdir = fullfile('results',name_root);
savedir_base = fullfile(masterdir,'schematics','fir_pca');
savedir_fir = fullfile(masterdir,'analyses','fir');
mkdir(savedir_base);
concTS = THRESHOLD(concTS,zdim);
%% fit FIR + PCA model to generate TS, TS_task and TS_rest

N=100; K = 8; TR =3; T = 50; fin=6;
nTRs = allScanTRs(1);


allstim = load(fullfile('data/task/idemo/stimulus/all_3col.txt'));
[~,boxcar] = CONV_STIM_BOXCAR(allstim,5.5,TR,nTRs);
%regressor = [ones(nTRs,1) GET_FIR_REGRESSOR(allstim(:,1),6,TR,nTRs)];

% load subject regressor
FIR_Design_Reg1 = load(fullfile(masterdir,'analyses','fir','design_matrices',['ThreatNonthreatAllStimuliStratified_FIRDesignMatrix_fin',num2str(fin),'st',num2str(st),'.mat']));
regressor = [ones(nTRs,1) FIR_Design_Reg1.X(subjInd == N,FIR_Design_Reg1.columnLabelsSubject==N)];
TS = concTS(subjInd==N,:)
%FIR_Design_Reg1 = load(fullfile(savedir_fir,'design_matrices',['AllStimulus_FIRDesignMatrix.mat']),'X');
%[coeff,scores,explained,TS_reconstruct] = FIR_PCA(FIR_Design_Reg1.X,concTS); FIRDesignMatrix_fin',num2str(fin),'st',num2str(st),'.mat'% use scores from input components

% fit FIR and PCA on single subject
[coeff,scores,explained,X_task] = FIR_PCA(regressor,TS); % use scores from input components

% select K regions with the most task variance for visualization
[~,idx] = sort(var(X_task),'descend');
P = idx(1:K);
%P = randperm(nparc,K);

% plot FC
f=figure;
subplot(1,3,1); imagesc(corr(TS - X_task,'rows','pairwise')); caxis([-1 1]); axis square
subplot(1,3,2); imagesc(corr(X_task,'rows','pairwise')); caxis([-1 1]); axis square
subplot(1,3,3); imagesc(corr(TS)); caxis([-1 1]); axis square

X = TS(1:T,P); X_task = X_task(1:T,P);
X_rest = X - X_task;
X_task = X_task.*repmat(0.5*(max(X)./max(X_task)),[T 1]);
boxcar = boxcar(1:T);
% only plot first T time points

%% plot BOLD time series = task + rest

offset = repmat(3 + 3*[1:K],[T 1]);
t_tr = linspace(0,T,4); % time in trs

MPL = load('data/colors/mpl_cmaps.mat')'
cmap = MPL.Blues; stim_col = [0.2 0.2 0.2]; lw = 2;
cmap = cmap(round(linspace(1,size(MPL.Blues,1),K)),:); % downsample color map to use discretely and evenly spaced
scl = 0.6; % scale all the time series if overlapping

% plot full time series
f=figure;
subplot(1,3,1);
hL = plot(scl*X + offset,'LineWidth',lw); hold on;
for i=1:length(hL),set(hL(i),'color',cmap(i,:));end

plot(boxcar+0.5,'Color',stim_col);
%yticks([1, 3 + 3*[1:K]]); 
%yticklabels([{'Stimulus'},cellfun(@(x) sprintf('Region %d',x),num2cell(1:K),'UniformOutput',false)]);
yticks([]); 
ylabel('Region');
xlim([0 T]); xticks(t_tr); xticklabels(TR*t_tr); xlabel('Time (s)');
title('BOLD Time Series');
prettifyEJC;
yl = get(gca,'YLim');
ylim(yl);
% plot task-related and task-unrelated variance separately


subplot(1,3,2);
hL = plot(scl*X_task + offset,'LineWidth',lw); hold on;
for i=1:length(hL),set(hL(i),'color',cmap(i,:));end
plot(boxcar+0.5,'Color',stim_col);
%yticks([1, 3 + 3*[1:K]]); 
%yticklabels([{'Stimulus'},cellfun(@(x) sprintf('Region %d',x),num2cell(1:K),'UniformOutput',false)]);
yticks([]);
xlim([0 T]); xticks(t_tr); xticklabels(TR*t_tr); xlabel('Time (s)');
title('Task-Related Signals');
prettifyEJC;
ylim(yl);

%{
% plot task-unrelated signals -- looks same as overall signal
subplot(1,3,3);
cmap = repmat([0.7 0.7 0.7],[K 1]);
for j = 1:K
    plot(X_rest(:,j) + offset(:,j),'Color',cmap(j,:)); hold on;
end
plot(boxcar+0.5,'Color',stim_col);
yticks([1, 3 + 3*[1:K]]); 
yticklabels([{'Stimulus'},cellfun(@(x) sprintf('Region %d',x),num2cell(1:K),'UniformOutput',false)]);
xlim([0 T]); xticks(t_tr); xticklabels(TR*t_tr); xlabel('Time (s)');
title('Task-Unrelated Signals');
prettifyEJC;
ylim(yl);
%}

%
K_pc = 3;
scores_plot = scores(1:T,1:K_pc);
scores_plot = scores_plot.*repmat((0.5*max(X(:,1:K_pc))./max(scores_plot)),[T 1]); 
offset = repmat(3 + 3*(7/8)*(K/K_pc)*[1:K_pc],[T 1]);
subplot(1,3,3);

hL=plot(scl*scores_plot + offset(:,fliplr(1:K_pc)),'LineWidth',lw); hold on;
for i=1:length(hL),set(hL(i),'color',cmap(i+2,:));end
plot(boxcar+0.5,'Color',stim_col);
%yticks([1, 3 + 3*(7/8)*(K/K_pc)*[1:K_pc]]); 
%yticklabels([{'Stimulus'},cellfun(@(x) sprintf('PC %d',x),fliplr(num2cell(1:K_pc)),'UniformOutput',false)]);
yticks([]);
ylabel('Principal Components');
xlim([0 T]); xticks(t_tr); xticklabels(TR*t_tr); xlabel('Time (s)');
title('Modes of Task Activity');
prettifyEJC;
ylim(yl)
%}

f=FIGURE_SIZE_CM(f,22,6);
saveas(f,fullfile(savedir_base,'FIRPCA_Schematic.pdf'));
%% plot Task BOLD time series = [6x column vector of brains fading into orthogonal time series]