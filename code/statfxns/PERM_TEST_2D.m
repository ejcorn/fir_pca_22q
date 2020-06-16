function [pvals_twotail,observedA_gt_B,expectedA_gt_B] = PERM_TEST_2D(A,B,nperms)

    % test null hypothesis that group average transition matrix elements in A do not differ from group average of B
    % by swapping 50% of the membership between A and B and recomputing difference between A and B

    % INPUT:
    % A and B: N_axK / N_bxK matrices of K features measured on N subjects for two conditions, A and B
    % nperms: number of permutations, pretty fast so 10000 is good to get a precise p-value
        
    % OUTPUT:
    % pvals_twotail: equal tail bootstrap p-value for mean(A,1) ~= mean(B,1). 
    % observedA_gt_B: observed value of mean(A,1) - mean(B,1)
    
    disp('start permutation testing')
    [N_a, K] = size(A);
    [N_b, K] = size(B);    
    N = N_a + N_b;
    X = [A;B]; % NxK matrix
    groupAssignmentsOriginal = [zeros(N_a,1); ones(N_b,1)]; % make group assignments for each of the N_a + N_b observations in A and B
    fractionB = N_b / N; % what percentage of the sample makes up B?
    %[~, permutationIndex] = sort(rand(N,nperms),1); % equivalent of randomly permuting without replacement, i.e. randperm(N,N) repeated nperms times ... don't preallocate this b/c it's massive and takes forever
    expectedA_gt_B = zeros(nperms,K); % difference under null hypothesis 
    for P = 1:nperms
        % Compute amount by which mean B exceeds mean A
        groupAssignmentsPermuted = groupAssignmentsOriginal(randperm(N,N)); % shuffle original group assignments without replacement
        expectedA_gt_B(P,:) = nanmean(X(groupAssignmentsPermuted==0,:),1) - nanmean(X(groupAssignmentsPermuted==1,:),1); % A - B null mean
        disp(['Perm ',num2str(P)])
    end

    %% compare null distribution to observed 
    disp('compare null to observed')

    % Compute actual difference between mean B and mean A
    observedA_gt_B = nanmean(A,1) - nanmean(B,1);

    % Compute equal tail p-value for (A - B) > null
    % this should be plotted on the actual A-B
    
    pvals_twotail = NP_TWOTAIL(repmat(observedA_gt_B,[nperms 1]) - expectedA_gt_B,0); % is the observed difference between mean of A and mean of B expected if observations in A and B were shuffled?