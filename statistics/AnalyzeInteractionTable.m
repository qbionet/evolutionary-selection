function [nr, Xr, nr_mod, Xr_mod, n, X] = AnalyzeInteractionTable(interaction_table, organism_name)
       
    %% RAW DATA
    M = interaction_table;
    for i = 1:size(M,1)
        regulator_list(i) = lower(string(M{i,1}));  % all regulators
        target_list(i)    = lower(string(M{i,2}));  % all targets
        effect_list(i)    = string(M{i,3});         % all effects
    end
    regulators = unique(regulator_list); % list of all TF regulators
    targets    = unique(target_list);    % list of all targets genes
    
    
    %% REMOVAL OF AMBIGUOUS INTERACTIONS
    index = find(strcmp(effect_list,'?'));
    removed_interactions = M(index,:);
    regulator_list(index) = [];
    target_list(index) = [];
    effect_list(index) = [];
    
    regulators = unique(regulator_list); % list of all TF regulators
    targets    = unique(target_list);    % list of all targets genes   
    
    
    %% REGULATORS 
    for i = 1:length(regulators)
        index = find(strcmp(regulators(i),regulator_list));
        if sum(strcmp(effect_list(index),'+'))
            pos_neg_regulator(i,1) = 1; % acts as an activator
        else
            pos_neg_regulator(i,1) = 0;
        end
        if sum(strcmp(effect_list(index),'-'))
            pos_neg_regulator(i,2) = 1; % acts as a repressor
        else
            pos_neg_regulator(i,2) = 0;
        end
        if sum(strcmp(effect_list(index),'d'))
            pos_neg_regulator(i,1:2) = 1; % acts as both
        end
    end
    regulators_pos  = regulators(pos_neg_regulator(:,1)==1); % TFs acting as activators
    regulators_neg  = regulators(pos_neg_regulator(:,2)==1); % TFs acting as repressors
    
    % Find the autoregulated regulators (auto_pos & auto_neg)
    pos_count = 0;
    neg_count = 0;
    for i = 1:length(regulator_list)
        if strcmp(regulator_list(i),target_list(i))
            if strcmp(effect_list(i),'+')
                pos_count = pos_count+1;
                auto_pos(pos_count) = regulator_list(i);
            elseif strcmp(effect_list(i),'-')
                neg_count = neg_count+1;
                auto_neg(neg_count) = regulator_list(i);
            end
        end
    end
    if ~exist('auto_pos')
        auto_pos = [];
    end
    if ~exist('auto_neg')
        auto_neg = [];
    end

    % auto_both: simultaneously +FB and -FB
    if isempty(auto_pos) | isempty(auto_neg)
        auto_both = 'NULL';
    else
        auto_both = intersect(auto_pos, auto_neg);
    end

    % naive reference sample size and frequencies 
    FBcount(1,:) = [length(auto_neg) length(regulators)-length(auto_pos)-length(auto_neg) length(auto_pos)];

    
    
    %% TARGETS

    % Find positively/negatively regulated targets
    for i = 1:length(targets)
        index = find(strcmp(targets(i),target_list));
        if sum(strcmp(effect_list(index),'+'))
            pos_neg_target(i,1) = 1;
        else
            pos_neg_target(i,1) = 0;
        end
        if sum(strcmp(effect_list(index),'-'))
            pos_neg_target(i,2) = 1;
        else
            pos_neg_target(i,2) = 0;
        end

    end
    targets_pos_strict  = targets(pos_neg_target(:,1)-pos_neg_target(:,2)==1); % strictly positively regulated
    targets_neg_strict  = targets(pos_neg_target(:,2)-pos_neg_target(:,1)==1); % strictly negatively regulated

    
    %% Samples of interest

    % Among those that are negatively regulated only, what percentage of
    % regulators is positively/negatively self-regulated
    regulators_of_purely_repressed_targets = [];
    negOnly_regulators = [];
    frac_neg_negFB = [];
    frac_neg_posFB = [];
    frac_neg_noFB  = [];
    for i = 1:length(targets_neg_strict)
        index = find(strcmp(targets_neg_strict(i),target_list)); % original list with target that is only repressed
        actual_regulators = regulator_list(index); % regulators of the target that is only repressed
        regulators_of_purely_repressed_targets = [regulators_of_purely_repressed_targets actual_regulators];

        % fraction of regulators of the target that are self-activated
        if ~isempty(auto_pos)
            frac_neg_posFB(i) = sum(ismember(actual_regulators,auto_pos))/length(actual_regulators);
        else
            frac_neg_posFB(i) = 0;
        end
    
        % fraction of regulators of the target that are self-repressed
        if ~isempty(auto_neg)
            frac_neg_negFB(i) = sum(ismember(actual_regulators,auto_neg))/length(actual_regulators);
        else
            frac_neg_negFB(i) = 0;
        end
    
        % fraction of regulators of the target that are not self-regulated
        frac_neg_noFB(i)  = 1 - frac_neg_posFB(i) - frac_neg_negFB(i);
    
        % number of regulators of the target
        negOnly_regulators(i) = length(actual_regulators);
    end


    % probability of negative/no/positive self-regulators among regulators of
    % targets with strict repression    
    if isempty(negOnly_regulators*frac_neg_negFB')
        feedback_matrix(2,1) = 0;
    else
        feedback_matrix(2,1) = negOnly_regulators*frac_neg_negFB';
    end

    if isempty(negOnly_regulators*frac_neg_noFB')
        feedback_matrix(2,2) = 0;
    else
        feedback_matrix(2,2) = negOnly_regulators*frac_neg_noFB';
    end

    if isempty(negOnly_regulators*frac_neg_posFB')
        feedback_matrix(2,3) = 0;
    else
        feedback_matrix(2,3) = negOnly_regulators*frac_neg_posFB';
    end
    FBcount(2,:) = feedback_matrix(2,:);
    feedback_matrix(2,:) = feedback_matrix(2,:)/sum(negOnly_regulators);

    
    % Among those that are positively regulated only, what percentage is self-regulated
    regulators_of_purely_activated_targets = [];
    posOnly_regulators = [];
    frac_pos_negFB = [];
    frac_pos_posFB = [];
    frac_pos_noFB  = [];
    for i = 1:length(targets_pos_strict)
        index = find(strcmp(targets_pos_strict(i),target_list)); % original list with target that is only activated
        actual_regulators = regulator_list(index); % regulators of the target that is only activated
        regulators_of_purely_activated_targets = [regulators_of_purely_activated_targets actual_regulators];
    
        % fraction of regulators of the target that are self-activated
        if ~isempty(auto_pos)
            frac_pos_posFB(i) = sum(ismember(actual_regulators,auto_pos))/length(actual_regulators);
        else
            frac_pos_posFB(i) = 0;
        end

        % fraction of regulators of the target that are self-repressed
        if ~isempty(auto_neg)
            frac_pos_negFB(i) = sum(ismember(actual_regulators,auto_neg))/length(actual_regulators);
        else
            frac_pos_negFB(i) = 0;
        end

        % fraction of regulators of the target that are not self-regulated
        frac_pos_noFB(i)  = 1 - frac_pos_posFB(i) - frac_pos_negFB(i);
    
        % number of regulators of the target
        posOnly_regulators(i) = length(actual_regulators);
    end
    
    % probability of negative/no/positive self-regulators among regulators of
    % targets with strict activation
    if isempty(posOnly_regulators*frac_pos_negFB')
        feedback_matrix(3,1) = 0;
    else
        feedback_matrix(3,1) = posOnly_regulators*frac_pos_negFB';
    end

    if isempty(posOnly_regulators*frac_pos_noFB')
        feedback_matrix(3,2) = 0;
    else
        feedback_matrix(3,2) = posOnly_regulators*frac_pos_noFB';
    end

    if isempty(posOnly_regulators*frac_pos_posFB')
        feedback_matrix(3,3) = 0;
    else
        feedback_matrix(3,3) = posOnly_regulators*frac_pos_posFB';
    end
    FBcount(3,:) = feedback_matrix(3,:);
    feedback_matrix(3,:) = feedback_matrix(3,:)/sum(posOnly_regulators);
    

    %% Modified baseline
    % Among those that are regulated, what percentage is self-regulated
    regulators_of_all_targets = [];
    all_regulators = [];
    frac_all_negFB = [];
    frac_all_posFB = [];
    frac_all_noFB  = [];
    for i = 1:length(targets)
        index = find(strcmp(targets(i),target_list)); % original list with target
        actual_regulators = regulator_list(index); % regulators of the target
        regulators_of_all_targets = [regulators_of_all_targets actual_regulators];
    
        % fraction of regulators of the target that are self-activated
        if ~isempty(auto_pos)
            frac_all_posFB(i) = sum(ismember(actual_regulators,auto_pos))/length(actual_regulators);
        else
            frac_all_posFB(i) = 0;
        end
    
        % fraction of regulators of the target that are self-repressed
        if ~isempty(auto_neg)
            frac_all_negFB(i) = sum(ismember(actual_regulators,auto_neg))/length(actual_regulators);
        else
            frac_all_negFB(i) = 0;
        end

        % fraction of regulators of the target that are not self-regulated
        frac_all_noFB(i)  = 1 - frac_all_posFB(i) - frac_all_negFB(i);
    
        % number of regulators of the target
        all_regulators(i) = length(actual_regulators);
    end
    
    % probability of negative/no/positive self-regulators among regulators of
    % all targets
    if isempty(all_regulators*frac_all_negFB')
        feedback_matrix(4,1) = 0;
    else
        feedback_matrix(4,1) = all_regulators*frac_all_negFB';
    end

    if isempty(all_regulators*frac_all_noFB')
        feedback_matrix(4,2) = 0;
    else
        feedback_matrix(4,2) = all_regulators*frac_all_noFB';
    end

    if isempty(all_regulators*frac_all_posFB')
        feedback_matrix(4,3) = 0;
    else
        feedback_matrix(4,3) = all_regulators*frac_all_posFB';
    end
    FBcount(4,:) = feedback_matrix(4,:);
    feedback_matrix(4,:) = feedback_matrix(4,:)/sum(all_regulators); 


    

    nr     = [sum(FBcount(1,:)) sum(FBcount(1,:))];
    Xr     = [FBcount(1,1) FBcount(1,3)];
    nr_mod = [sum(FBcount(4,:)) sum(FBcount(4,:))];
    Xr_mod = [FBcount(4,1) FBcount(4,3)];
    n      = [sum(FBcount(2,:)) sum(FBcount(3,:))];
    X      = [FBcount(2,1) FBcount(3,3)];
end