%% Read and analyze TRN data
clear all, close all, clc

LineWidth = 3;
FontSize = 15;
MarkerSize = 12;

data_dir = 'Data';

% Read raw data
organism_filenames = regexpi(ls(data_dir), '\s', 'split')';
organism_filenames = organism_filenames(~cellfun('isempty',organism_filenames)); % last entry is empty, remove it
organism_filenames = sort(organism_filenames);
N_org = length(organism_filenames); % number of organisms

% name of each organism
for organism_index = 1:N_org
    index = find(organism_filenames{organism_index}=='/');
    if isempty(index)
        start = 1;
    else
        start = index(end)+1;
    end
    index = strfind(organism_filenames{organism_index},'_regulations.csv');
    stop  = index-1;
    organism_names{organism_index} = organism_filenames{organism_index}(start:stop);
end

% split the data into regulator, target, and effect
effect_names = [];
for organism_index = 1:N_org
    T = readtable([data_dir '/' organism_filenames{organism_index}]); 
    column_headers = T.Properties.VariableNames;

    regulator_column = find(strcmp(column_headers,'TF_name'));
    target_column    = find(strcmp(column_headers,'TG_name'));
    effect_column    = find(strcmp(column_headers,'Role'));

    regulator_list = [];
    target_list    = [];
    effect_list    = []; 
    for i = 1:size(T,1)
        regulator_list{i} = char(T{i,regulator_column});  % all regulators
        target_list{i}    = char(T{i,target_column});     % all targets
        effect_list{i}    = char(T{i,effect_column});     % all effects
    end
    RList{organism_index} = regulator_list; 
    TList{organism_index} = target_list;
    EList{organism_index} = effect_list;

    effect_names = [effect_names unique(EList{organism_index})];
end


% compile the interaction matrices
for organism_index = 1:N_org
    effect_list_unified = {};

    % those that act as activators, repressors, or anything else
    act_index  = find(contains(EList{organism_index},["A","+"]));
    rep_index  = find(contains(EList{organism_index},["R","-"]));
    both_index = find(contains(EList{organism_index},["D"]));
    amb_index = setdiff([1:length(EList{organism_index})],[act_index, rep_index, both_index]);
    
    for i = 1:length(act_index)
        effect_list_unified{act_index(i)} = '+';
    end
    for i = 1:length(rep_index)
        effect_list_unified{rep_index(i)} = '-';
    end
    for i = 1:length(both_index)
        effect_list_unified{both_index(i)} = 'd';
    end
    for i = 1:length(amb_index)
        effect_list_unified{amb_index(i)} = '?';
    end

    interaction_table{organism_index} = table(RList{organism_index}', TList{organism_index}', effect_list_unified', 'VariableNames', {'regulator', 'target', 'effect'});
end

% analyze the interacion table
for i = 1:N_org
    [nr(i,:), Xr(i,:), nr_mod(i,:), Xr_mod(i,:), n(i,:), X(i,:)] = AnalyzeInteractionTable(interaction_table{i}, organism_names{i});
end