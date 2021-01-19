% FADE-SAME: Table S1

% clear
% clc

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set subjects file
subj_file = 'subjects/subjects.xls';
str_cohs  = {'young', 'older', 'middle-aged', 'replication'};

% load subjects file
[num, txt, raw] = xlsread(subj_file);
clear num txt
tab_hdr   = raw(1,:);
subj_data = raw(2:end,:);
num_subj  = size(subj_data,1);
CV_grp    = cell2mat(subj_data(:,end));
tab_hdr(2:4) = {'scanner', 'sex', 'age'};

% get subject groups
num_cohs = numel(str_cohs);
coh_inds = zeros(num_subj,1);
for i = 1:num_subj
    if strncmp(subj_data{i,1},'subA',4)
        if subj_data{i,4} < 50
            coh_inds(i) = 1;    % young AiA
        elseif subj_data{i,4} < 60
            coh_inds(i) = 3;    % middle-aged AiA
        else
            coh_inds(i) = 2;    % older AiA
        end;
    else
        coh_inds(i) = 4;        % yFADE
    end;
end;

% report statistics
for i = 1:num_cohs
    
    % display message
    fprintf('\n-> %s subjects:', str_cohs{i});
    
    % load both groups
    subj_info1 = subj_data(coh_inds==i & CV_grp==1,:);
    subj_info2 = subj_data(coh_inds==i & CV_grp==2,:);
    
    % Info 1: number of subjects
    num_subj1 = size(subj_info1,1);
    num_subj2 = size(subj_info2,1);
    fprintf('\n   - number of subjects: ');
    fprintf('\n     - G1: N = %d;', num_subj1);
    fprintf('\n     - G2: N = %d;', num_subj2);
    
    % Info 2: age range
    y1a = cell2mat(subj_info1(:,strcmp(tab_hdr,'age')));
    y1b = cell2mat(subj_info2(:,strcmp(tab_hdr,'age')));
    fprintf('\n   - age range: ');
    fprintf('\n     - G1: %d-%d yrs;', round(min(y1a)), round(max(y1a)));
    fprintf('\n     - G2: %d-%d yrs;', round(min(y1b)), round(max(y1b)));
    
    % Test 1: mean age (two-sample t-test)
    fprintf('\n   - mean age: ');
    y1a = cell2mat(subj_info1(:,strcmp(tab_hdr,'age')));
    y1b = cell2mat(subj_info2(:,strcmp(tab_hdr,'age')));
    [h, p1, ci, stats1] = ttest2(y1a, y1b, 'tail', 'both', 'vartype', 'unequal');
    fprintf('\n     - G1: %2.2f +- %1.3f yrs;', mean(y1a), std(y1a));
    fprintf('\n     - G2: %2.2f +- %1.3f yrs;', mean(y1b), std(y1b));
    fprintf('\n     - two-sample t-test: ');
    if p1 < 0.001
        fprintf('t = %1.2f, p < 0.001.', stats1.tstat);
    else
        fprintf('t = %1.2f, p = %0.3f.', stats1.tstat, p1);
    end;
    
    % Test 2: sex ratio (Fisher's exact test)
    fprintf('\n   - sex ratio: ');
    sex1 = cell2mat(subj_info1(:,strcmp(tab_hdr,'sex')));
    sex2 = cell2mat(subj_info2(:,strcmp(tab_hdr,'sex')));
    Y2 = [sum(sex1==1), sum(sex1==2);
          sum(sex2==1), sum(sex2==2)];
    [h, p2, stats2] = fishertest(Y2, 'tail', 'both');
    fprintf('\n     - G1: %d : %d m/f;', Y2(1,1), Y2(1,2));
    fprintf('\n     - G2: %d : %d m/f;', Y2(2,1), Y2(2,2));
    fprintf('\n     - Fisher''s exact test: ');
    if p2 < 0.001
        fprintf('OR = %1.2f, p < 0.001.', stats2.OddsRatio);
    else
        fprintf('OR = %1.2f, p = %0.3f.', stats2.OddsRatio, p2);
    end;
    
    % Test 3: scanner ratio (Fisher's exact test)
    fprintf('\n   - scanner ratio: ');
    scan1 = cell2mat(subj_info1(:,strcmp(tab_hdr,'scanner')));
    scan2 = cell2mat(subj_info2(:,strcmp(tab_hdr,'scanner')));
    Y3 = [sum(scan1==1), sum(scan1==2);
          sum(scan2==1), sum(scan2==2)];
    [h, p3, stats3] = fishertest(Y3, 'tail', 'both');
    fprintf('\n     - G1: %d : %d V/S;', Y3(1,1), Y3(1,2));
    fprintf('\n     - G2: %d : %d V/S;', Y3(2,1), Y3(2,2));
    fprintf('\n     - Fisher''s exact test: ');
    if p3 < 0.001
        fprintf('OR = %1.2f, p < 0.001.', stats3.OddsRatio);
    else
        fprintf('OR = %1.2f, p = %0.3f.', stats3.OddsRatio, p3);
    end;
    fprintf('\n');
    
end;
fprintf('\n');