% Data analysis for FADE-SAME
% _
% This script performs the entire data analysis for the FADE-SAME study.
% 
% written by Joram Soch <Joram.Soch@DZNE.de>, 11/01/2021, 16:53


clc
clear

%%% Step 0: analysis parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% directories and subjects
stud_dir   =  'C:\Joram\projects\DZNE\FADE\sharing\FADE_SAME\';
subj_file  =  'subjects/subjects.xls';
subj_files = {'subjects_FADE_young_G1.xls', 'subjects_FADE_young_G2.xls';
              'subjects_FADE_older_G1.xls', 'subjects_FADE_older_G2.xls';
              'subjects_yFADE_G1.xls',      'subjects_yFADE_G2.xls'};

% create analysis super-folders
if ~exist(strcat(stud_dir,'FADE_scores/'),'dir')
    mkdir(strcat(stud_dir,'FADE_scores/'));
end;
if ~exist(strcat(stud_dir,'SPM_analyses/'),'dir')
    mkdir(strcat(stud_dir,'SPM_analyses/'));
end;

% specify analyses to perform
thx2do = [1, 2, 3, 4];


%%% Step 1: Generate first-level contrasts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(1,thx2do)

% compute novelty contrast
contrast_meta_batch(stud_dir, subj_file, 'GLM_TD_1th-a', [1  0 -1], 'con_0001.nii');

% compute memory contrast
contrast_meta_batch(stud_dir, subj_file, 'GLM_TD_1th-a', [0  1  0], 'con_0002.nii');
    
end;


%%% Step 2: Estimate second-level models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(2,thx2do)

% perform one-sample t-tests
GLM_name  = 'GLM_TD_1th-a';
con_names ={'con_0001.nii', 'con_0002.nii'};
cov_names ={};
job_list1 = cell(numel(subj_files)*numel(con_names),1);
for i1 = 1:size(subj_files,1)
    for i2 = 1:size(subj_files,2)
        for j = 1:numel(con_names)
            subj_list = strcat('subjects/',subj_files{i1,i2});
            img_path  = strcat('beta_images/*_',GLM_name,'_',con_names{j});
            stat_suff = subj_files{i1,i2};
            stat_suff = stat_suff(strfind(stat_suff,'subjects_')+9:strfind(stat_suff,'.xls')-1);
            if j == 1, stat_suff = strcat(GLM_name,'_novelty_',stat_suff); end;
            if j == 2, stat_suff = strcat(GLM_name,'_memory_', stat_suff); end;
            job_list1{(i1-1)*4+(i2-1)*2+j*1} = ...
                fact_des_ttest1(stud_dir, subj_list, img_path, cov_names, stat_suff);
        end;
    end;
end;
job_list2 = cell(numel(con_names),1);
for j = 1:numel(con_names)
    subj_list = strcat('subjects/','subjects_FADE_young.xls');
    img_path  = strcat('beta_images/*_',GLM_name,'_',con_names{j});
    stat_suff = 'FADE_young_all';
    if j == 1, stat_suff = strcat(GLM_name,'_novelty_',stat_suff); end;
    if j == 2, stat_suff = strcat(GLM_name,'_memory_', stat_suff); end;    
    job_list2{j} = fact_des_ttest1(stud_dir, subj_list, img_path, cov_names, stat_suff);
end;
job_list = [job_list1; job_list2];
spm_exec_mult_jobs(job_list);

% create reference maps
con = [1, 2, 3];
FWE = true;
p   = 0.05;
k   = 10;
for i = 1:numel(job_list)
    SPM_dir = fileparts(job_list{i});
    SPM_mat = strcat(SPM_dir,'/','SPM.mat');
    load(SPM_mat);
    for j = 1:numel(con)
        if FWE
            filename = sprintf('con_%04d_FWE_%s_%d.nii', con(j), num2str(p), k);
        else
            filename = sprintf('con_%04d_unc_%s_%d.nii', con(j), num2str(p), k);
        end;
        spm_save_thr_SPM(SPM, con(j), FWE, p, k, filename);
        cd(strcat(stud_dir,'FADE_SAME/'));
    end;
end;
    
end;


%%% Step 3: Calculate FADE and SAME scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(3,thx2do)

% specify analyses
GLM_name  =  'GLM_TD_1th-a';
ref_dirs  = {'GLM_TD_1th-a_~_FADE_young_G1', 'GLM_TD_1th-a_~_FADE_young_G2';
             'GLM_TD_1th-a_~_yFADE_G1',      'GLM_TD_1th-a_~_yFADE_G2'};
ref_maps  = {'con_0002_FWE_0.05_10.nii',     'con_0001_FWE_0.05_10.nii'};
con_vecs  = {[+1, 0, -1], [0, 1, 0]};
contrasts = {'novelty', 'memory'};
scores    = {'FADE-pos', 'SAME'};

% calculate scores (FADE only)
subj_file = 'subjects/subjects_FADE.xls';
res_suff  = 'FADE_only';
calc_FADE_SAME_scores(stud_dir, subj_file, GLM_name, ref_dirs, ref_maps, con_vecs, contrasts, scores, res_suff)

% calculate scores (FADE & yFADE)
subj_file = 'subjects/subjects.xls';
res_suff  = 'FADE_yFADE';
calc_FADE_SAME_scores(stud_dir, subj_file, GLM_name, ref_dirs, ref_maps, con_vecs, contrasts, scores, res_suff)

% calculate scores (reference to FADE)
ref_dirsA = [ref_dirs(1,:); ref_dirs(1,:)];
res_suff  = 'ref_FADE';
calc_FADE_SAME_scores(stud_dir, subj_file, GLM_name, ref_dirsA, ref_maps, con_vecs, contrasts, scores, res_suff)

% calculate scores (reference to yFADE)
ref_dirsB = [ref_dirs(2,:); ref_dirs(2,:)];
res_suff  = 'ref_yFADE';
calc_FADE_SAME_scores(stud_dir, subj_file, GLM_name, ref_dirsB, ref_maps, con_vecs, contrasts, scores, res_suff)

% calculate scores (reference to FADE group 1)
subj_file = 'subjects/subjects_FADE_older.xls';
ref_dirs1 = [ref_dirs(:,1), ref_dirs(:,1)];
res_suff  = 'ref_FADE_G1';
calc_FADE_SAME_scores(stud_dir, subj_file, GLM_name, ref_dirs1, ref_maps, con_vecs, contrasts, scores, res_suff)

% calculate scores (reference to FADE group 2)
ref_dirs2 = [ref_dirs(:,2), ref_dirs(:,2)];
res_suff  = 'ref_FADE_G2';
calc_FADE_SAME_scores(stud_dir, subj_file, GLM_name, ref_dirs2, ref_maps, con_vecs, contrasts, scores, res_suff)

% calculate scores (reference to all FADE subjects)
ref_dirs12={'GLM_TD_1th-a_~_FADE_young_all', 'GLM_TD_1th-a_~_FADE_young_all';
            'GLM_TD_1th-a_~_yFADE_G1',      'GLM_TD_1th-a_~_yFADE_G2'};
res_suff  = 'ref_FADE_all';
calc_FADE_SAME_scores(stud_dir, subj_file, GLM_name, ref_dirs12, ref_maps, con_vecs, contrasts, scores, res_suff)

clear ref_dirsA ref_dirsB ref_dirs1 ref_dirs2 ref_dirs12

end;


%%% Step 4: Reproduce Figures and Tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(4,thx2do)
    
% Tables 2-3
% Figures 3-6
% Tables S3-S5
% Figures S2-S7

end;