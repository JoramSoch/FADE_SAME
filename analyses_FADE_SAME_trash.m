%%% novelty/memory reference maps, 26/11/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ../project_directories_BS_JS.mat
% SPM_dirs = {'MS_FADE_04_FADE_GLM_1a_novelty_subjects_all_2020_11_05_young_all', ...
%             'MS_FADE_04_FADE_GLM_1a_memory_subjects_all_2020_11_05_young_all'};
% con = [1, 2, 3];
% FWE = true;
% p   = 0.05;
% k   = 10;
% for i = 1:numel(SPM_dirs)
%     SPM_dir = strcat(stat_dir,'group_statistics','/',SPM_dirs{i},'/');
%     load(strcat(SPM_dir,'SPM.mat'));
%     for j = 1:numel(con)
%         if FWE
%             filename = sprintf('con_%04d_FWE_%s_%d.nii', con(j), num2str(p), k);
%         else
%             filename = sprintf('con_%04d_unc_%s_%d.nii', con(j), num2str(p), k);
%         end;
%         spm_save_thr_SPM(SPM, con(j), FWE, p, k, filename);
%     end;
% end;
% cd(tool_dir);


%%% novelty/memory reference maps, 26/11/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_11_05_young.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_04.mat';
% GLM_name  = 'GLM_TD_1th-a';
% con_names ={'con_0001.nii', 'con_0002.nii'};
% cov_names ={};
% for j = 1:numel(con_names)
%     img_path = strcat('group_statistics/+/*_s6w_',GLM_name,'_',con_names{j});
%     [folder, file, ext] = fileparts(subj_file);
%     if j == 1, stat_suff = strcat('FADE_GLM_1a_novelty_',file,'_all'); end;
%     if j == 2, stat_suff = strcat('FADE_GLM_1a_memory_',file,'_all');  end;
%     fact_des_ttest1(subj_file, MS_file, img_path, cov_names, stat_suff);
% end;
% spm_exec_mult_jobs;


%%% analyze subjects, 20/11/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file  = 'subjects_all_2020_11_05_young_very_old.xls';
% analyze_subjects(subj_file);
% delay_file = 'FADE_delay_all_subjects_2020_11_05_JS.xls';
% analyze_delays(subj_file, delay_file);


%%% scanner x gender x age ANOVA, 12/11/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_11_05_young_very_old.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_04.mat';
% img_paths ={'group_statistics/+/*_s6w_GLM_TD_1th-a_con_0001.nii', ...
%             'group_statistics/+/*_s6w_GLM_TD_1th-a_con_0002.nii'};
% var_names ={'scanner', 'sex', 'age_group'};
% cov_names ={};
% stat_suff ={'new_GLM_1a_novel-master_by_scanner_gender_age', ...
%             'new_GLM_1a_memory_by_scanner_gender_age'};
% for i = 1:numel(stat_suff)
%     fact_des_full_fact(subj_file, MS_file, img_paths(i), var_names, cov_names, stat_suff{i})
% end;
% spm_exec_mult_jobs;


%%% novelty/memory reference maps, 11/11/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ../project_directories_BS_JS.mat
% SPM_dirs = {'MS_FADE_04_FADE_GLM_1a_novelty_subjects_all_2020_11_05_young_G1', ...
%             'MS_FADE_04_FADE_GLM_1a_novelty_subjects_all_2020_11_05_young_G2', ...
%             'MS_FADE_04_FADE_GLM_1a_memory_subjects_all_2020_11_05_young_G1', ...
%             'MS_FADE_04_FADE_GLM_1a_memory_subjects_all_2020_11_05_young_G2'};
% con = [1, 2, 3];
% FWE = true;
% p   = 0.05;
% k   = 10;
% for i = 1:numel(SPM_dirs)
%     SPM_dir = strcat(stat_dir,'group_statistics','/',SPM_dirs{i},'/');
%     load(strcat(SPM_dir,'SPM.mat'));
%     for j = 1:numel(con)
%         if FWE
%             filename = sprintf('con_%04d_FWE_%s_%d.nii', con(j), num2str(p), k);
%         else
%             filename = sprintf('con_%04d_unc_%s_%d.nii', con(j), num2str(p), k);
%         end;
%         spm_save_thr_SPM(SPM, con(j), FWE, p, k, filename);
%     end;
% end;
% cd(tool_dir);


%%% novelty/memory reference maps, 10/11/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_files = {'subjects_all_2020_11_05_young_G1.xls', ...
%               'subjects_all_2020_11_05_young_G2.xls'};
% MS_file    =  'FADE_model_spaces/MS_FADE_04.mat';
% GLM_name   =  'GLM_TD_1th-a';
% con_names  = {'con_0001.nii', 'con_0002.nii'};
% cov_names  = {};
% for i = 1:numel(subj_files)
%     for j = 1:numel(con_names)
%         img_path = strcat('group_statistics/+/*_s6w_',GLM_name,'_',con_names{j});
%         [folder, file, ext] = fileparts(subj_files{i});
%         if j == 1, stat_suff = strcat('FADE_GLM_1a_novelty_',file); end;
%         if j == 2, stat_suff = strcat('FADE_GLM_1a_memory_',file);  end;
%         fact_des_ttest1(subj_files{i}, MS_file, img_path, cov_names, stat_suff);
%     end;
% end;
% spm_exec_mult_jobs;


%%% young vs. old comparisons, 10/11/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_11_05_young_very_old.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_04.mat';
% img_path  ={'group_statistics/+/*_s6w_GLM_TD_1th-a_con_0001.nii', ...
%             'group_statistics/+/*_s6w_GLM_TD_1th-a_con_0002.nii'};
% var_name  = 'age_group';
% cov_names ={};
% stat_suff ={'GLM_1a_memory_young_vs_very_old', 'GLM_1a_novelty-master_young_vs_very_old'};
% for i = 1:numel(img_path)
%     fact_des_ttest2(subj_file, MS_file, img_path{i}, var_name, cov_names, stat_suff{i});
% end;
% spm_exec_mult_jobs;


%%% partition subjects for FADE/SAME, 10/11/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_files = {'subjects_all_2020_11_05_young.xls', 
%               'subjects_all_2020_11_05_middle_aged.xls', 
%               'subjects_all_2020_11_05_very_old.xls'};
% for i = 1:numel(subj_files)
%     partition_subjects(subj_files{i})
% end;


%%% repeat fMRI analyses (259 subjects), 05/11/2020 %%%%%%%%%%%%%%%%%%%%%%%

% % subject files
% subj_files = {'subjects_all_2020_11_05.xls';
%               'subjects_all_2020_11_05_young.xls';
%               'subjects_all_2020_11_05_old.xls'};
% subj_suffs = {'', '_young', '_old'};
% stat_pref  =  'new_'; % '259_subj_';
% 
% % 0
% MS_file = 'FADE_model_spaces/MS_FADE_04_00_0.mat';
% family.mods = [1 1 1 1 2 2 2 2];
% family.fams = {'GLMs_PE', 'GLMs_TD'};
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat(stat_pref, '00_0_PE_vs_TD', subj_suffs{i})));
%     MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
%     load(strcat(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man.dir{1},'BMS.mat'));
%     MS_SMM_BMS_fams(BMS)
% end;
% 
% % 1
% MS_file = 'FADE_model_spaces/MS_FADE_04_00_0.mat';
% family.mods = [1 1 2 2 1 1 2 2];
% family.fams = {'GLMs_00', 'GLMs_0'};
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat(stat_pref, '00_0_00_vs_0', subj_suffs{i})));
%     MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
%     load(strcat(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man.dir{1},'BMS.mat'));
%     MS_SMM_BMS_fams(BMS)
% end;
% 
% % 2
% MS_file = 'FADE_model_spaces/MS_FADE_04_00_0.mat';
% family.mods = [1 2 1 2 1 2 1 2];
% family.fams = {'GLMs_x1', 'GLMs_x2'};
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat(stat_pref, '00_0_x1_vs_x2', subj_suffs{i})));
%     MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
%     load(strcat(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man.dir{1},'BMS.mat'));
%     MS_SMM_BMS_fams(BMS)
% end;
% 
% % 3
% % MS_file = 'FADE_model_spaces/MS_FADE_04_TD_0_1_2.mat';
% % cvLFE_meta_batch(subj_files{1}, MS_file, 's6w', [1 0 2 2 2 2 2 2 3 3 3], {'GLMs_TD_0', 'GLMs_TD_1', 'GLMs_TD_2'});
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_0_1_2_fams.mat';
% family.mods = [1 2 2];
% family.fams = {'GLMs_TD_0', 'GLMs_TD_12'};
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLFE.nii', strcat(stat_pref, 'TD_0_1_2_12_vs_0', subj_suffs{i})));
%     MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
%     load(strcat(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man.dir{1},'BMS.mat'));
%     MS_SMM_BMS_fams(BMS)
% end;
% 
% % 4a
% % subjects_all_2020_11_05_GLMs
% % subj_file = 'subjects_all_2020_11_05_GLM_TD_5.xls';
% % MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_3_5.mat';
% % cvLFE_meta_batch(subj_file, MS_file, 's6w', [1 2], {'GLMs_TD_3', 'GLMs_TD_5'});
% % MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_0-5_fams.mat';
% % cvLBF_meta_batch(subj_file, MS_file, 's6w', [2 1 0 0 0], 'MA_cvLFE.nii', 'GLMs_TD_1_vs_0');
% % cvLBF_meta_batch(subj_file, MS_file, 's6w', [2 0 1 0 0], 'MA_cvLFE.nii', 'GLMs_TD_2_vs_0');
% % cvLBF_meta_batch(subj_file, MS_file, 's6w', [2 0 0 1 0], 'MA_cvLFE.nii', 'GLMs_TD_3_vs_0');
% % cvLBF_meta_batch(subj_file, MS_file, 's6w', [2 0 0 0 1], 'MA_cvLFE.nii', 'GLMs_TD_5_vs_0');
% 
% % 4b
% subj_GLM5 ={'subjects_all_2020_11_05_GLM_TD_5.xls';
%             'subjects_all_2020_11_05_GLM_TD_5_young.xls';
%             'subjects_all_2020_11_05_GLM_TD_5_old.xls'};
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_0-5_fams.mat';
% img_paths ={'model_selection/+/*_s6w_GLMs_TD_1_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s6w_GLMs_TD_2_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s6w_GLMs_TD_3_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s6w_GLMs_TD_5_vs_0_MC_cvLBF.nii'};
% var_names ={};
% cov_names ={'sex'};
% stat_suff = 'LBF_GLMs_TD_1235_vs_0';
% for i = 1:numel(subj_GLM5)
%     fact_des_full_fact(subj_GLM5{i}, MS_file, img_paths, var_names, cov_names, strcat(stat_pref, stat_suff, subj_suffs{i}));
% end;
% spm_exec_mult_jobs;
% 
% % 5
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_1_2.mat';
% family.mods = [1 1 1 1 1 1 2 2 2];
% family.fams = {'GLMs_TD_1', 'GLMs_TD_2'};
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat(stat_pref, 'TD_1_2_1_vs_2', subj_suffs{i})));
%     MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
%     load(strcat(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man.dir{1},'BMS.mat'));
%     MS_SMM_BMS_fams(BMS)
% end;
% 
% % 6
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_1.mat';
% family.mods = [1 1 1 2 2 2];
% family.fams = {'GLMs_TD_1dd', 'GLMs_TD_1th'};
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat(stat_pref, 'TD_1_dd_vs_th', subj_suffs{i})));
%     MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
%     load(strcat(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man.dir{1},'BMS.mat'));
%     MS_SMM_BMS_fams(BMS)
% end;
% 
% % 7, 8, 9
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_2.mat';
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat(stat_pref, 'TD_2_nf_nr_ns', subj_suffs{i})));
%     spm_jobman('run', matlabbatch);
% end;
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_1dd.mat';
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat(stat_pref, 'TD_1dd_ip_cp_lr', subj_suffs{i})));
%     spm_jobman('run', matlabbatch);
% end;
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_1th.mat';
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat(stat_pref, 'TD_1th_l_a_s', subj_suffs{i})));
%     spm_jobman('run', matlabbatch);
% end;
% 
% % 10a, 10b
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_WM3.mat';
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat(stat_pref, 'TD_WM3', subj_suffs{i})));
%     spm_jobman('run', matlabbatch);
% end;
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_WM2.mat';
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat(stat_pref, 'TD_WM2', subj_suffs{i})));
%     spm_jobman('run', matlabbatch);
% end;
% 
% % contrast images
% subj_file = 'subjects_all_2020_11_05.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_04.mat';
% contrast_meta_batch(subj_file, MS_file, 's6w', 'GLM_TD_1th-a', [1  0 -1], 'con_0001.nii');
% contrast_meta_batch(subj_file, MS_file, 's6w', 'GLM_TD_1th-a', [0  1  0], 'con_0002.nii');
% contrast_meta_batch(subj_file, MS_file, 's6w', 'GLM_TD_2nf', [1  1 -2], 'con_0001.nii');
% contrast_meta_batch(subj_file, MS_file, 's6w', 'GLM_TD_2nf', [1 -1  0], 'con_0002.nii');
% 
% % A, B, C, D
% MS_file   = 'FADE_model_spaces/MS_FADE_04.mat';
% GLM_names ={'GLM_TD_1th-a', 'GLM_TD_1th-a', 'GLM_TD_2nf', 'GLM_TD_1dd-ip'};
% con_names ={'con_0001.nii', 'beta_0002.nii', 'con_0002.nii', 'beta_0002.nii'};
% cov_names ={};
% for i = 1:numel(subj_files)
%     for j = 1:numel(GLM_names)
%         img_path = strcat('group_statistics/+/*_s6w_',GLM_names{j},'_',con_names{j});
%         if j == 1, stat_suff = strcat(stat_pref,GLM_names{j},'_novelty-master',subj_suffs{i}); end;
%         if j >  1, stat_suff = strcat(stat_pref,GLM_names{j},'_memory', subj_suffs{i});  end;
%         fact_des_ttest1(subj_files{i}, MS_file, img_path, cov_names, stat_suff);
%     end;
% end;
% spm_exec_mult_jobs;


%%% repeat behavioral analyses (259 subjects), 05/11/2020 %%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_11_05.xls';
% analyze_subjects(subj_file);
% behavioral_data(subj_file);
% bhvr_logfile_data(subj_file);
% bhvr_memorability(subj_file);
% 
% delay_file = 'FADE_delay_all_subjects_2020_11_05_JS.xls';
% analyze_delays(subj_file, delay_file);


%%% one-sample t-tests (LOO), 03/11/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file  =  'subjects_all_2020_10_30_young.xls';
% MS_file    =  'FADE_model_spaces/MS_FADE_04.mat';
% GLM_name   =  'GLM_TD_1th-a';
% con_names  = {'con_0001.nii', 'con_0002.nii'};
% cov_names  = {};
% for j = 1:numel(con_names)
%     img_path = strcat('group_statistics/+/*_s6w_',GLM_name,'_',con_names{j});
%     [folder, file, ext] = fileparts(subj_file);
%     if j == 1, stat_suff = strcat('FADE_GLM_1a_novelty_',file,'_LOO'); end;
%     if j == 2, stat_suff = strcat('FADE_GLM_1a_memory_',file,'_LOO');  end;
%     fact_des_ttest1_LOO(subj_file, MS_file, img_path, cov_names, stat_suff);
% end;


%%% two-sample t-tests, 24/08/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_07_27.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_04.mat';
% img_path  ={'group_statistics/+/*_s6w_GLM_TD_1th-a_con_0001.nii', ...
%             'group_statistics/+/*_s6w_GLM_TD_1th-a_con_0002.nii'};
% var_name  = 'age_group';
% cov_names ={};
% stat_suff ={'GLM_1a_memory_young_vs_old', 'GLM_1a_novelty-master_young_vs_old'};
% for i = 1:numel(img_path)
%     fact_des_ttest2(subj_file, MS_file, img_path{i}, var_name, cov_names, stat_suff{i});
% end;
% spm_exec_mult_jobs;


%%% multiple regressions, 03/08/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file =  'subjects_all_2020_07_27_T1_imgs_FADE_JS.xls';
% MS_file   =  'FADE_model_spaces/MS_VBM_s8_GM.mat';
% img_path  =  'structural_analyses/+/rsmwp1*_T1-MPRAGE_1.nii';
% cov_names = {'novelty_FADE-pos', 'novelty_SAMe', 'memory_FADE-pos', 'memory_SAMe'};
% for i = 1:numel(cov_names)
%     var_names = {'young', strcat(cov_names{i},'_young'), strcat(cov_names{i},'_old')};
%     mean_cent = [false, false, false];
%     stat_suff = strcat('vs_',cov_names{i},'_with_age');
%     fact_des_mult_reg(subj_file, MS_file, img_path, var_names, mean_cent, stat_suff);
% end;
% spm_exec_mult_jobs;

% subj_file =  'subjects_all_2020_07_27_T1_imgs_FADE.xls';
% MS_file   =  'FADE_model_spaces/MS_VBM_s8_GM.mat';
% img_path  =  'structural_analyses/+/rsmwp1*_T1-MPRAGE_1.nii';
% var_names = {'novelty_FADE-pos', 'novelty_SAMe', 'memory_FADE-pos', 'memory_SAMe'};
% mean_cent = [true, true, true, true];
% for i = 1:numel(var_names)
%     stat_suff = strcat('vs_',var_names{i});
%     fact_des_mult_reg(subj_file, MS_file, img_path, var_names(i), mean_cent(i), stat_suff);
% end;
% spm_exec_mult_jobs;


%%% Hippocampus GM volume extraction, 03/08/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%

% load ../project_directories_BS_JS.mat
% subj_file = 'subjects_all_2020_07_27_T1_imgs.xls';
% MS_name   = 's0_GM';
% GM_img    = 'mwp1*_T1-MPRAGE_1.nii';
% ROI_imgs  ={'sHippo_CA_Sub.nii';
%             'sHippocampus_AAL.nii'};
% for j = 1:numel(ROI_imgs)
%     ROI_imgs{j} = strcat(tool_dir,'ROI_images/',ROI_imgs{j});
% end;
% calc_GM_group(subj_file, MS_name, GM_img, ROI_imgs, 'HC');


%%% FADE and SAME score extraction, 03/08/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FADE_SAMe_mega_batch; % 2x: 2020_08_03_a/b
% behavioral_data('subjects_all_2020_07_27.xls');


%%% multiple regressions for BPM, 29/07/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_07_27_T1_imgs.xls';
% var_names = {'age'};
% mean_cent = [false];
% % var_names = {'old', 'age'};
% % mean_cent = [false, false];
% 
% MS_file   = 'FADE_model_spaces/MS_FADE_04.mat';
% img_path  = 'group_statistics/+/*_s6w_GLM_TD_1th-a_con_0001.nii';
% stat_suff = 'BPM_GLM_1a_novelty_vs_age';
% fact_des_mult_reg(subj_file, MS_file, img_path, var_names, mean_cent, stat_suff);
% 
% MS_file   = 'FADE_model_spaces/MS_FADE_04.mat';
% img_path  = 'group_statistics/+/*_s6w_GLM_TD_1th-a_con_0002.nii';
% stat_suff = 'BPM_GLM_1a_memory_vs_age';
% fact_des_mult_reg(subj_file, MS_file, img_path, var_names, mean_cent, stat_suff);
% 
% MS_file   = 'FADE_model_spaces/MS_VBM_s8_GM.mat';
% img_path  = 'structural_analyses/+/rsmwp1*_T1-MPRAGE_1.nii';
% stat_suff = 'BPM_s8_GM_vs_age';
% fact_des_mult_reg(subj_file, MS_file, img_path, var_names, mean_cent, stat_suff);
% spm_exec_mult_jobs;


%%% novelty/memory reference maps, 23/07/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_files = {'subjects_all_2020_07_23_young_G1.xls', 'subjects_all_2020_07_23_young_G2.xls', ...
%               'subjects_all_2020_07_23_old_G1.xls',   'subjects_all_2020_07_23_old_G2.xls'};
% MS_file    =  'FADE_model_spaces/MS_FADE_04.mat';
% GLM_name   =  'GLM_TD_1th-a';
% con_names  = {'con_0001.nii', 'con_0002.nii'};
% cov_names  = {};
% for i = 1:numel(subj_files)
%     for j = 1:numel(con_names)
%         img_path = strcat('group_statistics/+/*_s6w_',GLM_name,'_',con_names{j});
%         [folder, file, ext] = fileparts(subj_files{i});
%         if j == 1, stat_suff = strcat('FADE_GLM_1a_novelty_',file); end;
%         if j == 2, stat_suff = strcat('FADE_GLM_1a_memory_',file);  end;
%         fact_des_ttest1(subj_files{i}, MS_file, img_path, cov_names, stat_suff);
%     end;
% end;
% spm_exec_mult_jobs;


%%% MS FADE 04: group statistics, 20/07/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_07_20.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_04.mat';
% contrast_meta_batch(subj_file, MS_file, 's6w', 'GLM_TD_1th-a', [1  0 -1], 'con_0001.nii');
% contrast_meta_batch(subj_file, MS_file, 's6w', 'GLM_TD_1th-a', [0  1  0], 'con_0002.nii');
% 
% subj_files = {'subjects_all_2020_04_10_JS.xls';
%               'subjects_all_2020_04_10_JS_young.xls';
%               'subjects_all_2020_04_10_JS_old.xls'};
% subj_suffs = {'', '_young', '_old'};
% 
% MS_file   = 'FADE_model_spaces/MS_FADE_04.mat';
% GLM_name  = 'GLM_TD_1th-a';
% con_name  = 'con_0001.nii';
% img_path  = strcat('group_statistics/+/*_s6w_',GLM_name,'_',con_name);
% cov_names = {};
% for i = 1:numel(subj_files)
%     stat_suff = strcat(GLM_name,'_novelty-master',subj_suffs{i});
%     fact_des_ttest1(subj_files{i}, MS_file, img_path, cov_names, stat_suff);
% end;
% spm_exec_mult_jobs;


%%% one-sample t-test across studies, 20/07/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_files{1} = 'subjects_all_2020_07_20.xls';
% subj_files{2} = 'C:\Users\Joram\ownCloud\DZNE\yFADE\tools\subjects_2020_07_20.xls';
% MS_file       = 'FADE_model_spaces/MS_FADE_04.mat';
% img_paths(1,:)= {'group_statistics/+/*_s6w_GLM_TD_1th-a_con_0001.nii', ...
%                  'group_statistics/+/*_s6_3w_GLM_TD_1th-a_con_0001.nii'};
% img_paths(2,:)= {'group_statistics/+/*_s6w_GLM_TD_1th-a_con_0002.nii', ...
%                  'group_statistics/+/*_s6_3w_GLM_TD_1th-a_con_0002.nii'};
% cov_names     = {'AiA', 'young', 'Verio', 'male'};
% stat_suff     = {'GLM_1a_novelty-master_study_age_scanner_sex_ttest', 'GLM_1a_memory_study_age_scanner_sex_ttest'};
% for i = 1:numel(stat_suff)
%     fact_des_ttest1_study(subj_files, MS_file, img_paths(i,:), cov_names, stat_suff{i});
% end;
% spm_exec_mult_jobs;


%%% full factorial across studies, 20/07/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_files{1} = 'subjects_all_2020_07_20.xls';
% subj_files{2} = 'C:\Users\Joram\ownCloud\DZNE\yFADE\tools\subjects_2020_07_20.xls';
% MS_file       = 'FADE_model_spaces/MS_FADE_04.mat';
% img_paths(1,:)= {'group_statistics/+/*_s6w_GLM_TD_1th-a_con_0001.nii', ...
%                  'group_statistics/+/*_s6_3w_GLM_TD_1th-a_con_0001.nii'};
% img_paths(2,:)= {'group_statistics/+/*_s6w_GLM_TD_1th-a_con_0002.nii', ...
%                  'group_statistics/+/*_s6_3w_GLM_TD_1th-a_con_0002.nii'};
% var_names     = {'scanner', 'sex'};
% cov_names     = {'AiA', 'age'};
% stat_suff     = {'GLM_1a_novelty-master_study_age_scanner_sex_ANOVA', 'GLM_1a_memory_study_age_scanner_sex_ANOVA'};
% for i = 1:numel(stat_suff)
%     fact_des_full_fact_study(subj_files, MS_file, img_paths(i,:), var_names, cov_names, stat_suff{i});
% end;
% spm_exec_mult_jobs;


%%% full factorial across studies, 15/07/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_files{1} = 'subjects_all_2020_07_15.xls';
% subj_files{2} = 'C:\Users\Joram\ownCloud\DZNE\yFADE\tools\subjects_2020_07_15.xls';
% MS_file       = 'FADE_model_spaces/MS_FADE_04.mat';
% img_paths(1,:)= {'group_statistics/+/*_s6w_GLM_TD_1th-a_beta_0001.nii', ...
%                  'group_statistics/+/*_s6_3w_GLM_TD_1th-a_beta_0001.nii'};
% img_paths(2,:)= {'group_statistics/+/*_s6w_GLM_TD_1th-a_beta_0002.nii', ...
%                  'group_statistics/+/*_s6_3w_GLM_TD_1th-a_beta_0002.nii'};
% img_paths(3,:)= {'group_statistics/+/*_s6w_GLM_TD_2nf_con_0002.nii', ...
%                  'group_statistics/+/*_s6_3w_GLM_TD_2nf_con_0002.nii'};
% img_paths(4,:)= {'group_statistics/+/*_s6w_GLM_TD_1dd-ip_beta_0002.nii', ...
%                  'group_statistics/+/*_s6_3w_GLM_TD_1dd-ip_beta_0002.nii'};
% var_names     = {'AiA_yFADE'};
% cov_names     = {'sex'};
% stat_suff     = {'GLM_1a_novelty_AiA_vs_yFADE', 'GLM_1a_memory_AiA_vs_yFADE', 'GLM_2nf_memory_AiA_vs_yFADE', 'GLM_1ip_memory_AiA_vs_yFADE'};
% for i = 1:numel(stat_suff)
%     fact_des_full_fact_study(subj_files, MS_file, img_paths(i,:), var_names, cov_names, stat_suff{i});
% end;
% spm_exec_mult_jobs;


%%% calculate average delays, 16/06/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file  = 'subjects_all_2020_04_10_JS.xls';
% delay_file = 'FADE_delay_249_subjects_2020_04_08_AR_JS_AR.xls';
% delay_file = 'FADE_delay_256_subjects_2020_06_16_JS.xls';
% analyze_delays(subj_file, delay_file);


%%% MS FADE 04: group statistics (245 subj), 16/06/2020 %%%%%%%%%%%%%%%%%%%

% subj_files = {'subjects_all_2020_04_10_JS.xls';
%               'subjects_all_2020_04_10_JS_young.xls';
%               'subjects_all_2020_04_10_JS_old.xls'};
% subj_suffs = {'', '_young', '_old'};
% 
% MS_file   = 'FADE_model_spaces/MS_FADE_04.mat';
% GLM_names ={'GLM_TD_1th-a', 'GLM_TD_1th-a', 'GLM_TD_2nf', 'GLM_TD_1dd-ip'};
% con_names ={'beta_0001.nii', 'beta_0002.nii', 'con_0002.nii', 'beta_0002.nii'};
% cov_names ={};
% for i = 1:numel(subj_files)
%     for j = 1:numel(GLM_names)
%         img_path = strcat('group_statistics/+/*_s6w_',GLM_names{j},'_',con_names{j});
%         if j == 1, stat_suff = strcat(GLM_names{j},'_novelty',subj_suffs{i}); end;
%         if j >  1, stat_suff = strcat(GLM_names{j},'_memory', subj_suffs{i});  end;
%         fact_des_ttest1(subj_files{i}, MS_file, img_path, cov_names, stat_suff);
%     end;
% end;
% spm_exec_mult_jobs;


%%% MS FADE 04: contrast_meta_batch, 16/06/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_04_10_JS.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_04.mat';
% contrast_meta_batch(subj_file, MS_file, 's6w', 'GLM_TD_2nf', [1  1 -2], 'con_0001.nii');
% contrast_meta_batch(subj_file, MS_file, 's6w', 'GLM_TD_2nf', [1 -1  0], 'con_0002.nii');


%%% MS FADE 04: further analyses (245 subj), 26/05/2020 %%%%%%%%%%%%%%%%%%%

% subj_files = {'subjects_all_2020_04_10_JS.xls';
%               'subjects_all_2020_04_10_JS_young.xls';
%               'subjects_all_2020_04_10_JS_old.xls'};
% subj_suffs = {'', '_young', '_old'};
% 
% % 10a, 10b
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_WM3.mat';
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat('TD_WM3', subj_suffs{i})));
%     spm_jobman('run', matlabbatch);
% end;
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_WM2.mat';
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat('TD_WM2', subj_suffs{i})));
%     spm_jobman('run', matlabbatch);
% end;


%%% behavior logfile data, 04/05/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_04_10_JS.xls';
% behavioral_data(subj_file);
% bhvr_logfile_data(subj_file);
% bhvr_memorability(subj_file);


%%% analyze subject demographics, 20/04/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze_subjects('subjects_all_2020_04_10_JS.xls')


%%% MS FADE 04: complete analysis (245 subj), 20/04/2020ff. %%%%%%%%%%%%%%%

% subj_files = {'subjects_all_2020_04_10_JS.xls';
%               'subjects_all_2020_04_10_JS_young.xls';
%               'subjects_all_2020_04_10_JS_old.xls'};
% subj_suffs = {'', '_young', '_old'};

% 0
% MS_file = 'FADE_model_spaces/MS_FADE_04_00_0.mat';
% family.mods = [1 1 1 1 2 2 2 2];
% family.fams = {'GLMs_PE', 'GLMs_TD'};
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat('00_0_PE_vs_TD', subj_suffs{i})));
%     MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
%     load(strcat(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man.dir{1},'BMS.mat'));
%     MS_SMM_BMS_fams(BMS)
% end;

% 1
% MS_file = 'FADE_model_spaces/MS_FADE_04_00_0.mat';
% family.mods = [1 1 2 2 1 1 2 2];
% family.fams = {'GLMs_00', 'GLMs_0'};
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat('00_0_00_vs_0', subj_suffs{i})));
%     MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
%     load(strcat(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man.dir{1},'BMS.mat'));
%     MS_SMM_BMS_fams(BMS)
% end;

% 2
% MS_file = 'FADE_model_spaces/MS_FADE_04_00_0.mat';
% family.mods = [1 2 1 2 1 2 1 2];
% family.fams = {'GLMs_x1', 'GLMs_x2'};
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat('00_0_x1_vs_x2', subj_suffs{i})));
%     MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
%     load(strcat(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man.dir{1},'BMS.mat'));
%     MS_SMM_BMS_fams(BMS)
% end;

% 3
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_0_1_2.mat';
% cvLFE_meta_batch(subj_files{1}, MS_file, 's6w', [1 0 2 2 2 2 2 2 3 3 3], {'GLMs_TD_0', 'GLMs_TD_1', 'GLMs_TD_2'});
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_0_1_2_fams.mat';
% family.mods = [1 2 2];
% family.fams = {'GLMs_TD_0', 'GLMs_TD_12'};
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLFE.nii', strcat('TD_0_1_2_12_vs_0', subj_suffs{i})));
%     MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
%     load(strcat(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man.dir{1},'BMS.mat'));
%     MS_SMM_BMS_fams(BMS)
% end;

% 4a
% subjects_all_2020_04_10_GLMs
% subj_file = 'subjects_all_2020_04_10_GLM_TD_5.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_3_5.mat';
% cvLFE_meta_batch(subj_file, MS_file, 's6w', [1 2], {'GLMs_TD_3', 'GLMs_TD_5'});
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_0-5_fams.mat';
% cvLBF_meta_batch(subj_file, MS_file, 's6w', [2 1 0 0 0], 'MA_cvLFE.nii', 'GLMs_TD_1_vs_0');
% cvLBF_meta_batch(subj_file, MS_file, 's6w', [2 0 1 0 0], 'MA_cvLFE.nii', 'GLMs_TD_2_vs_0');
% cvLBF_meta_batch(subj_file, MS_file, 's6w', [2 0 0 1 0], 'MA_cvLFE.nii', 'GLMs_TD_3_vs_0');
% cvLBF_meta_batch(subj_file, MS_file, 's6w', [2 0 0 0 1], 'MA_cvLFE.nii', 'GLMs_TD_5_vs_0');

% 4b
% subj_files={'subjects_all_2020_04_10_GLM_TD_5.xls';
%             'subjects_all_2020_04_10_GLM_TD_5_young.xls';
%             'subjects_all_2020_04_10_GLM_TD_5_old.xls'};
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_0-5_fams.mat';
% img_paths ={'model_selection/+/*_s6w_GLMs_TD_1_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s6w_GLMs_TD_2_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s6w_GLMs_TD_3_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s6w_GLMs_TD_5_vs_0_MC_cvLBF.nii'};
% var_names ={};
% cov_names ={'sex'};
% stat_suff = 'LBF_GLMs_TD_1235_vs_0';
% for i = 1:numel(subj_files)
%     fact_des_full_fact(subj_files{i}, MS_file, img_paths, var_names, cov_names, strcat(stat_suff, subj_suffs{i}));
% end;
% spm_exec_mult_jobs;

% 5
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_1_2.mat';
% family.mods = [1 1 1 1 1 1 2 2 2];
% family.fams = {'GLMs_TD_1', 'GLMs_TD_2'};
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat('TD_1_2_1_vs_2', subj_suffs{i})));
%     MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
%     load(strcat(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man.dir{1},'BMS.mat'));
%     MS_SMM_BMS_fams(BMS)
% end;

% 6
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_1.mat';
% family.mods = [1 1 1 2 2 2];
% family.fams = {'GLMs_TD_1dd', 'GLMs_TD_1th'};
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat('TD_1_dd_vs_th', subj_suffs{i})));
%     MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
%     load(strcat(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man.dir{1},'BMS.mat'));
%     MS_SMM_BMS_fams(BMS)
% end;

% 7, 8, 9
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_2.mat';
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat('TD_2_nf_nr_ns', subj_suffs{i})));
%     spm_jobman('run', matlabbatch);
% end;
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_1dd.mat';
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat('TD_1dd_ip_cp_lr', subj_suffs{i})));
%     spm_jobman('run', matlabbatch);
% end;
% MS_file = 'FADE_model_spaces/MS_FADE_04_TD_1th.mat';
% for i = 1:numel(subj_files)
%     load(cvBMS_analysis_batch(subj_files{i}, MS_file, 's6w', 'MA_cvLME.nii', strcat('TD_1th_l_a_s', subj_suffs{i})));
%     spm_jobman('run', matlabbatch);
% end;


%%% MS FADE 04: complete analysis (235 subj), 20/04/2020 %%%%%%%%%%%%%%%%%%

% 0, 1, 2
% subj_file = 'subjects_all_2020_04_10_JS_incl.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_04_00_0.mat';
% load(cvBMS_analysis_batch(subj_file, MS_file, 's6w', 'MA_cvLME.nii', '00_0_PE_vs_TD'));
% family.mods = [1 1 1 1 2 2 2 2];
% family.fams = {'GLMs_PE', 'GLMs_TD'};
% MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
% load(cvBMS_analysis_batch(subj_file, MS_file, 's6w', 'MA_cvLME.nii', '00_0_00_vs_0'));
% family.mods = [1 1 2 2 1 1 2 2];
% family.fams = {'GLMs_00', 'GLMs_0'};
% MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
% load(cvBMS_analysis_batch(subj_file, MS_file, 's6w', 'MA_cvLME.nii', '00_0_x1_vs_x2'));
% family.mods = [1 2 1 2 1 2 1 2];
% family.fams = {'GLMs_x1', 'GLMs_x2'};
% MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);

% 3
% subj_file = 'subjects_all_2020_04_10_JS_incl.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_0_1_2.mat';
% cvLFE_meta_batch(subj_file, MS_file, 's6w', [1 0 2 2 2 2 2 2 3 3 3], {'GLMs_TD_0', 'GLMs_TD_1', 'GLMs_TD_2'});
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_0_1_2_fams.mat';
% load(cvBMS_analysis_batch(subj_file, MS_file, 's6w', 'MA_cvLFE.nii', 'TD_0_1_2_12_vs_0'));
% family.mods = [1 2 2];
% family.fams = {'GLMs_TD_0', 'GLMs_TD_12'};
% MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);

% 5, 6
% subj_file = 'subjects_all_2020_04_10_JS_incl.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_1_2.mat';
% load(cvBMS_analysis_batch(subj_file, MS_file, 's6w', 'MA_cvLME.nii', 'TD_1_2_1_vs_2'));
% family.mods = [1 1 1 1 1 1 2 2 2];
% family.fams = {'GLMs_TD_1', 'GLMs_TD_2'};
% MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_1.mat';
% load(cvBMS_analysis_batch(subj_file, MS_file, 's6w', 'MA_cvLME.nii', 'TD_1_dd_vs_th'));
% family.mods = [1 1 1 2 2 2];
% family.fams = {'GLMs_TD_1dd', 'GLMs_TD_1th'};
% MS_BMS_group_fams(matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man, 'RFX-VB', family, true);

% 7, 8, 9
% subj_file = 'subjects_all_2020_04_10_JS_incl.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_2.mat';
% load(cvBMS_analysis_batch(subj_file, MS_file, 's6w', 'MA_cvLME.nii', 'TD_2_nf_nr_ns'));
% spm_jobman('run', matlabbatch);
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_1dd.mat';
% load(cvBMS_analysis_batch(subj_file, MS_file, 's6w', 'MA_cvLME.nii', 'TD_1dd_ip_cp_lr'));
% spm_jobman('run', matlabbatch);
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_1th.mat';
% load(cvBMS_analysis_batch(subj_file, MS_file, 's6w', 'MA_cvLME.nii', 'TD_1th_l_a_s'));
% spm_jobman('run', matlabbatch);

% 4a
% subjects_all_2020_04_10_GLMs
% subj_file = 'subjects_all_2020_04_10_GLM_TD_5.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_3_5.mat';
% cvLFE_meta_batch(subj_file, MS_file, 's6w', [1 2], {'GLMs_TD_3', 'GLMs_TD_5'});
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_0-5_fams.mat';
% cvLBF_meta_batch(subj_file, MS_file, 's6w', [2 1 0 0 0], 'MA_cvLFE.nii', 'GLMs_TD_1_vs_0');
% cvLBF_meta_batch(subj_file, MS_file, 's6w', [2 0 1 0 0], 'MA_cvLFE.nii', 'GLMs_TD_2_vs_0');
% cvLBF_meta_batch(subj_file, MS_file, 's6w', [2 0 0 1 0], 'MA_cvLFE.nii', 'GLMs_TD_3_vs_0');
% cvLBF_meta_batch(subj_file, MS_file, 's6w', [2 0 0 0 1], 'MA_cvLFE.nii', 'GLMs_TD_5_vs_0');

% 4b
% subj_file = 'subjects_all_2020_04_10_GLM_TD_5.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_04_TD_0-5_fams.mat';
% img_paths ={'model_selection/+/*_s6w_GLMs_TD_1_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s6w_GLMs_TD_2_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s6w_GLMs_TD_3_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s6w_GLMs_TD_5_vs_0_MC_cvLBF.nii'};
% var_names ={};
% cov_names ={'age_group0', 'sex'};
% stat_suff = 'LBF_GLMs_TD_1235_vs_0';
% fact_des_full_fact(subj_file, MS_file, img_paths, var_names, cov_names, stat_suff);
% spm_exec_mult_jobs;


%%% rest of subjects, 14/04/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stats_meta_batch('subjects_rest.xls', 'FADE_model_spaces/MS_FADE_04.mat', 's6w', 'prepare')
% stats_meta_batch('subjects_rest.xls', 'FADE_model_spaces/MS_FADE_04.mat', 's6w', 'perform')
% cvLME_meta_batch('subjects_rest.xls', 'FADE_model_spaces/MS_FADE_04.mat', 's6w')
% collect_meta_batch('subjects_rest.xls', 'FADE_model_spaces/MS_FADE_04.mat', 's6w', 'MA_cvLME.nii');
% 
% subj_file = 'subjects_rest.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_04_3regs.mat';
% collect_meta_batch(subj_file, MS_file, 's6w', 'beta_0001.nii');
% collect_meta_batch(subj_file, MS_file, 's6w', 'beta_0002.nii');
% collect_meta_batch(subj_file, MS_file, 's6w', 'beta_0003.nii');
% collect_meta_batch(subj_file, MS_file, 's6w', 'mask.nii');
% collect_meta_batch(subj_file, MS_file, 's6w', 'ResMS.nii');
% collect_meta_batch(subj_file, MS_file, 's6w', 'SPM.mat');


%%% collect meta batch, 14/04/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collect_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_04.mat', 's6w', 'MA_cvLME.nii');

% subj_file = 'subjects_test.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_04_3regs.mat';
% collect_meta_batch(subj_file, MS_file, 's6w', 'beta_0001.nii');
% collect_meta_batch(subj_file, MS_file, 's6w', 'beta_0002.nii');
% collect_meta_batch(subj_file, MS_file, 's6w', 'beta_0003.nii');
% collect_meta_batch(subj_file, MS_file, 's6w', 'mask.nii');
% collect_meta_batch(subj_file, MS_file, 's6w', 'ResMS.nii');
% collect_meta_batch(subj_file, MS_file, 's6w', 'SPM.mat');


%%% cvLME meta batch, 14/04/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cvLME_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_04.mat', 's6w')


%%% check for logfiles, 11/04/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check_for_logfiles('subjects_all_2020_01_30.txt');
% check_for_logfiles('subjects_all_2020_04_10.txt');
% check_for_logfiles('subjects_excluded_AR_JS.xls');


%%% stats meta batch, 08/04/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stats_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_04.mat', 's6w', 'prepare')
% stats_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_04.mat', 's6w', 'perform')


%%% preproc meta batch, 04-05/04/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preproc_meta_batch('subjects_test.xls', 'prepare')
% preproc_meta_batch('subjects_test.xls', 'perform')

% preproc_meta_batch('subjects_rest.xls', 'prepare')
% preproc_meta_batch('subjects_rest.xls', 'perform')


%%% calculate GM volumes, 31/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ../project_directories_BS_JS.mat
% subj_file = 'subjects_all_2020_01_30.txt';
% MS_name   = 's0_GM';
% GM_img    = 'mwp1*_T1-MPRAGE_1.nii';
% ROI_imgs  ={'sHippo_CA_Sub.nii';
%             'sHippo_CA_Sub_left.nii';
%             'sHippo_CA_Sub_right.nii';
%             'sHippocampus_AAL.nii';
%             'sHippocampus_AAL_left.img';
%             'sHippocampus_AAL_right.img'};
% for j = 1:numel(ROI_imgs)
%     ROI_imgs{j} = strcat(tool_dir,'ROI_images/',ROI_imgs{j});
% end;
% calc_GM_group(subj_file, MS_name, GM_img, ROI_imgs, 'hippocampus')


%%% skull-strip meta batch, 27/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_test.xls';
% tp_thr    = 0.95;
% dest_dir  = 'C:\Joram\projects\DZNE\FADE\analyses_BS\structural_analyses\T1_images\';
% skull_strip_meta_batch(subj_file, tp_thr, dest_dir);


%%% logistic regression (X2y), 26/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load C:\Joram\projects\DZNE\FADE\analyses_BS\model_selection\MS_FADE_02_LBF_GLMs_1x1_vs_0_vs_age\SPM.mat
% fact_des_log_reg_X2y(SPM, 'lodds', true, 2);
% load C:\Joram\projects\DZNE\FADE\analyses_BS\model_selection\MS_FADE_02_LBF_GLMs_1x1_vs_2_vs_age\SPM.mat
% fact_des_log_reg_X2y(SPM, 'lodds', true, 2);


%%% logistic regression (y2X), 26/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load C:\Joram\projects\DZNE\FADE\analyses_BS\group_statistics\MS_FADE_02_GLM_1x1a_memory_young_vs_old\SPM.mat
% fact_des_log_reg_y2X(SPM, true, 2);
% load C:\Joram\projects\DZNE\FADE\analyses_BS\group_statistics\MS_FADE_02_GLM_1x1a_novel-master_young_vs_old\SPM.mat
% fact_des_log_reg_y2X(SPM, true, 2);


%%% behavioral data analysis, 25/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% behavioral_data(subj_file);


%%% calculate SAMe scores, 25/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ../project_directories_BS_JS.mat
% subj_file = 'subjects_all_2020_01_30.txt';
% MS_name   = 'MS_FADE_02';
% GLM_names ={'GLM_1x1a', 'GLM_1x1l', 'GLM_1x2nf', 'GLM_1x2nr'};
% img_pref  = 's8';
% con_vecs  ={[0 1 0], [0 1 0], [+1 -1 0], [+1 -1 0]};
% ref_maps  = cell(size(con_vecs));
% for j = 1:numel(GLM_names)
%     ref_maps{j} = strcat(stat_dir,'group_statistics/',MS_name,'_',GLM_names{j},...
%                          '_memory_young_no_gender/','con_0001_FWE_0.001_10.nii');
% end;
% calc_SAMe_group(subj_file, MS_name, GLM_names, img_pref, con_vecs, ref_maps, 'memory');


%%% one-sample t-tests, 25/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30_young.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% GLM_names ={'GLM_1x1l', 'GLM_1x1a', 'GLM_1x2nf', 'GLM_1x2nr'};
% con_names ={'beta_0002.nii', 'beta_0002.nii', 'con_0002.nii', 'con_0002.nii'};
% cov_names ={};
% for i = 1:numel(GLM_names)
%     img_path  = strcat('group_statistics/+/*_s8_',GLM_names{i},'_',con_names{i});
%     stat_suff = strcat(GLM_names{i},'_memory_young_no_gender');
%     fact_des_ttest1(subj_file, MS_file, img_path, cov_names, stat_suff);
% end;
% spm_exec_mult_jobs;


%%% stats meta batch, 24/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stats_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_03.mat', 's8', 'prepare')
% stats_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_03.mat', 's8', 'perform')

% stats_meta_batch('subjects_rest.xls', 'FADE_model_spaces/MS_FADE_03.mat', 's8', 'prepare')
% stats_meta_batch('subjects_rest.xls', 'FADE_model_spaces/MS_FADE_03.mat', 's8', 'perform')


%%% model space test, 23/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% memory = repmat([1, 2, 4, 5]', [25 1]);    
% num12  = sum(memory==1 | memory==2);
% num45  = sum(memory==4 | memory==5);
% ind3   = find(memory==3);
% ind3i  = randperm(numel(ind3));
% ind3a  = ind3i(1:ceil((num12/(num12+num45))*numel(ind3)));
% ind3b  = ind3i(ceil((num12/(num12+num45))*numel(ind3))+1:end);
% ind123 = sort([find(ismember(memory,[1 2])); ind3a]);
% ind345 = sort([find(ismember(memory,[4 5])); ind3b]);


%%% behavioral data analysis, 18/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% behavioral_data(subj_file);
% bhvr_memorability(subj_file);
% bhvr_logfile_data(subj_file);


%%% behavior/logfile data analysis, 16/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% bhvr_logfile_data(subj_file);


%%% multiple regression, 09/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30_FADE_old.xls';
% MS_file   = 'FADE_model_spaces/MS_VBM_s8_GM.mat';
% img_path  = 'structural_analyses/+/smwp1*_T1-MPRAGE_1.nii';
% mean_cent = [true, true, false, false];
% var_names ={'FADE_GLM_1x1a_pos_FWE_0.001_10', 'age', 'gender', 'scanner'};
% stat_suff = 'old_vs_FADE_1x1a_pos_0.001';
% fact_des_mult_reg(subj_file, MS_file, img_path, var_names, mean_cent, stat_suff)
% var_names ={'FADE_GLM_1x1a_neg_FWE_0.001_10', 'age', 'gender', 'scanner'};
% stat_suff = 'old_vs_FADE_1x1a_neg_0.001';
% fact_des_mult_reg(subj_file, MS_file, img_path, var_names, mean_cent, stat_suff)
% spm_exec_mult_jobs;


%%% full factorial analysis, 27/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_VBM_s8_GM.mat';
% img_paths ={'structural_analyses/+/smwp1*_T1-MPRAGE_1.nii'};
% var_names ={'scanner', 'gender', 'age_group'};
% cov_names ={};
% stat_suff = 'by_scanner_gender_age';
% fact_des_full_fact(subj_file, MS_file, img_paths, var_names, cov_names, stat_suff)


%%% two-sample t-test, 09/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_VBM_s8_GM.mat';
% img_path  = 'structural_analyses/+/smwp1*_T1-MPRAGE_1.nii';
% var_name  = 'age_group';
% cov_names ={'gender'};
% stat_suff = 'young_vs_old';
% fact_des_ttest2(subj_file, MS_file, img_path, var_name, cov_names, stat_suff)


%%% behavior/memorability analysis, 05/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% bhvr_memorability(subj_file);


%%% calculate FADE scores, 03/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ../project_directories_BS_JS.mat
% subj_file = 'subjects_all_2020_01_30.txt';
% MS_name   = 'MS_FADE_02';
% GLM_names ={'GLM_1x1a', 'GLM_1x1a'};
% img_pref  = 's8';
% con_vecs  ={[0 +1 0], [0 +1 0];
%             [0 -1 0], [0 -1 0]};
% con_types =['tt'; 'tt'];
% ref_maps  = cell(size(con_vecs));
% for j = 1:numel(GLM_names)
%     for k = 1:size(con_vecs,1)
%         if j == 1
%             ref_maps{k,j} = strcat(stat_dir,'group_statistics/',MS_name,'_',GLM_names{j},'_memory_young/',...
%                                    'con_000',num2str(k+1),'_FWE_0.001_10.nii');
%         else
%             ref_maps{k,j} = strcat(stat_dir,'group_statistics/',MS_name,'_',GLM_names{j},'_memory_young/',...
%                                    'con_000',num2str(k+1),'_FWE_0.05_10.nii');
%         end;
%     end;
% end;
% calc_FADE_group(subj_file, MS_name, GLM_names, img_pref, con_vecs, con_types, ref_maps, 'memory');


%%% behavioral data analysis, 20/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% behavioral_data('subjects_all_2020_01_30.txt');
% behavioral_data('subjects_all_2020_02_17.txt');


%%% second-level cvBMS, 02/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MS_file   = 'FADE_model_spaces/MS_FADE_02_only_1x2nf2nr.mat';
% subj_file = 'subjects_all_2020_01_30.txt';
% cvBMS_analysis_batch(subj_file, MS_file, 's8', 'MA_cvLME.nii', 'only_1x2nf2nr')
% subj_file = 'subjects_all_2020_01_30_young.txt';
% cvBMS_analysis_batch(subj_file, MS_file, 's8', 'MA_cvLME.nii', 'only_1x2nf2nr_young')
% subj_file = 'subjects_all_2020_01_30_old.txt';
% cvBMS_analysis_batch(subj_file, MS_file, 's8', 'MA_cvLME.nii', 'only_1x2nf2nr_old')
% spm_exec_mult_jobs;


%%% second-level cvBMS, 02/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MS_file   = 'FADE_model_spaces/MS_FADE_02_fams_1x1_2.mat';
% subj_file = 'subjects_all_2020_01_30.txt';
% cvBMS_analysis_batch(subj_file, MS_file, 's8', 'MA_cvLFE.nii', 'fams_1x1_2')
% subj_file = 'subjects_all_2020_01_30_young.txt';
% cvBMS_analysis_batch(subj_file, MS_file, 's8', 'MA_cvLFE.nii', 'fams_1x1_2_young')
% subj_file = 'subjects_all_2020_01_30_old.txt';
% cvBMS_analysis_batch(subj_file, MS_file, 's8', 'MA_cvLFE.nii', 'fams_1x1_2_old')
% spm_exec_mult_jobs;


%%% second-level cvBMS, 02/03/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MS_file   = 'FADE_model_spaces/MS_FADE_02_fams_1x0_12.mat';
% subj_file = 'subjects_all_2020_01_30.txt';
% cvBMS_analysis_batch(subj_file, MS_file, 's8', 'MA_cvLFE.nii', 'fams_1x0_12')
% subj_file = 'subjects_all_2020_01_30_young.txt';
% cvBMS_analysis_batch(subj_file, MS_file, 's8', 'MA_cvLFE.nii', 'fams_1x0_12_young')
% subj_file = 'subjects_all_2020_01_30_old.txt';
% cvBMS_analysis_batch(subj_file, MS_file, 's8', 'MA_cvLFE.nii', 'fams_1x0_12_old')
% spm_exec_mult_jobs;


%%% calculate FADE scores, 24/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ../project_directories_BS_JS.mat
% subj_file = 'subjects_all_2020_01_30.txt';
% MS_name   = 'MS_FADE_02';
% GLM_names ={'GLM_1x1a', 'GLM_1x1l', 'GLM_1x2nf', 'GLM_1x2nr'};
% img_pref  = 's8';
% con_vecs  ={[0  1 0], [0  1 0], [ 1 -1 0], [ 1 -1 0];
%             [0 +1 0], [0 +1 0], [+1 -1 0], [+1 -1 0];
%             [0 -1 0], [0 -1 0], [-1 +1 0], [-1 +1 0]};
% con_types =['FFFF'; 'tttt'; 'tttt'];
% ref_maps  = cell(size(con_vecs));
% for j = 1:numel(GLM_names)
%     for k = 1:size(con_vecs,1)
%         ref_maps{k,j} = strcat(stat_dir,'group_statistics/',MS_name,'_',GLM_names{j},'_memory_young/',...
%                                'con_000',num2str(k),'_FWE_0.05_10.nii');
%     end;
% end;
% calc_FADE_group(subj_file, MS_name, GLM_names, img_pref, con_vecs, con_types, ref_maps, 'memory');


%%% full factorial analysis, 27/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% img_paths ={'group_statistics/+/*_s8_GLM_1x1a_con_0001.nii'};
% var_names ={'scanner', 'gender', 'age_group'};
% cov_names ={};
% stat_suff = 'GLM_1x1a_novel-master_by_scanner_gender_age';
% fact_des_full_fact(subj_file, MS_file, img_paths, var_names, cov_names, stat_suff)
% img_paths ={'group_statistics/+/*_s8_GLM_1x1a_beta_0002.nii'};
% stat_suff = 'GLM_1x1a_memory_by_scanner_gender_age';
% fact_des_full_fact(subj_file, MS_file, img_paths, var_names, cov_names, stat_suff)
% spm_exec_mult_jobs;


%%% calculate FADE scores, 24/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_JS_sa.txt'; % 'subjects_JS.xls';
% MS_name   = 'MS_FADE_02';       GLM_name  = 'GLM_1x1a';
% img_pref  = 's8';               con_vec   = [0 1 0];
% ref_map   = 'C:\Joram\projects\DZNE\FADE\analyses_BS\group_statistics\MS_FADE_02_GLM_1x1a_memory_young\con_0001_FWE_0.05_10.nii';
% calc_FADE_group(subj_file, MS_name, GLM_name, img_pref, con_vec, ref_map);


%%% one-sample t-tests, 24/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30_young.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% GLM_names ={'GLM_1x1l', 'GLM_1x1a', 'GLM_1x2nf', 'GLM_1x2nr'};
% con_names ={'beta_0002.nii', 'beta_0002.nii', 'con_0002.nii', 'con_0002.nii'};
% cov_names ={'gender'};
% for i = 1:numel(GLM_names)
%     img_path  = strcat('group_statistics/+/*_s8_',GLM_names{i},'_',con_names{i});
%     stat_suff = strcat(GLM_names{i},'_memory_young');
%     fact_des_ttest1(subj_file, MS_file, img_path, cov_names, stat_suff);
% end;
% spm_exec_mult_jobs;


%%% full factorial analysis, 17/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30_GLM_1x5.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% img_paths ={'model_selection/+/*_s8_GLMs_1x1_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s8_GLMs_1x2_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s8_GLM_1x3_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s8_GLM_1x5_vs_0_MC_cvLBF.nii'};
% var_names ={'age_group'};
% cov_names ={'gender'};
% stat_suff = 'LBF_GLMs_1x1235_vs_0_young_vs_old';
% fact_des_full_fact(subj_file, MS_file, img_paths, var_names, cov_names, stat_suff)
% img_paths ={'model_selection/+/*_s8_GLM_1x1l_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s8_GLM_1x1a_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s8_GLM_1x2nf_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s8_GLM_1x2nr_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s8_GLM_1x3_vs_0_MC_cvLBF.nii';
%             'model_selection/+/*_s8_GLM_1x5_vs_0_MC_cvLBF.nii'};
% stat_suff = 'LBF_GLMs_1x1la2nfnr35_vs_0_young_vs_old';
% fact_des_full_fact(subj_file, MS_file, img_paths, var_names, cov_names, stat_suff)
% spm_exec_mult_jobs;


%%% log Bayes factors, 21/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MS_file   = 'FADE_model_spaces/MS_FADE_02_only_1x.mat';
% subj_file = 'subjects_all_2020_01_30_GLM_1x3.xls';
% cvLBF_meta_batch(subj_file, MS_file, 's8', [0 2 0 0 0 0 1 0], 'MA_cvLME.nii', 'GLM_1x3_vs_0');
% subj_file = 'subjects_all_2020_01_30_GLM_1x5.xls';
% cvLBF_meta_batch(subj_file, MS_file, 's8', [0 2 0 0 0 0 0 1], 'MA_cvLME.nii', 'GLM_1x5_vs_0');


%%% paired t-test, 17/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% img_paths ={'group_statistics/+/*_s8_GLM_1x2nf_con_0002.nii';
%             'group_statistics/+/*_s8_GLM_1x2nr_con_0002.nii'};
% cov_names ={'gender'};
% stat_suff = 'GLMs_1x2nf_vs_2nf_memory';
% fact_des_ttestp(subj_file, MS_file, img_paths, cov_names, stat_suff)


%%% full factorial analysis, 17/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% img_paths ={'group_statistics/+/*_s8_GLM_1x2nf_con_0002.nii';
%             'group_statistics/+/*_s8_GLM_1x2nr_con_0002.nii'};
% var_names ={'age_group'};
% cov_names ={'gender'};
% stat_suff = 'GLMs_1x2nf_vs_2nr_memory_young_vs_old';
% fact_des_full_fact(subj_file, MS_file, img_paths, var_names, cov_names, stat_suff)


%%% full factorial analysis, 20/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% img_paths ={'group_statistics/+/*_s8_GLM_1x1l_con_0002.nii',  'group_statistics/+/*_s8_GLM_1x1a_con_0002.nii';
%             'group_statistics/+/*_s8_GLM_1x2nf_con_0002.nii', 'group_statistics/+/*_s8_GLM_1x2nr_con_0002.nii'};
% var_names ={'age_group'};
% cov_names ={'gender'};
% stat_suff = 'GLMs_1x_vs_2x_memory_adjusted_young_vs_old';
% fact_des_full_fact(subj_file, MS_file, img_paths, var_names, cov_names, stat_suff)


%%% contrast meta batch, 20/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% contrast_meta_batch(subj_file, MS_file, 's8', 'GLM_1x1a', [0 2 0], 'con_0002.nii')
% contrast_meta_batch(subj_file, MS_file, 's8', 'GLM_1x1l', [0 2 0], 'con_0002.nii')


%%% full factorial analysis, 20/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% img_paths ={'group_statistics/+/*_s8_GLM_1x1l_beta_0002.nii', 'group_statistics/+/*_s8_GLM_1x1a_beta_0002.nii';
%             'group_statistics/+/*_s8_GLM_1x2nf_con_0002.nii', 'group_statistics/+/*_s8_GLM_1x2nr_con_0002.nii'};
% var_names ={'age_group'};
% cov_names ={'gender'};
% stat_suff = 'GLMs_1x_vs_2x_memory_young_vs_old';
% fact_des_full_fact(subj_file, MS_file, img_paths, var_names, cov_names, stat_suff)


%%% full factorial analysis, 20/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% GLM_names ={'GLM_1x2nf', 'GLM_1x2nr'};
% for i = 1:numel(GLM_names)
%     % full factorial
%     img_paths ={strcat('group_statistics/+/*_s8_',GLM_names{i},'_beta_0001.nii');
%                 strcat('group_statistics/+/*_s8_',GLM_names{i},'_beta_0002.nii');
%                 strcat('group_statistics/+/*_s8_',GLM_names{i},'_beta_0003.nii')};
%     stat_suff = strcat(GLM_names{i},'_rem_for_master_young_vs_old');
%     fact_des_full_fact(subj_file, MS_file, img_paths, {'age_group'}, {'gender'}, stat_suff)
% end;
% spm_exec_mult_jobs;


%%% contrast meta batch, 20/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% contrast_meta_batch(subj_file, MS_file, 's8', 'GLM_1x2nf', [1 1 -2], 'con_0001.nii')
% contrast_meta_batch(subj_file, MS_file, 's8', 'GLM_1x2nf', [1 -1 0], 'con_0002.nii')
% contrast_meta_batch(subj_file, MS_file, 's8', 'GLM_1x2nr', [1 1 -2], 'con_0001.nii')
% contrast_meta_batch(subj_file, MS_file, 's8', 'GLM_1x2nr', [1 -1 0], 'con_0002.nii')


%%% behavioral data analysis, 20/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% behavioral_data('subjects_all_2020_02_17.txt');


%%% paired t-test, 17/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% img_paths ={'group_statistics/+/*_s8_GLM_1x1l_beta_0002.nii';
%             'group_statistics/+/*_s8_GLM_1x1a_beta_0002.nii'};
% cov_names ={'gender'};
% stat_suff = 'GLMs_1x1l_vs_1a_memory';
% fact_des_ttestp(subj_file, MS_file, img_paths, cov_names, stat_suff)


%%% full factorial analysis, 17/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% img_paths ={'group_statistics/+/*_s8_GLM_1x1l_beta_0002.nii';
%             'group_statistics/+/*_s8_GLM_1x1a_beta_0002.nii'};
% var_names ={'age_group'};
% cov_names ={'gender'};
% stat_suff = 'GLMs_1x1l_vs_1a_memory_young_vs_old';
% fact_des_full_fact(subj_file, MS_file, img_paths, var_names, cov_names, stat_suff)


%%% all group statistics, 17/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% GLM_names ={'GLM_1x1l', 'GLM_1x1a'};
% for i = 1:numel(GLM_names)
%     % full factorial
%     img_paths ={strcat('group_statistics/+/*_s8_',GLM_names{i},'_beta_0001.nii');
%                 strcat('group_statistics/+/*_s8_',GLM_names{i},'_beta_0003.nii')};
%     stat_suff = strcat(GLM_names{i},'_novel_vs_master_young_vs_old');
%     fact_des_full_fact(subj_file, MS_file, img_paths, {'age_group'}, {'gender'}, stat_suff)
%     % two-sample t-tests
%     img_path  = strcat('group_statistics/+/*_s8_',GLM_names{i},'_con_0001.nii');
%     stat_suff = strcat(GLM_names{i},'_novel-master_young_vs_old');
%     fact_des_ttest2(subj_file, MS_file, img_path, 'age_group', {'gender'}, stat_suff)
%     img_path  = strcat('group_statistics/+/*_s8_',GLM_names{i},'_beta_0002.nii');
%     stat_suff = strcat(GLM_names{i},'_memory_young_vs_old');
%     fact_des_ttest2(subj_file, MS_file, img_path, 'age_group', {'gender'}, stat_suff)
% end;
% spm_exec_mult_jobs;


%%% two-sample t-test, 17/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% img_path  = 'group_statistics/+/*_s8_GLM_1x1l_con_0001.nii';
% var_name  = 'age_group';
% cov_names ={'gender'};
% stat_suff = 'GLM_1x1a_novel-master_young_vs_old';
% fact_des_ttest2(subj_file, MS_file, img_path, var_name, cov_names, stat_suff)


%%% contrast meta batch, 17/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% contrast_meta_batch(subj_file, MS_file, 's8', 'GLM_1x1a', [1 0 -1], 'con_0001.nii')
% contrast_meta_batch(subj_file, MS_file, 's8', 'GLM_1x1l', [1 0 -1], 'con_0001.nii')


%%% two-sample t-test, 17/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% img_path  = 'group_statistics/+/*_s8_GLM_1x1l_beta_0002.nii';
% var_name  = 'age_group';
% cov_names ={'gender'};
% stat_suff = 'GLM_1x1a_memory_young_vs_old';
% fact_des_ttest2(subj_file, MS_file, img_path, var_name, cov_names, stat_suff)


%%% full factorial analysis, 17/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% img_paths ={'group_statistics/+/*_s8_GLM_1x1l_beta_0001.nii';
%             'group_statistics/+/*_s8_GLM_1x1l_beta_0003.nii'};
% var_names ={'age_group'};
% cov_names ={'gender'};
% stat_suff = 'GLM_1x1a_novel_vs_master_young_vs_old';
% fact_des_full_fact(subj_file, MS_file, img_paths, var_names, cov_names, stat_suff)


%%% cvLBF regression analysis, 12/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30_memory.xls';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% var_names ={'Aprime'};
% LBF_names ={'GLMs_1x1_vs_0', 'GLMs_1x2_vs_0', 'GLMs_1x1_vs_2', 'GLMs_1x12_vs_0'};
% for i = 1:numel(LBF_names)
%     img_path  = strcat('model_selection/+/*_s8_',LBF_names{i},'_MC_cvLBF.nii');
%     stat_suff = strcat('LBF_',LBF_names{i},'_vs_Aprime');
%     fact_des_mult_reg(subj_file, MS_file, img_path, var_names, stat_suff);
% end;
% spm_exec_mult_jobs;


%%% behavioral data analysis, 10/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% behavioral_data('subjects_all_2020_01_30.txt');


%%% cvLBF regression analysis, 10/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% 
% % create cvLBF images
% cvLFE_meta_batch(subj_file, 'FADE_model_spaces/MS_FADE_02_no_x35.mat', 's8', ...
%                  [0 1 2 2 3 3, 0 0 0 0 0 0], {'GLMs_1x0', 'GLMs_1x1', 'GLMs_1x2'});
% cvLFE_meta_batch(subj_file, 'FADE_model_spaces/MS_FADE_02_no_x35.mat', 's8', ...
%                  [0 0 1 1 1 1, 0 0 0 0 0 0], {'GLMs_1x12'});
% cvLBF_meta_batch(subj_file, 'FADE_model_spaces/MS_FADE_02_only_1x012.mat', 's8', ...
%                  [2 1 0 0], 'MA_cvLFE.nii', 'GLMs_1x1_vs_0');
% cvLBF_meta_batch(subj_file, 'FADE_model_spaces/MS_FADE_02_only_1x012.mat', 's8', ...
%                  [2 0 1 0], 'MA_cvLFE.nii', 'GLMs_1x2_vs_0');
% cvLBF_meta_batch(subj_file, 'FADE_model_spaces/MS_FADE_02_only_1x012.mat', 's8', ...
%                  [2 0 0 1], 'MA_cvLFE.nii', 'GLMs_1x12_vs_0');
% cvLBF_meta_batch(subj_file, 'FADE_model_spaces/MS_FADE_02_only_1x012.mat', 's8', ...
%                  [0 1 2 0], 'MA_cvLFE.nii', 'GLMs_1x1_vs_2');
% 
% % perform multiple regression
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% var_names ={'age'};
% LBF_names ={'GLM_1x1l_vs_0', 'GLM_1x1a_vs_0', 'GLM_1x2nf_vs_0', 'GLM_1x2nr_vs_0', 'GLMs_1x1_vs_0', 'GLMs_1x2_vs_0', 'GLMs_1x1_vs_2', 'GLMs_1x12_vs_0'};
% for i = 1:numel(LBF_names)
%     img_path  = strcat('model_selection/+/*_s8_',LBF_names{i},'_MC_cvLBF.nii');
%     stat_suff = strcat('LBF_',LBF_names{i},'_vs_age');
%     fact_des_mult_reg(subj_file, MS_file, img_path, var_names, stat_suff);
% end;
% spm_exec_mult_jobs;


%%% cvLBF regression analysis, 07/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% MS_file   = 'FADE_model_spaces/MS_FADE_02.mat';
% var_names ={'age'};
% lbf_names ={'GLM_1x1l_vs_0', 'GLM_1x1a_vs_0', 'GLM_1x2nf_vs_0', 'GLM_1x2nr_vs_0'};
% for i = 1:numel(lbf_names)
%     img_path  = strcat('model_selection/+/*_s8_',lbf_names{i},'_MC_cvLBF.nii');
%     stat_suff = strcat('LBF_',lbf_names{i},'_vs_age');
%     fact_des_mult_reg(subj_file, MS_file, img_path, var_names, stat_suff);
% end;
% spm_exec_mult_jobs;


%%% cvLBF meta batch, 07/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subj_file = 'subjects_all_2020_01_30.txt';
% cvLBF_meta_batch(subj_file, 'FADE_model_spaces/MS_FADE_02_no_x35.mat', 's8', ...
%                  [0 2 1 0 0 0, 0 0 0 0 0 0], 'GLM_1x1l_vs_0');
% cvLBF_meta_batch(subj_file, 'FADE_model_spaces/MS_FADE_02_no_x35.mat', 's8', ...
%                  [0 2 0 1 0 0, 0 0 0 0 0 0], 'GLM_1x1a_vs_0');
% cvLBF_meta_batch(subj_file, 'FADE_model_spaces/MS_FADE_02_no_x35.mat', 's8', ...
%                  [0 2 0 0 1 0, 0 0 0 0 0 0], 'GLM_1x2nf_vs_0');
% cvLBF_meta_batch(subj_file, 'FADE_model_spaces/MS_FADE_02_no_x35.mat', 's8', ...
%                  [0 2 0 0 0 1, 0 0 0 0 0 0], 'GLM_1x2nr_vs_0');


%%% second-level BMS, 05/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cvLMEs for test model space and subject
% cvLME_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_test.mat', 's8')

% cvBMS across 1st and 2nd model space
% subj_file = 'subjects_all_2020_01_30.txt';
% MS_names  = {'MS_FADE_01', 'MS_FADE_02'}';
% GLM_names = {'GLM3', 'GLM_1x1a'};
% cvBMS_across_spaces(subj_file, MS_names, GLM_names, 's8', 'MA_cvLME.nii', 'MS_FADEs_01_GLM3_vs_02_GLM_1x1a');
% spm_exec_mult_jobs;


%%% second-level BMS, 03/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cvBMS_analysis_batch('subjects_JS.xls', 'FADE_model_spaces/MS_FADE_02_fams_12x.mat', 's8', 'MA_cvLFE.nii', 'fams_12x');
% cvBMS_analysis_batch('subjects_JS.xls', 'FADE_model_spaces/MS_FADE_02_fams_x-1012.mat', 's8', 'MA_cvLFE.nii', 'fams_x-1012');
% cvBMS_analysis_batch('subjects_JS.xls', 'FADE_model_spaces/MS_FADE_02_fams_x1l1a.mat', 's8', 'MA_cvLFE.nii', 'fams_x1l1a');
% cvBMS_analysis_batch('subjects_JS.xls', 'FADE_model_spaces/MS_FADE_02_fams_x2nf2nr.mat', 's8', 'MA_cvLFE.nii', 'fams_x2nf2nr');


%%% cvLFE meta batch, 03/02/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cvLFE_meta_batch('subjects_JS.xls', 'FADE_model_spaces/MS_FADE_02_no_x35.mat', 's8', ...
%                  [1 1 1 1 1 1, 2 2 2 2 2 2], {'GLMs_1x', 'GLMs_2x'});
% cvLFE_meta_batch('subjects_JS.xls', 'FADE_model_spaces/MS_FADE_02_no_x35.mat', 's8', ...
%                  [1 2 3 3 4 4, 1 2 3 3 4 4], {'GLMs_x-1', 'GLMs_x0', 'GLMs_x1', 'GLMs_x2'});
% cvLFE_meta_batch('subjects_JS.xls', 'FADE_model_spaces/MS_FADE_02_no_x35.mat', 's8', ...
%                  [0 0 1 2 0 0, 0 0 1 2 0 0], {'GLMs_x1l', 'GLMs_x1a'});
% cvLFE_meta_batch('subjects_JS.xls', 'FADE_model_spaces/MS_FADE_02_no_x35.mat', 's8', ...
%                  [0 0 0 0 1 2, 0 0 0 0 1 2], {'GLMs_x2nf', 'GLMs_x2nr'});


%%% second-level BMS, 30/01/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collect_meta_batch('subjects_JS.xls', 'FADE_model_spaces/MS_FADE_02.mat', 's8', 'MA_cvLME.nii');
% cvBMS_analysis_batch('subjects_JS.xls', 'FADE_model_spaces/MS_FADE_02_only_1x-10.mat', 's8', 'MA_cvLME.nii', 'only_1x-10');
% cvBMS_analysis_batch('subjects_JS.xls', 'FADE_model_spaces/MS_FADE_02_only_1x01l.mat', 's8', 'MA_cvLME.nii', 'only_1x01l');
% cvBMS_analysis_batch('subjects_JS.xls', 'FADE_model_spaces/MS_FADE_02_only_1x1l1a.mat', 's8', 'MA_cvLME.nii', 'only_1x1l1a');


%%% behavioral data, 28/01/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% behavioral_data('subjects_BS_JS.xls');


%%% stats meta batch, 24/01/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stats_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_02.mat', 's8', 'prepare')
% stats_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_02.mat', 's8', 'perform')
% cvLME_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_02.mat', 's8');

% stats_meta_batch('subjects_rest.xls', 'FADE_model_spaces/MS_FADE_02.mat', 's8', 'prepare')
% stats_meta_batch('subjects_rest.xls', 'FADE_model_spaces/MS_FADE_02.mat', 's8', 'perform')
% cvLME_meta_batch('subjects_rest.xls', 'FADE_model_spaces/MS_FADE_02.mat', 's8');


%%% correlate estimates, 23/01/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GLM_dirs = {'C:\Joram\projects\DZNE\FADE\subjects_verio\af47\MS_FADE_test\s8_GLM3\';
%             'C:\Joram\projects\DZNE\FADE\subjects_verio\af47\MS_FADE_test\s8_GLM3n\';
%             'C:\Joram\projects\DZNE\FADE\subjects_verio\af47\MS_FADE_test\s8_GLM3o\';
%             'C:\Joram\projects\DZNE\FADE\subjects_verio\af47\MS_FADE_test\s8_GLM3on\'};
% beta_imgs = strcat(GLM_dirs,'beta_0002.nii');
% [R, P] = MF_correlate(beta_imgs);
% 
% load C:\Joram\projects\DZNE\FADE\subjects_verio\af47\MS_FADE_test\s8_GLM3\SPM.mat
% X  = SPM.xX.X(:,[1:3]);
% load C:\Joram\projects\DZNE\FADE\subjects_verio\af47\MS_FADE_test\s8_GLM3o\SPM.mat
% Xo = SPM.xX.X(:,[1:3]);
% figure; plot(X, Xo, '.');


%%% stats meta batch, 23/01/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stats_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_test.mat', 's8', 'prepare')
% stats_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_test.mat', 's8', 'perform')


%%% cvBMS results, 13/01/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BMS_dirs = {'C:\Joram\projects\DZNE\FADE\analyses_BS\model_selection\MS_FADE_01_only_n_young\';
%             'C:\Joram\projects\DZNE\FADE\analyses_BS\model_selection\MS_FADE_01_only_n_old\';
%             'C:\Joram\projects\DZNE\FADE\analyses_BS\model_selection\MS_FADE_01_only_n_23_all\'};
% ROI_imgs = {'C:\Joram\projects\DZNE\FADE\tools\FADE_ROI_images\rHippo_CA_Sub.nii';
%             'C:\Joram\projects\DZNE\FADE\tools\FADE_ROI_images\rHippocampus_AAL.nii'};
% for i = 1:numel(BMS_dirs)
%     cvBMS_results_report(BMS_dirs{i});
%     cvBMS_results_report(BMS_dirs{i}, ROI_imgs{1});
%     cvBMS_results_report(BMS_dirs{i}, ROI_imgs{2});
% end;


%%% stats meta batch, 08/01/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stats_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_01.mat', 's8', 'prepare')
% stats_meta_batch('subjects_test.xls', 'FADE_model_spaces/MS_FADE_01.mat', 's8', 'perform')
% 
% stats_meta_batch('subjects_rest.xls', 'FADE_model_spaces/MS_FADE_01.mat', 's8', 'prepare')
% stats_meta_batch('subjects_rest.xls', 'FADE_model_spaces/MS_FADE_01.mat', 's8', 'perform')


%%% first subject test, 08/01/2020 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scan_id = 1;
% subj_id = 'af47';
% MS_name = 'MS_FADE_01';
% create_onset_files(scan_id, subj_id, MS_name)
% 
% GLM_name = 'GLM1n';
% img_pref = 's8';
% create_stats_batch(scan_id, subj_id, MS_name, GLM_name, img_pref)