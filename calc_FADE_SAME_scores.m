function calc_FADE_SAME_scores(stud_dir, subj_file, GLM_name, ref_dirs, ref_maps, con_vecs, contrasts, scores, res_suff)
% _
% This script calculates (i) FADE and SAMe score for (ii) novelty and
% memory contrast from (iii) AiA and yFADE subjects using (iv) cross-
% validated reference maps.
% 
% written by Joram Soch <Joram.Soch@DZNE.de>, 27/07/2020, 15:51 (V1);
% adapted: 03/08/2020, 08:36 (V2); 17/08/2020, 11:16 (V3);
% edited for upload: 11/01/2021, 16:34


%%% Step 1: load subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load subjects
[num, txt, raw] = xlsread(subj_file);
tab_hdr   = raw(1,:);
subj_ids  = raw(2:end,1);
num_subj  = numel(subj_ids);
clear num txt

% get covariates
subj_info = raw(2:end,:);
stud_ids  = 2-strncmp(subj_ids,'subA',4);
grp_ids   = cell2mat(subj_info(:,strcmp(tab_hdr,'CV group')));


%%% Step 2: calculate scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preallocate scores
FADE_SAME = zeros(num_subj,numel(contrasts)*numel(scores));
res_hdr   = cell(1,numel(contrasts)*numel(scores));
fprintf('\n-> Calculate FADE/SAME scores:\n');

% for all subjects
for i = 1:num_subj
    
    fprintf('   - Subject %s (%d out of %d):\n', subj_ids{i}, i, num_subj);
    
    % for all contrasts
    for j = 1:numel(contrasts)
        
        fprintf('     - Contrast "%s" (%d out of %d):\n', contrasts{j}, j, numel(contrasts));
        
        % for all scores
        for k = 1:numel(scores)
            
            fprintf('       - Score "%s" (%d out of %d) ...', scores{k}, k, numel(scores));
        
            % prepare FADE/SAME score
            SPM_mat = strcat(stud_dir,'beta_images','/',subj_ids{i},'_',GLM_name,'_','SPM.mat');
            con_vec = con_vecs{strcmp(contrasts,contrasts{j})};
            if size(ref_dirs,1) > 1         % study-dependent relative directory names
                SPM_ref = strcat(stud_dir,'SPM_analyses','/',...
                                 ref_dirs{stud_ids(i),3-grp_ids(i)},'/','SPM.mat');
            else                            % study-independent absolute directory names
                SPM_ref = strcat(ref_dirs{1,3-grp_ids(i)},'/','SPM.mat');
            end;
            SPM_ref = strcat(SPM_ref(1:strfind(SPM_ref,'~')-1),contrasts{j},...
                             SPM_ref(strfind(SPM_ref,'~')+1:end));
                         
            % calculate FADE/SAME score
            switch scores{k}
                case 'FADE-pos'
                    ref_map = strcat(fileparts(SPM_ref),'/',ref_maps{k});
                    score   = calc_FADE_score(SPM_mat, +1*con_vec, 't', ref_map, false);
                    FADE_pos= score;
                    ref_pos = ref_map;
                case 'FADE-neg'
                    ref_map = strcat(fileparts(SPM_ref),'/',ref_maps{k});
                    score   = calc_FADE_score(SPM_mat, -1*con_vec, 't', ref_map, false);
                    FADE_neg= score;
                    ref_neg = ref_map;
                case 'FADE-avg'
                    score   = mean([FADE_pos, FADE_neg]);
                    clear  FADE_pos FADE_neg
                case 'FADE-two'
                    score   = calc_FADE2_score(SPM_mat, con_vec, ref_pos, ref_neg, false);
                    clear  ref_pos ref_neg
                case 'SAME'
                    ref_map = ref_maps{k};
                    score   = calc_SAME_score(SPM_mat, con_vec, SPM_ref, ref_map, false);
            end;
            
            % store FADE/SAME score
            if i == 1
                res_hdr{(j-1)*numel(scores)+k} = strcat(contrasts{j},'_',scores{k});
            end;
            FADE_SAME(i,(j-1)*numel(scores)+k) = score;
            fprintf('successful!\n');
            clear score
            
        end;
        
    end;
    
end;


%%% Step 3: save scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create results table
fprintf('\n');
res_tab = [ [tab_hdr; subj_info], [res_hdr; num2cell(FADE_SAME)] ];
res_dir = strcat(stud_dir,'FADE_scores/');
if ~exist(res_dir,'dir'), mkdir(res_dir); end;

% save FADE/SAME scores
today    = datevec(datetime);
if ~isempty(res_suff), res_suff = strcat('_',res_suff); end;
filename = strcat(res_dir,'FADE_SAME_scores_',...
                  num2str(today(1)),'_',sprintf('%02d',today(2)),'_',sprintf('%02d',today(3)),res_suff,'.xls');
xlswrite(filename, res_tab);

% save extraction settings
filename = strcat(res_dir,'FADE_SAME_scores_',...
                  num2str(today(1)),'_',sprintf('%02d',today(2)),'_',sprintf('%02d',today(3)),res_suff,'.mat');
save(filename, 'stud_dir', 'subj_file', 'GLM_name', 'ref_dirs', 'ref_maps', ...
               'con_vecs', 'contrasts', 'scores',   'res_suff');