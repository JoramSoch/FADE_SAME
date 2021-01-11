function contrast_meta_batch(stud_dir, subj_file, GLM_name, con_vec, con_img)
% _
% Contrast meta batch for analysis of the FADE data set
%     stud_dir  - study directory in data/analyses are located
%     subj_file - an XLS file containing subject IDs and covariates
%     GLM_name  - a string indicating the GLM of interest
%     con_vec   - a  1 x p vector of contrast weights
%     con_img   - a string indicating the contrast image name
% 
% written by Joram Soch <Joram.Soch@DZNE.de>, 17/02/2020, 12:26;
% edited for upload: 11/01/2021, 12:23


% load subjects
[num, txt, raw] = xlsread(subj_file);
subj_ids = raw(2:end,1);
num_subj = numel(subj_ids);
clear num txt

% specify directory
beta_dir = strcat(stud_dir,'beta_images');
fprintf('\n');

% for all subjects
for i = 1:num_subj
    
    % display subject
    fprintf('-> Estimate contrast for subject "%s" (%d out of %d) ... ', subj_ids{i}, i, num_subj);
    
    % get image header
    beta_imgs = dir(strcat(beta_dir,'/',subj_ids{i},'_',GLM_name,'_beta_*.nii'));
    beta_hdr  = spm_vol(strcat(beta_dir,'/',beta_imgs(1).name));
    clear beta_imgs
    
    % load beta images
    B = zeros(numel(con_vec),prod(beta_hdr.dim));
    for j = 1:numel(con_vec)
        filename = strcat(beta_dir,'/',subj_ids{i},'_',GLM_name,'_beta_',sprintf('%04d',j),'.nii');
        if exist(filename,'file')
            b_hdr = spm_vol(filename);
            b_img = spm_read_vols(b_hdr);
            B(j,:)= reshape(b_img, [1 prod(beta_hdr.dim)]);
        end;
    end;
    clear b_hdr b_img
    
    % write contrast image
    C = con_vec*B;
    con_hdr         = beta_hdr;
    con_hdr.fname   = strcat(beta_dir,'/',subj_ids{i},'_',GLM_name,'_',con_img);
    con_hdr.descrip = sprintf('contrast_meta_batch: Subject "%s", GLM "%s", c = [%s]', subj_ids{i}, GLM_name, num2str(con_vec));
    spm_write_vol(con_hdr, reshape(C,beta_hdr.dim));
    fprintf('successful!\n');
    
end;

% change to tools directory
fprintf('\n');
cd(strcat(stud_dir,'FADE_SAME/'));