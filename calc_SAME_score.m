function SAMe = calc_SAME_score(SPM_mat, con_vec, SPM_ref, ref_map, do_plot)
% _
% Calculate SAMe score for a single subject using reference map
%     SPM_mat  - filepath to SPM.mat of estimated GLM
%     con_vec  - the memory contrast vector (e.g. [+1 -1 0])
%     SPM_ref  - filepath to SPM.mat of reference GLM
%     ref_map  - relative filepath of reference SPM [1]
%     do_plot  - logical indicating plotting of histograms
% 
%     SAMe     - the SAMe score, a scalar quantifying the
%               "successful aging in memory (processing)"
% 
% [1] Düzel E, Schütze H, Yonelinas AP, Heinze HJ (2010).
%     Functional Phenotyping of Successful Again in Long-Term Memory.
%     Hippocampus, vol. 21, pp. 803-814.
% 
% written by Joram Soch <Joram.Soch@DZNE.de>, 23/03/2020, 23:11


% SPM_mat  = 'C:\Joram\projects\DZNE\FADE\analyses_BS\group_statistics\MS_FADE_04\ad30_s6w_GLM_TD_1th-a_SPM.mat';
% con_vec  = [ 0  1  0];
% SPM_ref  = 'C:\Joram\projects\DZNE\FADE\analyses_BS\group_statistics\MS_FADE_04_FADE_GLM_1a_novelty_subjects_all_2020_11_05_young_G1\SPM.mat';
% ref_map = 'con_0001_FWE_0.05_10.nii';
% do_plot  = true;


%%% Step 1a: load estimated GLM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load SPM.mat
GLM_dir = SPM_mat(1:strfind(SPM_mat,'SPM.mat')-1);
load(SPM_mat);

% get model properties
n = numel(SPM.xY.VY);
p = numel(SPM.Vbeta);
v = prod(SPM.VM.dim);
c = [con_vec, zeros(1,p-numel(con_vec))]';

% load mask image
SPM.VM.fname = strcat(GLM_dir,SPM.VM.fname);
m_img = spm_read_vols(SPM.VM);
m_img = reshape(m_img,[1 v]);
M     = m_img;
clear m_img

% load beta estimates
B = zeros(p,v);
for i = find(c)'
    SPM.Vbeta(i).fname = strcat(GLM_dir,SPM.Vbeta(i).fname);
    b_img = spm_read_vols(SPM.Vbeta(i));
    B(i,:)= reshape(b_img,[1 v]);
end;
clear b_img


%%% Step 1b: load reference GLM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load SPM.mat
load(SPM_ref);
ref_dir = strcat(SPM.swd,'/');

% load reference map
ref_hdr = spm_vol(strcat(ref_dir,ref_map));
ref_img = spm_read_vols(ref_hdr);
R       = reshape(ref_img,[1 v]);
clear ref_hdr ref_img

% load beta estimate
SPM.Vbeta(1).fname = strcat(ref_dir,SPM.Vbeta(1).fname);
b1_img = spm_read_vols(SPM.Vbeta(1));
B1     = reshape(b1_img,[1 v]);
clear b1_img

% load residual variance
SPM.VResMS.fname = strcat(ref_dir,SPM.VResMS.fname);
s2_img = spm_read_vols(SPM.VResMS);
S2     = reshape(s2_img,[1 v]);
clear s2_img


%%% Step 2: calculate SAMe score %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get voxel indices
Jp = find(M~=0 & R~=0 & B1>0);         % positive memory effects
Jm = find(M~=0 & R~=0 & B1<0);         % negative memory effects

% calculate contrast
C1 = c'*B;

% calculate SAMe score
SAMe_pos = (C1(Jp)-B1(Jp))./sqrt(S2(Jp));
SAMe_neg = (B1(Jm)-C1(Jm))./sqrt(S2(Jm));
SAMe = 1/numel(Jp)*sum(SAMe_pos) + 1/numel(Jm)*sum(SAMe_neg);

% plot positive/negative
if do_plot
    figure('Name', 'calc_SAMe_score', 'Color', [1 1 1], 'Position', [50 50 1280 720]);
    % positive effect voxels
    subplot(1,2,1);
    hist(SAMe_pos,100);
    title(sprintf('pos: numel = %d, mean = %2.2f', numel(Jp), mean(SAMe_pos)), 'FontSize', 16);
    % negative effect voxels
    subplot(1,2,2);
    hist(SAMe_neg,100);
    title(sprintf('neg: numel = %d, mean = %2.2f', numel(Jm), mean(SAMe_neg)), 'FontSize', 16);
end;