function FADE = calc_FADE2_score(SPM_mat, con_vec, ref_pos, ref_neg, do_plot)
% _
% Calculate two-sided FADE score for a single subject using reference maps
%     SPM_mat  - filepath to SPM.mat of estimated GLM
%     con_vec  - the memory contrast vector (e.g. [+1 -1 0])
%     ref_pos  - filename of reference map for positive effects
%     ref_neg  - filename of reference map for negative effects
%     do_plot  - logical indicating plotting of histograms
% 
%     FADE     - the FADE score, a scalar quantifying the
%               "functional activity deviation during encoding"
% 
% [1] Düzel E, Schütze H, Yonelinas AP, Heinze HJ (2010).
%     Functional Phenotyping of Successful Again in Long-Term Memory.
%     Hippocampus, vol. 21, pp. 803-814.
% 
% written by Joram Soch <Joram.Soch@DZNE.de>, 03/08/2020, 08:33


% SPM_mat  = 'C:\Joram\projects\DZNE\FADE\analyses_BS\group_statistics\MS_FADE_04\ae09_s6w_GLM_TD_1th-a_SPM.mat';
% con_vec  = [0 1 0];
% ref_pos  = 'C:\Joram\projects\DZNE\FADE\analyses_BS\group_statistics\MS_FADE_04_FADE_GLM_1a_memory_subjects_all_2020_07_23_young_G1\con_0002_FWE_0.05_10.nii';
% ref_neg  = 'C:\Joram\projects\DZNE\FADE\analyses_BS\group_statistics\MS_FADE_04_FADE_GLM_1a_memory_subjects_all_2020_07_23_young_G1\con_0003_FWE_0.05_10.nii';
% do_plot  = true;


%%% Step 1: load estimates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% load residual variance
SPM.VResMS.fname = strcat(GLM_dir,SPM.VResMS.fname);
s2_img = spm_read_vols(SPM.VResMS);
S2     = reshape(s2_img,[1 v]);
clear s2_img

% load reference maps
ref_hdr = spm_vol(ref_pos);
ref_img = spm_read_vols(ref_hdr);
Rp      = reshape(ref_img,[1 v]);
ref_hdr = spm_vol(ref_neg);
ref_img = spm_read_vols(ref_hdr);
Rn      = reshape(ref_img,[1 v]);
clear ref_hdr ref_img


%%% Step 2: calculate FADE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get voxel indices
VOIp_ind = find(M~=0 & Rp~=0);  % positive volume of interest indices
VOIn_ind = find(M~=0 & Rn~=0);  % negative volume of interest indices
VOnI_ind = find(M~=0 & Rp==0 & Rn==0);      % volume of no interest

% calculate t-map
covB = SPM.xX.Bcov;
%    = inv(SPM.xX.X'*inv(SPM.xVi.V)*SPM.xX.X);
t    = c'*B ./ sqrt( S2 * (c'*covB*c) );
clear covB

% calculate FADE score
statp_in = t(VOIp_ind);         % t-statistics inside positive VOI
statn_in =-t(VOIn_ind);         % t-statistics inside negative VOI
stat_out = t(VOnI_ind);         % t-statistics outside VOIs
FADE     = mean(stat_out) - mean(statp_in) - mean(statn_in);

% plot t-distribution
if do_plot
    figure('Name', 'calc_FADE2_score', 'Color', [1 1 1], 'Position', [50 50 1280 720]);
    % outside VOI
    subplot(1,3,1);
    hist(stat_out,100);
    title(sprintf('out: mean = %2.2f', mean(stat_out)), 'FontSize', 16);
    % inside positive VOI
    subplot(1,3,2);
    hist(statp_in,100);
    title(sprintf('in (pos): mean = %2.2f', mean(statp_in)), 'FontSize', 16);
    % inside positive VOI
    subplot(1,3,3);
    hist(statn_in,100);
    title(sprintf('in (neg): mean = %2.2f', mean(statn_in)), 'FontSize', 16);
end;