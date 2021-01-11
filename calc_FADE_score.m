function FADE = calc_FADE_score(SPM_mat, con_vec, con_type, ref_map, do_plot)
% _
% Calculate FADE score for a single subject using reference map [1]
%     SPM_mat  - filepath to SPM.mat of estimated GLM
%     con_vec  - the memory contrast vector (e.g. [+1 -1 0])
%     con_type - the memory contrast type (i.e. 't' or 'F')
%     ref_map  - filename of reference map specifying VOI
%     do_plot  - logical indicating plotting of histograms
% 
%     FADE     - the FADE score, a scalar quantifying the
%               "functional activity deviation during encoding"
% 
% [1] Düzel E, Schütze H, Yonelinas AP, Heinze HJ (2010).
%     Functional Phenotyping of Successful Again in Long-Term Memory.
%     Hippocampus, vol. 21, pp. 803-814.
% 
% written by Joram Soch <Joram.Soch@DZNE.de>, 24/02/2020, 14:57 (V1);
% adapted: 27/02/2020, 20:27; 28/02/2020, 12:40 (V2)


% SPM_mat  = 'C:\Joram\projects\DZNE\FADE\analyses_BS\group_statistics\MS_FADE_04\ad30_s6w_GLM_TD_1th-a_SPM.mat';
% con_vec  = [ 0  1  0];
% con_type = 't';
% ref_map  = 'C:\Joram\projects\DZNE\FADE\analyses_BS\group_statistics\MS_FADE_04_FADE_GLM_1a_novelty_subjects_all_2020_11_05_young_G1\con_0002_FWE_0.05_10.nii';
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

% load reference map
ref_hdr = spm_vol(ref_map);
ref_img = spm_read_vols(ref_hdr);
R       = reshape(ref_img,[1 v]);
clear ref_hdr ref_img


%%% Step 2: calculate FADE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get voxel indices
VOI_ind = find(M~=0 & R~=0);                % volume of interest indices
VOnI_ind= find(M~=0 & R==0);                % volume of no interest inds

% calculate t-map
if strcmp(con_type,'t')
    covB = SPM.xX.Bcov;
    %    = inv(SPM.xX.X'*inv(SPM.xVi.V)*SPM.xX.X);
    t    = c'*B ./ sqrt( S2 * (c'*covB*c) );
end;
clear covB

% calculate F-map
if strcmp(con_type,'F')
    C    = c;
    q    = size(C,2);
    CB   = C'*B;
    covB = SPM.xX.Bcov;
    invC = inv(C'*covB*C);
    % F-contrast vector
    if q == 1
        F = (CB * invC .* CB) ./ S2;
    % F-contrast matrix
    elseif q > 1
        numF = zeros(1,v);
        for j = 1:v
            numF(j) = CB(:,j)' * invC * CB(:,j);
        end;
        F = (1/q) * numF./S2;
    end;
end;
clear C q CB covB invC numF

% calculate FADE score
if strcmp(con_type,'t')
    stat_in  = t(VOI_ind);      % t-statistics inside VOI
    stat_out = t(VOnI_ind);     % t-statistics outside VOI
    FADE     = mean(stat_out) - mean(stat_in);
elseif strcmp(con_type,'F')
    stat_in  = F(VOI_ind);      % F-statistics inside VOI
    stat_out = F(VOnI_ind);     % F-statistics outside VOI
    FADE     = median(stat_out) - median(stat_in);
end;

% plot t/F-distribution
if do_plot
    figure('Name', 'calc_FADE_score', 'Color', [1 1 1], 'Position', [50 50 1280 720]);
    % outside VOI
    subplot(1,2,1);
    hist(stat_out,100);
    if strcmp(con_type,'t'), title(sprintf('out: mean = %2.2f', mean(stat_out)), 'FontSize', 16); end;
    if strcmp(con_type,'F'), title(sprintf('out: median = %2.2f', median(stat_out)), 'FontSize', 16); end;
    % inside VOI
    subplot(1,2,2);
    hist(stat_in,100);
    if strcmp(con_type,'t'), title(sprintf('in: mean = %2.2f', mean(stat_in)), 'FontSize', 16); end;
    if strcmp(con_type,'F'), title(sprintf('in: median = %2.2f', median(stat_in)), 'FontSize', 16); end;
end;