% FADE-SAME: Figure S3

% clear
% close all

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set results files
FADE_file = '../FADE_scores/FADE_SAME_scores_2021_01_11_FADE_only.xls';
covs_file = 'covariates/covariates_FADE.xls';

% define cohorts
cohs = {'young', 'middle-aged' ,'older', 'all'};
gens = {'E2/E2', 'E2/E3', 'E2/E4', 'E3/E3', 'E3/E4', 'E4/E4'};
cohs_cols = [1, 0, 0; 0, 1, 1; 0, 0, 1; 0, 0, 0];
gens_cols = ['bcmgyr'];

% define proportions
Ns   = [245, 300, 244, 206]';
perc = [0.8, 14.7, 3.6, 58.4, 18.4, 4.0; ...
        0.3, 13.9, 2.0, 55.2, 26.9, 1.6; ...
        0.8, 10.2, 2.0, 63.5, 21.3, 2.0; ...
        0.5, 10.0, 3.9, 63.7, 20.9, 0.5];
Na   = round(perc.*repmat(Ns,[1 numel(gens)])./100);
prop = sum(Na,1)./(sum(sum(Na)));
prop_old = [0.005, 0.13, 0.013, 0.6, 0.235, 0.017];
% Source: Li et al., Behavioral Genetics, 2019, vol. 49, pp. 455-468.
% Source (old): https://de.wikipedia.org/wiki/Apolipoprotein_E#Genetische_Varianten_und_ihr_Einfluss_auf_Erkrankungsrisiken.

% define variables
Y_vars = {'nov_FADE','nov_SAMe','mem_FADE','mem_SAMe'};
Y_labs = {'nov-FADE','nov-SAMe','mem-FADE','mem-SAMe'};
X_vars = {'age', 'ApoE'};
alpha  = 0.1;

% load results files
[num, txt, raw1] = xlsread(FADE_file);
[num, txt, raw2] = xlsread(covs_file);
clear num txt

% get results data
FADE_data = raw1(2:end,:);
covs_data = raw2(2:end,:);

% extract age data
num_subj  = size(FADE_data,1);
age       = cell2mat(FADE_data(:,4));
age_coh{1}= (age<50);
age_coh{2}= (age>=50 & age < 60);
age_coh{3}= (age>=60);
age_coh{4}= true(size(age));

% extract FADE/SAME scores
FADE_inds = 5+[1:4];
FADE_SAME = cell2mat(FADE_data(:,FADE_inds));
FADE_vars = raw1(1,FADE_inds);

% extract ApoE genotypes
ApoE = covs_data(:,6);


%%% Step 2: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract numbers
N = zeros(numel(cohs),numel(gens));
for j = 1:numel(cohs)
    for k = 1:numel(gens)
        N(j,k) = sum(strcmp(ApoE(age_coh{j}),gens{k}));
    end;
end;

% calculte frequencies
f = N./repmat(sum(N,2),[1 size(N,2)]);
f = [f; prop];

% perform chi^2 GOF tests
p     = zeros(1,numel(cohs));
bins  = [1:numel(gens)];
for j = 1:numel(cohs)
    [h, p(j), stats(j)] = chi2gof(bins, 'Ctrs', bins, 'Frequency', N(j,:), ...
                                        'Expected', prop*sum(N(j,:)), 'NParams', 0);
end;
clear h

% prepare between-subject ANOVAs
age_gr  = 1*(age<50) + 2*(age>=50 & age < 60) + 3*(age>=60);
ApoE_gt = 1*strcmp(ApoE,'E3/E3') + 2*strcmp(ApoE,'E3/E4') + 2*strcmp(ApoE,'E4/E4');
Y_inds  = (age_gr>0 & ApoE_gt>0);
Y = FADE_SAME(Y_inds,:);
X = [age_gr(Y_inds), ApoE_gt(Y_inds)];
T = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), X(:,1), X(:,2), 'VariableNames', [Y_vars, X_vars]);
n = size(T,1);
c = eye(3*2); % (size(Y,2));

% run between-subject ANOVAs
for i = 1:numel(Y_vars)
    % estimate means
    m(i).y =  Y(:,i);
    m(i).X = [X(:,1)==1 & X(:,2)==1, X(:,1)==1 & X(:,2)==2, ...
              X(:,1)==2 & X(:,2)==1, X(:,1)==2 & X(:,2)==2, ...
              X(:,1)==3 & X(:,2)==1, X(:,1)==3 & X(:,2)==2];
   [m(i).b_est, m(i).s2_est] = ME_GLM(m(i).y, m(i).X, eye(n));
    m(i).Bcov = (m(i).X'*m(i).X)^(-1);
    m(i).SEs  = sqrt(diag(c' * (m(i).s2_est*m(i).Bcov) * c));
    m(i).CIs  = m(i).SEs * norminv(mean([(1-alpha), 1]), 0, 1);
    % perform ANOVA
    m(i).model = fitlm(T, sprintf('%s ~ age*ApoE', Y_vars{i}));
    m(i).anova = anova(m(i).model);
end;

% generate ANOVA table
col = Y_labs;
row = m(1).anova.Properties.RowNames(1:end-1);
Res = cell(numel(row),numel(col));
for i = 1:numel(col)
    F = table2array(m(i).anova(1:end-1,4));
    q = table2array(m(i).anova(1:end-1,5));
    for j = 1:numel(row)
        % store F-/p-value
        if q(j) < 0.001
            Res{j,i} = sprintf('F = %1.2f, p < 0.001', F(j));
        else
            Res{j,i} = sprintf('F = %1.2f, p = %0.3f', F(j), q(j));
        end;
    end;
end;

% save ANOVA tables
Res = [cell(1,1), col; row, Res];
xlswrite('Figure_S3.xls', Res);


%%% Step 3: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot genotypes
figure('Name', 'FADE/SAME: ApoE genotypes', 'Color', [1 1 1], 'Position', [50 50 1600 900]);

for j = 1:size(f,1)
    if j < 4, subplot(2,3,3+j); end;
    if j > 3, subplot(2,3,j-2); end;
    hold on;
    for k = 1:numel(gens)
        bar(k, f(j,k), gens_cols(k));
    end;
    axis([(1-1), (numel(gens)+1), 0, 7/10]);
    set(gca,'Box','On');
    set(gca,'XTick',[1:numel(gens)],'XTickLabel',gens);
    if j == 5, legend(gens, 'Location', 'West'); end;
    xlabel('genotype', 'FontSize', 12);
    ylabel('frequency', 'FontSize', 12)
    if j  < 5, title(sprintf('%s subjects (N = %d)', cohs{j}, sum(N(j,:))), 'FontSize', 16); end;
    if j == 5, title('population distribution (Li et al., 2019)', 'FontSize', 16); end;
    for k = 1:numel(gens)
        if j < 5
            text(k, 2/3, num2str(N(j,k)), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        else
            text(k, 2/3, sprintf('%0.3f', prop(k)), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        end;
    end;
    if j < 5
        if p(j) < 0.001
            text((1-0.5), 1/2, sprintf('chi^2 = %2.2f, p < 0.001', stats(j).chi2stat), ...
                 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Middle');
        else
            text((1-0.5), 1/2, sprintf('chi^2 = %2.2f, p = %0.3f', stats(j).chi2stat, p(j)), ...
                 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Middle');
        end;
    end;
end;

% plot ANOVAs
figure('Name', 'FADE/SAME: ApoE ANOVAs', 'Color', [1 1 1], 'Position', [50 50 1600 900]);

for i = 1:numel(m)
    subplot(1,numel(m),i); hold on;
    for j = 1:max(age_gr)
        for k = 1:max(ApoE_gt)
            bar((j-1)*2+k, m(i).b_est((j-1)*2+k), 'FaceColor', (3/(2+k))*cohs_cols(j,:));
        end;
    end;
    errorbar([1:size(c,2)], m(i).b_est, m(i).CIs, '.k', 'LineWidth', 2, 'CapSize', 15);
    xlim([(1-0.5), (size(c,2)+0.5)]);
    ylim([-2.5, +0.25]);
    set(gca,'Box','On');
    set(gca,'XTick',[mean([1,2]), mean([3,4]), mean([5,6])],'XTickLabel',cohs(1:3),'XTickLabelRotation',0);
    xlabel('age group', 'FontSize', 12);
    if i == 1, ylabel(sprintf('parameter estimate and %d%% confidence intervals', ...
                      round((1-alpha)*100)), 'FontSize', 12);                    end;
    if i == 1, title(sprintf('novelty contrast:\nFADE score'), 'FontSize', 16);  end;
    if i == 2, title(sprintf('novelty contrast:\nSAME score'), 'FontSize', 16);  end;
    if i == 3, title(sprintf('memory contrast:\nFADE score'), 'FontSize', 16);   end;
    if i == 4, title(sprintf('memory contrast:\nSAME score'), 'FontSize', 16);   end;
    if i == 4, legend(repmat({'E3 homozygotes', 'E4 carriers'}, [1 3]), 'Location', 'SouthEast'); end;
end;