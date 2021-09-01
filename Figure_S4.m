% FADE-SAME: Figure S4

clear
% close all

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set results files
FADE_file = '../FADE_scores/FADE_SAME_scores_2021_01_11_FADE_only.xls';
covs_file = 'covariates/covariates_FADE.xls';

% define cohorts
cohs = {'middle-aged' ,'older'};
cohs_cols = [0, 1, 1; 0, 0, 1];
% cohs = {'young', 'middle-aged' ,'older'};
% cohs_cols = [1, 0, 0; 0, 1, 1; 0, 0, 1];

% define variables
Y_vars = {'nov_FADE','nov_SAMe','mem_FADE','mem_SAMe'};
Y_labs = {'nov-FADE','nov-SAMe','mem-FADE','mem-SAMe'};
X_vars = {'age', 'Abi'};
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

% extract FADE/SAME scores
FADE_inds = 5+[1:4];
FADE_SAME = cell2mat(FADE_data(:,FADE_inds));
FADE_vars = raw1(1,FADE_inds);

% extract ApoE genotypes
MWTB = cell2mat(covs_data(:,8));
MWTB(isnan(MWTB)) = 0;
Abi  = covs_data(:,7);
Abi  = 2-strcmp(Abi,'yes');


%%% Step 2: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compile data
Y = FADE_SAME(MWTB~=0,:);
X = MWTB(MWTB~=0);
x = 1*(age<50) + 2*(age>=50 & age < 60) + 3*(age>=60);
x = x(MWTB~=0);
p = max(x);
v = size(Y,2);

% prellocate data
x1 = cell(p,v);
y1 = cell(p,v);
mn1= cell(p,v);

% run correlation analyses
r1 = zeros(p,v);
p1 = zeros(p,v);
for i = 1:v
    for j = 1:max(x)
        y1{j,i} = Y(x==j,i);
        x1{j,i} = X(x==j);
        [r1(j,i), p1(j,i)] = corr(y1{j,i}, x1{j,i});
        mn1{j,i} = polyfit(x1{j,i}, y1{j,i}, 1);
    end;
end;

% prepare categorical analyses
age_gr = 1*(age>=50 & age < 60) + 2*(age>=60);
% age_gr = 1*(age<50) + 2*(age>=50 & age < 60) + 3*(age>=60);
Y_inds = (age_gr>0 & Abi>0);
Y = FADE_SAME(Y_inds,:);
X = [age_gr(Y_inds), Abi(Y_inds)];
T = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), X(:,1), X(:,2), 'VariableNames', [Y_vars, X_vars]);
n = size(T,1);
c = eye(2*2);
% c = eye(3*2);

% run two-sample t-tests
p2 = zeros(p,numel(Y_vars));
for i = 1:numel(Y_vars)
    for j = 1:max(age_gr)
        [h, p2(j,i)] = ttest2(Y(X(:,1)==j & X(:,2)==1,i), Y(X(:,1)==j & X(:,2)==2,i), 'Tail','both','Vartype','unequal');
    end;
end;

% run between-subject ANOVAs
for i = 1:numel(Y_vars)
    % estimate means
    m(i).y =  Y(:,i);
    m(i).X = [X(:,1)==1 & X(:,2)==1, X(:,1)==1 & X(:,2)==2, ...
              X(:,1)==2 & X(:,2)==1, X(:,1)==2 & X(:,2)==2];
  % m(i).X = [X(:,1)==1 & X(:,2)==1, X(:,1)==1 & X(:,2)==2, ...
  %           X(:,1)==2 & X(:,2)==1, X(:,1)==2 & X(:,2)==2, ...
  %           X(:,1)==3 & X(:,2)==1, X(:,1)==3 & X(:,2)==2];
   [m(i).b_est, m(i).s2_est] = ME_GLM(m(i).y, m(i).X, eye(n));
    m(i).Bcov = (m(i).X'*m(i).X)^(-1);
    m(i).SEs  = sqrt(diag(c' * (m(i).s2_est*m(i).Bcov) * c));
    m(i).CIs  = m(i).SEs * norminv(mean([(1-alpha), 1]), 0, 1);
    % perform ANOVA
    m(i).model = fitlm(T, sprintf('%s ~ age*Abi', Y_vars{i}));
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
xlswrite('Figure_S4.xls', Res);


%%% Step 3: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot correlations
figure('Name', 'FADE/SAME: MWT-B hits', 'Color', [1 1 1], 'Position', [50 50 1280 440]);
p1_thr = [0.05, 0.01, 0.001];
x_lim  = [17, 37];
x_lab  = {'MWT-B hits'};

for i = 1:v
    
    % young and older subjects
    subplot(2,v,i); hold on;
    plot(x1{1,i}, y1{1,i}, '.r', 'Color', [1, (p1(1,i)>p1_thr(1))*1/2, (p1(1,i)>p1_thr(1))*1/2]);
    plot(x1{3,i}, y1{3,i}, '.b', 'Color', [(p1(3,i)>p1_thr(1))*1/2, (p1(3,i)>p1_thr(1))*1/2, 1]);
    plot(x_lim, mn1{1,i}(1)*x_lim + mn1{1,i}(2), '-r', 'Color', [1, (p1(1,i)>p1_thr(1))*1/2, (p1(1,i)>p1_thr(1))*1/2]);
    plot(x_lim, mn1{3,i}(1)*x_lim + mn1{3,i}(2), '-b', 'Color', [(p1(3,i)>p1_thr(1))*1/2, (p1(3,i)>p1_thr(1))*1/2, 1]);
    xlim(x_lim);
    ylim([min(y1{3,i})-1/20*range(y1{3,i}), max(y1{1,i})+1/20*range(y1{1,i})]);
    set(gca,'Box','On');
    xlabel(x_lab{1}, 'FontSize', 12);
    ylabel('score', 'FontSize', 12);
    if i == 1, title(sprintf('novelty contrast: FADE score\n '), 'FontSize', 12); end;
    if i == 2, title(sprintf('novelty contrast: SAME score\n '), 'FontSize', 12); end;
    if i == 3, title(sprintf('memory contrast: FADE score\n '), 'FontSize', 12);  end;
    if i == 4, title(sprintf('memory contrast: SAME score\n '), 'FontSize', 12);  end;
    if p1(1,i) < p1_thr(end)
        txt1 = sprintf('r = %0.2f, p < 0.001', r1(1,i));
    else
        txt1 = sprintf('r = %0.2f, p = %0.3f', r1(1,i), p1(1,i));
    end;
    if p1(3,i) < p1_thr(end)
        txt2 = sprintf('r = %0.2f, p < 0.001', r1(3,i));
    else
        txt2 = sprintf('r = %0.2f, p = %0.3f', r1(3,i), p1(3,i));
    end;
    text(x_lim(1), max(y1{1,i})+1/20*range(y1{1,i}), txt1, ...
         'FontSize', 8, 'FontWeight', 'Bold', 'Color', [1, (p1(1,i)>p1_thr(1))*1/2, (p1(1,i)>p1_thr(1))*1/2], ...
         'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Left');
    text(x_lim(2), max(y1{1,i})+1/20*range(y1{1,i}), txt2, ...
         'FontSize', 8, 'FontWeight', 'Bold', 'Color', [(p1(3,i)>p1_thr(1))*1/2, (p1(3,i)>p1_thr(1))*1/2, 1], ...
         'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Right');
     
    % middle-aged subjects
    subplot(2,v,v+i); hold on;
    plot(x1{2,i}, y1{2,i}, '.b', 'Color', [(p1(2,i)>p1_thr(1))*2/3, 1, 1]);
    plot(x_lim, mn1{2,i}(1)*x_lim + mn1{2,i}(2), '-b', 'Color', [(p1(2,i)>p1_thr(1))*2/3, 1, 1]);
    xlim(x_lim);
    ylim([min(y1{2,i})-1/20*range(y1{2,i}), max(y1{2,i})+1/20*range(y1{2,i})]);
    set(gca,'Box','On');
    xlabel(x_lab{1}, 'FontSize', 12);
    ylabel('score', 'FontSize', 12);
    if p1(2,i) < p1_thr(end)
        txt2 = sprintf('r = %0.2f, p < 0.001', r1(2,i));
    else
        txt2 = sprintf('r = %0.2f, p = %0.3f', r1(2,i), p1(2,i));
    end;
    text(x_lim(2), max(y1{2,i})+1/20*range(y1{2,i}), txt2, ...
         'FontSize', 8, 'FontWeight', 'Bold', 'Color', [(p1(2,i)>p1_thr(1))*2/3, 1, 1], ... % [(p1(2,i)>p1_thr(1))*0, 1, 1], ...
         'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Right');
    
end;

% plot ANOVAs
figure('Name', 'FADE/SAME: Abitur ANOVAs', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
p2_thr  = [0.05, 0.05/(numel(cohs)), 0.05/(numel(cohs)*numel(Y_vars))];

for i = 1:numel(m)
    
    subplot(1,numel(m),i); hold on;
    for j = 1:max(age_gr)
        for k = 1:max(Abi)
            bar((j-1)*2+k, m(i).b_est((j-1)*2+k), 'FaceColor', (3/(2+k))*cohs_cols(j,:));
        end;
        if p2(j,i) < p2_thr(3)
            text((j-1)*2+1.5, 1/8, '***', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        elseif p2(j,i) < p2_thr(2)
            text((j-1)*2+1.5, 1/8, '**', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        elseif p2(j,i) < p2_thr(1)
            text((j-1)*2+1.5, 1/8, '*', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        else
            text((j-1)*2+1.5, 1/8, 'n.s.', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        end;
    end;
    errorbar([1:size(c,2)], m(i).b_est, m(i).CIs, '.k', 'LineWidth', 2, 'CapSize', 15);
    xlim([(1-0.5), (size(c,2)+0.5)]);
    ylim([-2.5, +0.25]);
    set(gca,'Box','On');
    set(gca,'XTick',[mean([1,2]), mean([3,4])],'XTickLabel',cohs,'XTickLabelRotation',0);
  % set(gca,'XTick',[mean([1,2]), mean([3,4]), mean([5,6])],'XTickLabel',cohs,'XTickLabelRotation',0);
    xlabel('age group', 'FontSize', 12);
    if i == 1, ylabel(sprintf('parameter estimate and %d%% confidence intervals', ...
                      round((1-alpha)*100)), 'FontSize', 12);                    end;
    if i == 1, title(sprintf('novelty contrast:\nFADE score'), 'FontSize', 16);  end;
    if i == 2, title(sprintf('novelty contrast:\nSAME score'), 'FontSize', 16);  end;
    if i == 3, title(sprintf('memory contrast:\nFADE score'), 'FontSize', 16);   end;
    if i == 4, title(sprintf('memory contrast:\nSAME score'), 'FontSize', 16);   end;
    if i == 4, legend(repmat({'with Abitur', 'w/o Abitur'}, [1 2]), 'Location', 'SouthEast'); end;
    
end;