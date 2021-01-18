% FADE-SAME: Figure 3

% clear
% close all

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load results file
FADE_file = '../FADE_scores/FADE_SAME_scores_2021_01_11_FADE_only.xls';
[num, txt, raw] = xlsread(FADE_file);
clear num txt

% specify inference
alpha = 0.1;

% extract variables
Y = cell2mat(raw(2:end,5+[1:4]));
X = cell2mat(raw(2:end,[2,3,4]));
X(:,3) = 1*double(X(:,3)<50) + 2*double(X(:,3)>=60);
Y = Y(X(:,3)~=0,:);
X = X(X(:,3)~=0,:);
[n, v] = size(Y);
[n, p] = size(X);

% create table
Y_vars = {'nov_FADE','nov_SAMe','mem_FADE','mem_SAMe'};
Y_labs = {'nov-FADE','nov-SAMe','mem-FADE','mem-SAMe'};
X_vars = {'scanner','gender','age'};
T = table(Y(:,1), Y(:,2), Y(:,3), Y(:,4), X(:,1), X(:,2), X(:,3), ...
    'VariableNames', [Y_vars, X_vars]);


%%% Step 2: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% within-subject design matrix
X2 = double([X(:,3)==1, X(:,3)==2]);
X2 = blkdiag(X2,X2);
X2 = [X2, repmat(eye(n),[2 1])];
l2 = {'FADE/young', 'FADE/older', 'SAME/young', 'SAME/older'};
c2 = [eye(2^2); zeros(size(X2,2)-2^2,2^2)];

% within-subject ANOVAs
for i = 1:2
    m2(i).Y_ind = (i-1)*2+[1:2];
    m2(i).Y_var = Y_vars(m2(i).Y_ind);
    m2(i).X_var = X_vars(3);
    % estimate parameters
    y = reshape(Y(:,m2(i).Y_ind),[2*n 1]);
   [m2(i).b_est, m2(i).s2_est] = ME_GLM_pinv(y, X2, eye(size(X2,1)));
    m2(i).Bcov  = pinv(X2)*pinv(X2)';
    m2(i).SEs   = sqrt(diag(c2' * (m2(i).s2_est*m2(i).Bcov) * c2));
    m2(i).CIs   = m2(i).SEs * norminv(1-alpha/2, 0, 1);
    % perform ANOVA
    m2(i).wsf   = table([1 2]', 'VariableNames', {'score'});
    m2(i).model = fitrm(T, sprintf('%s,%s ~ age', m2(i).Y_var{1}, m2(i).Y_var{2}), 'WithinDesign', m2(i).wsf);
    m2(i).anova = ranova(m2(i).model, 'WithinModel', 'score');
end;

% two-sample t-tests
for i = 1:2
    for j = 1:2
        y  = Y(:, m2(i).Y_ind(j));
        Xj = X2(1:size(X2,1)/2, 1:2);
        m2(i).t2(j).b = ME_GLM(y, Xj, eye(size(Xj,1)));
        y1 = y(Xj(:,1)==1);
        y2 = y(Xj(:,2)==1);
    [h, m2(i).t2(j).p, ci, stats] = ttest2(y1, y2);
        n1 = numel(y1);
        n2 = numel(y2);
        m2(i).t2(j).dp = (m2(i).t2(j).b(1)-m2(i).t2(j).b(2))/sqrt(((n1-1)*var(y1)+(n2-1)*var(y2))./(n1+n2-1));
    end;
end;
clear h ci stats


%%% Step 3: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% within-subject ANOVAs
figure('Name', 'FADE and SAME: ANOVAs', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

for i = 1:numel(m2)
    subplot(1,numel(m2),i); hold on;
    bar([1:4], [m2(i).b_est(1), 0, m2(i).b_est(3), 0], 'r');
    bar([1:4], [0, m2(i).b_est(2), 0, m2(i).b_est(4)], 'b');
    errorbar([1:size(c2,2)], m2(i).b_est(1:size(c2,2)), m2(i).CIs, '.k', 'LineWidth', 2, 'CapSize', 15);
    xlim([(1-0.5), (size(c2,2)+0.5)]);
    ylim([min([0, min(m2(i).b_est(1:size(c2,2))-m2(i).CIs)])-(1/10)*range(m2(i).b_est(1:size(c2,2))), ...
          max([0, max(m2(i).b_est(1:size(c2,2))+m2(i).CIs)])+(1/10)*range(m2(i).b_est(1:size(c2,2)))]);
    set(gca,'Box','On');
    set(gca,'XTick',[mean([1,2]), mean([3,4])],'XTickLabel',{'FADE', 'SAME'});
    legend('young subjects', 'older subjects', 'Location', 'SouthEast');
    xlabel('fMRI score', 'FontSize', 12);
    ylabel(sprintf('parameter estimate and %d%% confidence intervals', round((1-alpha)*100)), 'FontSize', 12);
    if i == 1, title('novelty contrast', 'FontSize', 16); end;
    if i == 2, title('memory contrast', 'FontSize', 16); end;
    for j = 1:numel(m2(i).t2)
        if m2(i).t2(j).p >= 0.001
            text((j-1)*2+1.5, (1/2)*(max([0, max(m2(i).b_est(1:size(c2,2))+m2(i).CIs)])+(1/10)*range(m2(i).b_est(1:size(c2,2)))), ...
                 sprintf('d'' = %0.2f, p = %0.3f', m2(i).t2(j).dp, m2(i).t2(j).p), ...
                 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        else
            text((j-1)*2+1.5, (1/2)*(max([0, max(m2(i).b_est(1:size(c2,2))+m2(i).CIs)])+(1/10)*range(m2(i).b_est(1:size(c2,2)))), ...
                 sprintf('d'' = %0.2f, p < 0.001', m2(i).t2(j).dp), ...
                 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        end;
    end;
end;