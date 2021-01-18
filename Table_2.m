% FADE-SAME: Table 2

% clear
% close all

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load results file
FADE_file = '../FADE_scores/FADE_SAME_scores_2021_01_11_FADE_only.xls';
[num, txt, raw] = xlsread(FADE_file);
clear num txt

% specify inference
conf_lvl = 0.9;

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
lvls2 = {'FADE',  'SAMe';
         'young', 'old'}';
X2 = double([X(:,3)==1, X(:,3)==2]);
X2 = blkdiag(X2,X2);
X2 = [X2, repmat(eye(n),[2 1])];
l2 = {'FADE/young', 'FADE/old', 'SAMe/young', 'SAMe/old'};
c2 = [eye(2^2); zeros(size(X2,2)-2^2,2^2)];

% within-subject ANOVAs
for i = 1:2
    m2(i).Y_ind = (i-1)*2+[1:2];
    m2(i).Y_var = Y_vars(m2(i).Y_ind);
    m2(i).X_var = X_vars(3);
    % estimate parameters
    y = reshape(Y(:,m2(i).Y_ind),[2*n 1]);
    m2(i).b_est = pinv(X2)*y;
    m2(i).s2_est= mean((y-X2*m2(i).b_est).^2);
    m2(i).Bcov  = pinv(X2)*pinv(X2)';
    m2(i).SEs   = sqrt(diag(c2' * (m2(i).s2_est*m2(i).Bcov) * c2));
    m2(i).CIs   = m2(i).SEs * norminv(mean([conf_lvl 1]), 0, 1);
    % perform ANOVA
    m2(i).wsf   = table([1 2]', 'VariableNames', {'score'});
    m2(i).model = fitrm(T, sprintf('%s,%s ~ age', m2(i).Y_var{1}, m2(i).Y_var{2}), 'WithinDesign', m2(i).wsf);
    m2(i).anova = ranova(m2(i).model, 'WithinModel', 'score');
end;


%%% Step 3: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate table (fitlm)
col = {'novelty', 'memory'};
row = m2(1).anova.Properties.RowNames([2,4,5]);
Res = cell(numel(row),numel(col));
for i = 1:numel(col)
    F = table2array(m2(i).anova([2,4,5],4));
    p = table2array(m2(i).anova([2,4,5],5));
    for j = 1:numel(row)
        % store F-/p-value
        if p(j) >= 0.001
            Res{j,i} = sprintf('F = %1.2f, p = %0.3f', F(j), p(j));
        else
            Res{j,i} = sprintf('F = %1.2f, p < 0.001', F(j));
        end;
    end;
end;

% save table
Res = [cell(1,1), col; row, Res];
xlswrite('Table_2.xls', Res);

% within-subject ANOVAs
figure('Name', 'FADE and SAMe ANOVAs', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

for i = 1:numel(m2)
    
    subplot(numel(m2),2,(i-1)*2+1); hold on;
    bar(m2(i).b_est(1:size(c2,2)), 'b');
    errorbar([1:size(c2,2)], m2(i).b_est(1:size(c2,2)), m2(i).CIs, '.k', 'LineWidth', 2, 'CapSize', 10);
    xlim([(1-0.5), (size(c2,2)+0.5)]);
    ylim([min([0, min(m2(i).b_est(1:size(c2,2))-m2(i).CIs)])-(1/10)*range(m2(i).b_est(1:size(c2,2))), ...
          max([0, max(m2(i).b_est(1:size(c2,2))+m2(i).CIs)])+(1/10)*range(m2(i).b_est(1:size(c2,2)))]);
    set(gca,'Box','On');
    set(gca,'XTick',[1:size(c2,2)],'XTickLabel',l2);
    xlabel('group', 'FontSize', 12);
    ylabel('estimate', 'FontSize', 12);
    if i == 1, title('novelty', 'FontSize', 16); end;
    if i == 2, title('memory', 'FontSize', 16); end;
    
    subplot(numel(m2),2,(i-1)*2+2);
    l = m2(i).anova.Properties.RowNames([1,2,4,5]);
    p = table2array(m2(i).anova([1,2,4,5],5));
    bar(p, 'r');
    axis([(1-0.5), (numel(p)+0.5), 0, 1]);
    set(gca,'XTick',[1:numel(p)],'XTickLabel',l);
    xlabel('effect', 'FontSize', 12);
    ylabel('p-value', 'FontSize', 12);
    for j = 1:numel(p)
        if p(j) <= 0.05/(numel(m2)*numel(p))
            text(j, 0.5, '***', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        elseif p(j) <= 0.05/numel(m2)
            text(j, 0.5, '**', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        elseif p(j)<=0.05
            text(j, 0.5, '*', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        end;
    end;
    
end;