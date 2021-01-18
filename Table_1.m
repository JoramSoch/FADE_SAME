% FADE-SAME: Table 1

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

% between-subject design matrix
lvls1 = {'Verio', 'Skyra';
         'male',  'female';
         'young', 'old'}';
X1 = zeros(n,2^p);
l1 = cell(1,2^p);
for j1 = 1:2
    for j2 = 1:2
        for j3 = 1:2
            X1(:,(j1-1)*4+(j2-1)*2+(j3)) = double(X(:,1)==j1 & X(:,2)==j2 & X(:,3)==j3);
            l1{1,(j1-1)*4+(j2-1)*2+(j3)} = strcat(lvls1{j1,1}(1),'/',lvls1{j2,2}(1),'/',lvls1{j3,3}(1));
        end;
    end;
end;
c1 = eye(2^p);

% between-subject ANOVAs
for i = 1:numel(Y_vars)
    m1(i).Y_var = Y_vars(i);
    m1(i).X_var = X_vars;
    % estimate parameters
   [m1(i).b_est, m1(i).s2_est] = ME_GLM(Y(:,i), X1, eye(size(X1,1)));
    m1(i).Bcov = (X1'*X1)^(-1);
    m1(i).SEs  = sqrt(diag(c1' * (m1(i).s2_est*m1(i).Bcov) * c1));
    m1(i).CIs  = m1(i).SEs * norminv(mean([conf_lvl 1]), 0, 1);
    % perform ANOVA
    m1(i).model = fitlm(T, sprintf('%s ~ scanner*gender*age', Y_vars{i}));
    m1(i).anova = anova(m1(i).model);
   [m1(i).p, m1(i).t, m1(i).stats, m1(i).terms] = ...
        anovan(Y(:,i), {X(:,1), X(:,2), X(:,3)}, ...
               'varnames', X_vars, 'model', 'full', 'display', 'off');
end;


%%% Step 3: save results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate table (fitlm)
col = Y_labs;
row = m1(1).anova.Properties.RowNames(1:end-1);
Res = cell(numel(row),numel(col));
for i = 1:numel(col)
    F = table2array(m1(i).anova(1:end-1,4));
    p = table2array(m1(i).anova(1:end-1,5));
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
xlswrite('Table_1.xls', Res);

% between-subject ANOVAs
figure('Name', 'FADE and SAMe ANOVAs', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

for i = 1:numel(m1)
    
    subplot(numel(m1),2,(i-1)*2+1); hold on;
    bar(m1(i).b_est, 'b');
    errorbar([1:size(m1(i).b_est,1)], m1(i).b_est, m1(i).CIs, '.k', 'LineWidth', 2, 'CapSize', 10);
    xlim([(1-0.5), (size(m1(i).b_est,1)+0.5)]);
    ylim([min([0, min(m1(i).b_est-m1(i).CIs)])-(1/10)*range(m1(i).b_est), ...
          max([0, max(m1(i).b_est+m1(i).CIs)])+(1/10)*range(m1(i).b_est)]);
    set(gca,'Box','On');
    set(gca,'XTick',[1:size(m1(i).b_est,1)],'XTickLabel',l1);
    xlabel('group', 'FontSize', 12);
    ylabel('estimate', 'FontSize', 12);
    title(Y_labs{i}, 'FontSize', 16);
    
    subplot(numel(m1),2,(i-1)*2+2);
    l = cellstr(num2str(m1(1).terms));
    p = table2array(m1(i).anova(1:end-1,5));
    bar(p, 'r');
    axis([(1-0.5), (numel(p)+0.5), 0, 1]);
    set(gca,'XTick',[1:numel(p)],'XTickLabel',l);
    xlabel('effect', 'FontSize', 12);
    ylabel('p-value', 'FontSize', 12);
    for j = 1:numel(p)
        if p(j) <= 0.05/(numel(m1)*numel(p))
            text(j, 0.5, '***', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        elseif p(j) <= 0.05/numel(m1)
            text(j, 0.5, '**', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        elseif p(j)<=0.05
            text(j, 0.5, '*', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        end;
    end;
    
end;