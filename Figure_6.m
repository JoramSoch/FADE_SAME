% FADE-SAME: Figure 6

% clear
% close all

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set results files
FADE_files = {'../FADE_scores/FADE_SAME_scores_2021_01_11_ref_FADE.xls';
              '../FADE_scores/FADE_SAME_scores_2021_01_11_ref_yFADE.xls'};

% load results files
num_files = numel(FADE_files);
for j = 1:num_files
    [num, txt, raw{j}] = xlsread(FADE_files{j});
end;
clear num txt

% get results data
FADE_inds = 5+[1:4];
num_subj  = size(raw{1},1)-1;
num_vars  = numel(FADE_inds);
FADE_data = zeros(num_subj,num_vars,num_files);
for j = 1:num_files
    FADE_data(:,:,j) = cell2mat(raw{j}(2:end,FADE_inds));
end;

% restrict results data
age       = cell2mat(raw{j}(2:end,4));
FADE_data = FADE_data(age>=60,:,:);
num_subj  = size(FADE_data,1);

% rename FADE/SAMe scores
FADE_vars = raw{1}(1,FADE_inds);
for i = 1:numel(FADE_vars)
    FADE_vars{i}(strfind(FADE_vars{i},'_')) = '-';
end;


%%% Step 2: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% correlation analysis
n = num_subj;
Y = FADE_data(:,:,2);
X = FADE_data(:,:,1);
V = eye(n);
b = zeros(2,size(Y,2));
r = zeros(1,num_vars);
p = zeros(1,num_vars);
for i = 1:num_vars
    b(:,i) = ME_GLM(Y(:,i), [X(:,i), ones(n,1)], V);
    [r(i), p(i)] = corr(X(:,i), Y(:,i));
end;


%%% Step 3: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot correlation analyses
figure('Name', 'FADE and SAME scores', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

for i = 1:numel(FADE_vars)
    subplot(sqrt(num_vars),sqrt(num_vars),i); hold on;
    Y = [FADE_data(:,i,1), FADE_data(:,i,2)];
    y_min_max = [min(min(Y)), max(max(Y))];
    plot(y_min_max, y_min_max, '-k', 'LineWidth', 1);
    plot(y_min_max, y_min_max*b(1,i)+b(2,i), '--k', 'LineWidth', 1);
    plot(Y(:,1), Y(:,2), '.b', 'MarkerSize', 10);
    axis([y_min_max, y_min_max]);
    set(gca,'Box','On');
    xlabel('reference map from young AiA', 'FontSize', 12);
    ylabel('reference map from yFADE', 'FontSize', 12);
    if i == 1, title('novelty contrast: FADE score', 'FontSize', 16); end;
    if i == 2, title('novelty contrast: SAME score', 'FontSize', 16); end;
    if i == 3, title('memory contrast: FADE score', 'FontSize', 16);  end;
    if i == 4, title('memory contrast: SAME score', 'FontSize', 16);  end;
    if p(i) >= 0.001
        text(mean(y_min_max), max(max(Y))-(1/20)*diff(y_min_max), ...
             sprintf('r = %0.2f, p = %0.3f, y = %0.2f x + %0.2f', r(i), p(i), b(1,i), b(2,i)), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    else
        text(mean(y_min_max), max(max(Y))-(1/20)*diff(y_min_max), ...
             sprintf('r = %0.2f, p < 0.001, y = %0.2f x + %0.2f', r(i), b(1,i), b(2,i)), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    end;
end;