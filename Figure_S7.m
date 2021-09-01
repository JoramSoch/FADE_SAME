% FADE-SAME: Figure S7

% clear
% close all

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set results files
FADE_files = {'../FADE_scores/FADE_SAME_scores_2021_01_11_ref_FADE_G1.xls';
              '../FADE_scores/FADE_SAME_scores_2021_01_11_ref_FADE_G2.xls';
              '../FADE_scores/FADE_SAME_scores_2021_01_11_ref_FADE_all.xls'};

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
age       = cell2mat(raw{1}(2:end,4));
FADE_data = FADE_data(age>=60,:,:);
num_subj  = size(FADE_data,1);

% rename FADE/SAMe scores
FADE_vars = raw{1}(1,FADE_inds);
for i = 1:num_vars
    FADE_vars{i}(strfind(FADE_vars{i},'_')) = '-';
end;


%%% Step 2: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% correlation analysis
n  = num_subj;
Y  = FADE_data(:,:,[1:2]);
X  = FADE_data(:,:,3);
V  = eye(n);
b1 = zeros(2,num_vars,2);
r1 = zeros(2,num_vars);
p1 = zeros(2,num_vars);
r2 = zeros(1,num_vars);
p2 = zeros(1,num_vars);
for i = 1:num_vars
    for j = 1:2
        b1(:,i,j) = ME_GLM(Y(:,i,j), [X(:,i), ones(n,1)], V);
        [r1(j,i), p1(j,i)] = corr(X(:,i), Y(:,i,j));
    end;
    M = squeeze(FADE_data(:,i,:));
    [r2(i), LB, UB, F, df1, df2, p2(i)] = ICC(M, 'C-1', 0.05, 0);
end;
clear M


%%% Step 3: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot correlation analyses
figure('Name', 'FADE and SAME scores', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

for i = 1:num_vars
    for j = 1:2
    
        subplot(sqrt(num_vars),sqrt(num_vars)*2,(i-1)*2+j); hold on;    
        Y = squeeze(FADE_data(:,i,[j,3]));
        y_min_max = [min(min(Y)), max(max(Y))];
        plot(y_min_max, y_min_max, '-k', 'LineWidth', 1);
        plot(y_min_max, y_min_max*b1(1,i,j)+b1(2,i,j), '--k', 'LineWidth', 1);
        plot(Y(:,2), Y(:,1), '.b', 'MarkerSize', 10);
        axis([y_min_max, y_min_max]);
        set(gca,'Box','On');
        xlabel('reference map from all subjects', 'FontSize', 12);
        ylabel(sprintf('reference map from CV group %d', j), 'FontSize', 12);
        if j == 1
            if i == 1, title('novelty contrast: FADE score', 'FontSize', 16, 'HorizontalAlignment', 'left'); end;
            if i == 2, title('novelty contrast: SAME score', 'FontSize', 16, 'HorizontalAlignment', 'left'); end;
            if i == 3, title('memory contrast: FADE score', 'FontSize', 16, 'HorizontalAlignment', 'left');  end;
            if i == 4, title('memory contrast: SAME score', 'FontSize', 16, 'HorizontalAlignment', 'left');  end;
        end;
        if j < 3
            if p1(j,i) >= 0.001
                text(y_min_max(1)+(1/20)*diff(y_min_max), y_min_max(2)-(1/20)*diff(y_min_max), ...
                     sprintf('r = %0.2f, p = %0.3f\ny = %0.2f x + %0.2f', r1(j,i), p1(j,i), b1(1,i,j), b1(2,i,j)), 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
            else
                text(y_min_max(1)+(1/20)*diff(y_min_max), y_min_max(2)-(1/20)*diff(y_min_max), ...
                     sprintf('r = %0.2f, p < 0.001\ny = %0.2f x + %0.2f', r1(j,i), b1(1,i,j), b1(2,i,j)), 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
            end;
        end;
        if j == 2
            if p2(i) >= 0.001
                text(y_min_max(2)-(1/20)*diff(y_min_max), y_min_max(1)+(1/20)*diff(y_min_max), ...
                     sprintf('ICC = %0.2f\np = %0.3f', r2(i), p2(i)), 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom');
            else
                text(y_min_max(2)-(1/20)*diff(y_min_max), y_min_max(1)+(1/20)*diff(y_min_max), ...
                     sprintf('ICC = %0.2f\np < 0.001', r2(i)), 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom');
            end;
        end;
            
        
    end;
end;