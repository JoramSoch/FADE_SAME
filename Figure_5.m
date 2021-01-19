% FADE-SAME: Figure 5

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

% get subject groups
AiA_inds = zeros(num_subj,1);
for i = 1:num_subj
    if strncmp(raw{1}{1+i,1},'subA',4)
        if raw{1}{1+i,4} < 50
            AiA_inds(i) = 1;    % young AiA
        elseif raw{1}{1+i,4} < 60
            AiA_inds(i) = 3;    % middle-aged AiA
        else
            AiA_inds(i) = 2;    % older AiA
        end;
    else
        AiA_inds(i) = 4;        % yFADE
    end;
end;
num_subj = [sum(AiA_inds==1), sum(AiA_inds==2), sum(AiA_inds==3), sum(AiA_inds==4)];

% rename FADE/SAMe scores
FADE_vars = raw{1}(1,FADE_inds);
for i = 1:numel(FADE_vars)
    FADE_vars{i}(strfind(FADE_vars{i},'_')) = '-';
end;


%%% Step 2: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% two-sample t-tests
n1 = num_subj(1);               % AiA subjects (young)
n2 = num_subj(4);               % yFADE subjects (young)
Y  = [FADE_data(AiA_inds==1,:,1); FADE_data(AiA_inds==4,:,2)];
X  = blkdiag(ones(n1,1), ones(n2,1));
V  = eye(n1+n2);
b  = ME_GLM(Y, X, V);
dp = (b(1,:)-b(2,:))./sqrt(((n1-1)*var(Y(1:n1,:))+(n2-1)*var(Y((n1+1):end,:)))./(n1+n2-1));
se = [std(Y(1:n1,:))./sqrt(n1); std(Y((n1+1):end,:))./sqrt(n2)];
p  = zeros(1,size(Y,2));
[h, p_alt, stats] = ME_GLM_con(Y, X, V, [+1 -1], 'F', 0.05);
for i = 1:numel(FADE_vars)
    [h, p(i), ci, stats] = ttest2(Y(1:n1,i), Y((n1+1):end,i));
end;


%%% Step 3: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot two-sample t-tests
figure('Name', 'FADE and SAMe scores', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

for i = 1:numel(FADE_vars)
    subplot(1,numel(FADE_vars),i); hold on;
    bar(1, b(1,i), 'r');
    bar(2, b(2,i), 'm');
    errorbar([1:2], b(:,i)', se(:,i)', '.k', 'LineWidth', 2, 'CapSize', 15);
    xlim([(1-0.5), (2+0.5)]);
    ylim([(11/10)*min(min(b-se)), -(1/10)*min(min(b-se))]);
    set(gca,'Box','On');
    set(gca,'XTick',[1:2],'XTickLabel',{'young AiA', 'yFADE'});
    xlabel('subject group', 'FontSize', 12);
    if i == 1, ylabel('parameter estimates and standard errors', 'FontSize', 12); end;
    if i == 1, title(sprintf('novelty contrast:\nFADE score'), 'FontSize', 16);   end;
    if i == 2, title(sprintf('novelty contrast:\nSAME score'), 'FontSize', 16);   end;
    if i == 3, title(sprintf('memory contrast:\nFADE score'), 'FontSize', 16);    end;
    if i == 4, title(sprintf('memory contrast:\nSAME score'), 'FontSize', 16);    end;
    if p(i) >= 0.001
        text(1.5, -(1/20)*min(min(b-se)), sprintf('d'' = %0.2f, p = %0.3f', dp(i), p(i)), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    else
        text(1.5, -(1/20)*min(min(b-se)), sprintf('d'' = %0.2f, p < 0.001', dp(i)), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    end;
end;