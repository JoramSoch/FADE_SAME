% FADE-SAME: Figure 4

% clear
% close all

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set results files
FADE_file = '../FADE_scores/FADE_SAME_scores_2021_01_11_FADE_yFADE.xls';
covs_file = 'covariates/covariates.xls';

% load results files
[num, txt, raw1] = xlsread(FADE_file);
[num, txt, raw2] = xlsread(covs_file);
clear num txt

% get results data
FADE_data = raw1(2:end,:);
covs_data = raw2(2:end,:);

% extract FADE/SAME scores
FADE_inds = 5+[1:4];
FADE_SAME = cell2mat(FADE_data(:,FADE_inds));
FADE_vars = raw1(1,FADE_inds);
num_subj  = size(FADE_SAME,1);
age       = cell2mat(FADE_data(:,4));
age_grp   = 1*(age<50) + 2*(age>=60);

% extract covariate values
bhvr = cell2mat(covs_data(:,3));
GMD  = cell2mat(covs_data(:,[4 5]));


%%% Step 2: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compile data
Y = FADE_SAME;
X = [age, bhvr, GMD];

% get numbers
v  = size(Y,2);
p  = size(X,2);

% prellocate data
x0 = cell(p,v);
x1 = cell(p,v);
x2 = cell(p,v);
y0 = cell(p,v);
y1 = cell(p,v);
y2 = cell(p,v);
mn0= cell(p,v);
mn1= cell(p,v);
mn2= cell(p,v);

% preallocate results
r0 = zeros(p,v);
r1 = zeros(p,v);
r2 = zeros(p,v);
p0 = zeros(p,v);
p1 = zeros(p,v);
p2 = zeros(p,v);

% perform analyses
for i = 1:v
    for j = 1:p
        
        % correlation for all subjects
        y = Y(age_grp>0,i);
        x = X(age_grp>0,j);
        y0{j,i} = y(~isnan(x));
        x0{j,i} = x(~isnan(x));
        [r0(j,i), p0(j,i)] = corr(y0{j,i}, x0{j,i});
        mn0{j,i} = polyfit(x0{j,i}, y0{j,i}, 1);
        
        % correlation within young subjects
        y = Y(age_grp==1,i);
        x = X(age_grp==1,j);
        y1{j,i} = y(~isnan(x));
        x1{j,i} = x(~isnan(x));
        [r1(j,i), p1(j,i)] = corr(y1{j,i}, x1{j,i});
        mn1{j,i} = polyfit(x1{j,i}, y1{j,i}, 1);
        
        % correlation within older subjects
        y = Y(age_grp==2,i);
        x = X(age_grp==2,j);
        y2{j,i} = y(~isnan(x));
        x2{j,i} = x(~isnan(x));
        [r2(j,i), p2(j,i)] = corr(y2{j,i}, x2{j,i});
        mn2{j,i} = polyfit(x2{j,i}, y2{j,i}, 1);
        
    end;
end;


%%% Step 3: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify significance
p_thr = [0.05, 0.01, 0.001];
x_lab = {'age [yrs]', 'A'' [AUC]', 'V_{HC-left} [mm^3]', 'V_{HC-right} [mm^3]'};

% plot correlations
figure('Name', 'FADE and SAME scores', 'Color', [1 1 1], 'Position', [50 50 1280 800]);

for i = 1:v
    for j = 1:p
        
        subplot(p,v,(j-1)*v+i); hold on;
        plot(x1{j,i}, y1{j,i}, '.r', 'Color', [1, (p1(j,i)>p_thr(1))*1/2, (p1(j,i)>p_thr(1))*1/2]);
        plot(x2{j,i}, y2{j,i}, '.b', 'Color', [(p2(j,i)>p_thr(1))*1/2, (p2(j,i)>p_thr(1))*1/2, 1]);
        if j == 1
            x_lim = [min(x0{j,i})-5, max(x0{j,i})+5];
            x1_lim= [min(x1{j,i})-10, max(x1{j,i})+10];
            x2_lim= [min(x2{j,i})-10, max(x2{j,i})+10];
            plot(x1_lim, mn1{j,i}(1)*x1_lim + mn1{j,i}(2), '-r', 'Color', [1, (p1(j,i)>p_thr(1))*1/2, (p1(j,i)>p_thr(1))*1/2]);
            plot(x2_lim, mn2{j,i}(1)*x2_lim + mn2{j,i}(2), '-b', 'Color', [(p2(j,i)>p_thr(1))*1/2, (p2(j,i)>p_thr(1))*1/2, 1]);
        elseif j == 2
            x_lim = [0.5, 1];
        else
            x_lim = [min(x0{j,i})-1/20*range(x0{j,i}), max(x0{j,i})+1/20*range(x0{j,i})];
        end;
        if j ~= 1
            plot(x_lim, mn1{j,i}(1)*x_lim + mn1{j,i}(2), '-r', 'Color', [1, (p1(j,i)>p_thr(1))*1/2, (p1(j,i)>p_thr(1))*1/2]);
            plot(x_lim, mn2{j,i}(1)*x_lim + mn2{j,i}(2), '-b', 'Color', [(p2(j,i)>p_thr(1))*1/2, (p2(j,i)>p_thr(1))*1/2, 1]);
        end;
        xlim(x_lim);
        ylim([min(y0{j,i})-1/20*range(y0{j,i}), max(y0{j,i})+1/20*range(y0{j,i})]);
        set(gca,'Box','On');
        xlabel(x_lab{j}, 'FontSize', 12);
        ylabel('score', 'FontSize', 12);
        if j == 1
            if i == 1, title(sprintf('novelty contrast: FADE score\n '), 'FontSize', 12); end;
            if i == 2, title(sprintf('novelty contrast: SAME score\n '), 'FontSize', 12); end;
            if i == 3, title(sprintf('memory contrast: FADE score\n '), 'FontSize', 12);  end;
            if i == 4, title(sprintf('memory contrast: SAME score\n '), 'FontSize', 12);  end;
        end;
        if p1(j,i) < p_thr(end)
            txt1 = sprintf('r = %0.2f, p < 0.001', r1(j,i));
        else
            txt1 = sprintf('r = %0.2f, p = %0.3f', r1(j,i), p1(j,i));
        end;
        if p2(j,i) < p_thr(end)
            txt2 = sprintf('r = %0.2f, p < 0.001', r2(j,i));
        else
            txt2 = sprintf('r = %0.2f, p = %0.3f', r2(j,i), p2(j,i));
        end;
        text(x_lim(1), max(y0{j,i})+1/20*range(y0{j,i}), txt1, ...
             'FontSize', 8, 'FontWeight', 'Bold', 'Color', [1, (p1(j,i)>p_thr(1))*1/2, (p1(j,i)>p_thr(1))*1/2], ...
             'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Left');
        text(x_lim(2), max(y0{j,i})+1/20*range(y0{j,i}), txt2, ...
             'FontSize', 8, 'FontWeight', 'Bold', 'Color', [(p2(j,i)>p_thr(1))*1/2, (p2(j,i)>p_thr(1))*1/2, 1], ...
             'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Right');
        
    end;
end;