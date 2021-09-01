% FADE-SAME: Table S5

% clear
% close all

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set results files
FADE_file = '../FADE_scores/FADE_SAME_scores_2021_01_11_FADE_only.xls';
vols_file = 'covariates/HC_volumes_FADE.xls';

% load results files
[num, txt, raw1] = xlsread(FADE_file);
[num, txt, raw2] = xlsread(vols_file);
clear num txt

% get results data
FADE_data = raw1(2:end,:);
vols_data = raw2(2:end,:);

% extract FADE/SAME scores
FADE_inds = 5+[1:4];
FADE_SAME = cell2mat(FADE_data(:,FADE_inds));
FADE_vars = raw1(1,FADE_inds);
num_subj  = size(FADE_SAME,1);
age       = cell2mat(FADE_data(:,4));

% get subject groups
AiA_inds = zeros(num_subj,1);
for i = 1:num_subj
    if strncmp(raw1{1+i,1},'subA',4)
        if age(i) < 50
            AiA_inds(i) = 1;    % young AiA
        elseif age(i) < 60
            AiA_inds(i) = 3;    % middle-aged AiA
        else
            AiA_inds(i) = 2;    % older AiA
        end;
    else
        AiA_inds(i) = 4;        % yFADE
    end;
end;

% extract covariate values
FADE_SAME = FADE_SAME(AiA_inds<3,:);
age       = age(AiA_inds<3);
age_grp   = 1*(age<50) + 2*(age>=60);
GMD       = cell2mat(vols_data(AiA_inds<3,2:11));


%%% Step 2: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compile data
Y = FADE_SAME;
X = GMD;

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
p_lim = 0.001;
p_thr = [0.05, 0.05/p, 0.05/(p*v)];
x_lab = {'Tail-left',  'Sub-left',  'CA1-left',  'CA3-left',  'CA4-left', ...
         'Tail-right', 'Sub-right', 'CA1-right', 'CA3-right', 'CA4-right'};

% collect correlations
row_hdr = {'HC', 'Subfield', 'Cohort', ...
           'novelty contrast: FADE score', 'novelty contrast: SAME score', ...
           'memory contrast: FADE score',  'memory contrast: SAME score'};
col_hdr = cell(p*2, 3);
results = cell(p*2, v);
for i = 1:v
    for j = 1:p
        % edit column header
        if i == 1
            hem = x_lab{j}(strfind(x_lab{j},'-')+1:end);
            reg = x_lab{j}(1:strfind(x_lab{j},'-')-1);
            if j == 1 || j == 6
                col_hdr{(j-1)*2+1,1} = hem;
            end;
            col_hdr{(j-1)*2+1,2} = reg;
            col_hdr{(j-1)*2+1,3} = 'young';
            col_hdr{(j-1)*2+2,3} = 'older';
        end;
        % store correlations
        if p1(j,i) < p_lim
            txt1 = sprintf('r = %0.2f, p < 0.001', r1(j,i));
        else
            txt1 = sprintf('r = %0.2f, p = %0.3f', r1(j,i), p1(j,i));
        end;
        for k = 1:numel(p_thr)
            if p1(j,i) < p_thr(k), txt1 = sprintf('%s*', txt1); end;
        end;
        if p2(j,i) < p_lim
            txt2 = sprintf('r = %0.2f, p < 0.001', r2(j,i));
        else
            txt2 = sprintf('r = %0.2f, p = %0.3f', r2(j,i), p2(j,i));
        end;
        for k = 1:numel(p_thr)
            if p2(j,i) < p_thr(k), txt2 = sprintf('%s*', txt2); end;
        end;
        results{(j-1)*2+1,i} = txt1;
        results{(j-1)*2+2,i} = txt2;
    end;
end;

% save correlations
res_tab = [row_hdr; col_hdr, results];
xlswrite('Table_S5.xls', res_tab);