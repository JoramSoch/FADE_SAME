% FADE-SAME: Figure S5

% clear
% close all

%%% Step 1: load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set results files
FADE_file = '../FADE_scores/FADE_SAME_scores_2021_01_11_FADE_only.xls';

% load results files
[num, txt, raw] = xlsread(FADE_file);
FADE_data = raw(2:end,:);
clear num txt

% extract FADE/SAMe scores
FADE_inds = 5+[1:4];
x = cell2mat(FADE_data(:,4));
Y = cell2mat(FADE_data(:,FADE_inds));
FADE_labs = raw(1,FADE_inds);
for i = 1:numel(FADE_labs)
    FADE_labs{i}(strfind(FADE_labs{i},'_')) = ':';
end;


%%% Step 2: analyze data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get numbers
n  = size(Y,1);
v  = size(Y,2);
n1 = sum(x< 50);
n2 = sum(x>=50);

% prepare smoothing
xw  = 32;
xp  = [(min(x)-xw/4):(max(x)+xw/4)]';
xb  = [(min(x)-0.5), (35+0.5), 50, (60-0.5), (max(x)+0.5)];

% preallocate results
m   = numel(xp); 
mun = zeros(m,v);
s2n = zeros(m,v);

% for each score
for i = 1:v
    % for each value
    for j = 1:m
        yij = Y(x>=(xp(j)-xw/2) & x<=(xp(j)+xw/2), i);
        nij = numel(yij);
        mun(j,i) = mean(yij);
        s2n(j,i) = var(yij);
    end;
end;


%%% Step 3: visualize results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot smoothing
figure('Name', 'Smoothed FADE and SAMe scores', 'Color', [1 1 1], 'Position', [50 50 1280 720]);

for i = 1:v
    subplot(floor(sqrt(v)),ceil(v/floor(sqrt(v))),i); hold on;
    if i == 4
        plot(x, Y(:,i), '.b', 'MarkerSize', 10);
        plot(xp, mun(:,i), '-r', 'LineWidth', 2);
        plot(xp, mun(:,i)+s2n(:,i), ':r', 'LineWidth', 2);
    end;
    plot([min(xp)-1, max(xp)+1], [0, 0], '-k', 'LineWidth', 1);
    plot([0, 0], [min(xp)-1, max(xp)+1], '--k', 'LineWidth', 1);
    plot(x, Y(:,i), '.b', 'MarkerSize', 10);
    plot(xp, mun(:,i), '-r', 'LineWidth', 2);
    plot(xp, mun(:,i)+s2n(:,i), ':r', 'LineWidth', 2);
    plot(xp, mun(:,i)-s2n(:,i), ':r', 'LineWidth', 2);
    for k = 1:numel(xb)
        plot([xb(k), xb(k)], [-100, +100], '--k', 'LineWidth', 1);
    end;
    xlim([min(xp)-1, max(xp)+1]);
    ylim([median(Y(:,i))-3*std(Y(:,i)), median(Y(:,i))+3*std(Y(:,i))]);
    set(gca,'Box','On');
    xlabel('age [yrs]', 'FontSize', 12);
    ylabel('score', 'FontSize', 12);
    if i == 1, title('novelty contrast: FADE score', 'FontSize', 16); end;
    if i == 2, title('novelty contrast: SAME score', 'FontSize', 16); end;
    if i == 3, title('memory contrast: FADE score', 'FontSize', 16);  end;
    if i == 4, title('memory contrast: SAME score', 'FontSize', 16);  end;
    if i == 4, legend('observations', 'smooth mean', 'smooth variance', 'Location', 'SouthWest'); end;
    if i == 2
        text((18+35)/2, -1.9, 'young subjects', 'FontSize', 10, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        text((50+59)/2, -1.9, sprintf('middle-\naged'), 'FontSize', 10, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        text((60+80)/2, -1.9, 'older subjects', 'FontSize', 10, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
    end;
end;