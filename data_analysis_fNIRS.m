%% fNIRS data analysis
%% Required to run in MATLAB version R2017b
%% Required to load Homer3 processed data files
%% Programmed by Feng Xiao (updated on 2025.1.2)
clear all,
clc,
%% Parameter settings
samplingrate = 5.85;
Channels = 1:34;
MediCond = [2 3]; %2 for mindfulness; 3 for endpoint-focus meditation
Mindfulness_first_subjs = [1 3 7 9 11 13 15 17 19 29 31 33 35 37 39];
Endpoint_first_subjs = [4 6 8 10 12 14 16 18 20 21 22 23 24 25 26 27 28 30 32 34 36 38 40];
subj = [Mindfulness_first_subjs, Endpoint_first_subjs];
sclConc = 1e6; %convert Conc from Molar to uMolar
%% Data analysis
cd fNIRS_preprocessed_data\derivatives\homer\
interval1 = [0 390]; %in second, for calculating means
ses_mindfulness_hbo = [];
ses_endpoint_hbo = [];
for i_subj = subj
load(num2str(i_subj), '-mat')
channelnum=size(output.dc.dataTimeSeries,3)/34;
setlength=channelnum*34; 
t=output.dcAvg.time;
temp_chPrune = cell2mat(output.misc.mlActAuto);
chPrune = temp_chPrune(1:34,3); %1:valid channels; 0: bad channels
subj_mindfulness_hbo = [];
subj_endpoint_hbo = [];
HbO_mindfulness = [];
HbO_endpoint = [];
    for j_Cond = MediCond
        for k_Ch = Channels
            SigHbO = output.dcAvg.dataTimeSeries(:,(j_Cond-1)*setlength*3+(k_Ch-1)*3+1)*sclConc;
            intind = [find(t >= interval1(1),1,'first') find(t <= interval1(end),1,'last')]; %extract data from 20s before stimulus onset to 390s after it
            HbO_bas_mean = mean(SigHbO(1:intind(1)));
            HbO_bas_std = std(SigHbO(1:intind(1)));
            HbO_act = SigHbO(intind(1):intind(end));
            HbO = (HbO_act - HbO_bas_mean) / HbO_bas_std; %calculate the relative activation based on the baseline (z-transformation)
            if chPrune(k_Ch,1) == 1
                if ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==2 || ismember(i_subj,Endpoint_first_subjs)==1&&j_Cond==2
                    HbO_mindfulness = [HbO_mindfulness, HbO];
                elseif ismember(i_subj,Endpoint_first_subjs)==1&&j_Cond==3 || ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==3
                    HbO_endpoint = [HbO_endpoint, HbO];
                end
            else
                if ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==2 || ismember(i_subj,Endpoint_first_subjs)==1&&j_Cond==2
                    HbO_mindfulness = [HbO_mindfulness, zeros(size(HbO,1),1)];
                elseif ismember(i_subj,Endpoint_first_subjs)==1&&j_Cond==3 || ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==3
                    HbO_endpoint = [HbO_endpoint, zeros(size(HbO,1),1)];
                end
            end
        end
    end
    subj_mindfulness_hbo = [subj_mindfulness_hbo, HbO_mindfulness];
    subj_endpoint_hbo = [subj_endpoint_hbo, HbO_endpoint];
    
    ses_mindfulness_hbo = [ses_mindfulness_hbo, subj_mindfulness_hbo];
    ses_endpoint_hbo = [ses_endpoint_hbo, subj_endpoint_hbo]; %col: channels * subj
end

clear subj_mindfulness_hbo subj_endpoint_hbo HbO_mindfulness HbO_endpoint

Ses_mindfulness_hbo = [];
Ses_endpoint_hbo = [];
Ses_mindfulness_hbo_se = [];
Ses_endpoint_hbo_se = [];
Ses_channel_hbo_t = ones(size(ses_mindfulness_hbo, 1), 34);
Ses_channel_hbo_p = ones(size(ses_mindfulness_hbo, 1), 34);
for i = Channels
    i_temp_mindfulness_hbo = [];
    i_temp_endpoint_hbo = [];
    i_mean_mindfulness_hbo = [];
    i_mean_endpoint_hbo = [];
    i_se_mindfulness_hbo = [];
    i_se_endpoint_hbo = [];
    ch_temp_t = ones(size(ses_mindfulness_hbo, 1), 1);
    ch_temp_p = ones(size(ses_mindfulness_hbo, 1), 1);
    for j = 1:size(subj, 2)
        temp_mindfulness_hbo = ses_mindfulness_hbo(:, i+size(Channels, 2)*(j-1));
        temp_endpoint_hbo = ses_endpoint_hbo(:, i+size(Channels, 2)*(j-1));
        
        i_temp_mindfulness_hbo = [i_temp_mindfulness_hbo, temp_mindfulness_hbo];
        i_temp_endpoint_hbo = [i_temp_endpoint_hbo, temp_endpoint_hbo];
    end
    i_temp_mindfulness_hbo = i_temp_mindfulness_hbo(:, any(i_temp_mindfulness_hbo)); 
    i_temp_endpoint_hbo = i_temp_endpoint_hbo(:, any(i_temp_endpoint_hbo)); %delete column with all zeros
    
    i_mean_mindfulness_hbo = mean(i_temp_mindfulness_hbo, 2);%average subjects' data for each channel
    i_mean_endpoint_hbo = mean(i_temp_endpoint_hbo, 2); %average subjects' data for each channel
    
    i_se_mindfulness_hbo = std(i_temp_mindfulness_hbo, 0, 2)./sqrt(size(i_temp_mindfulness_hbo, 2)); %calculate standard error
    i_se_endpoint_hbo = std(i_temp_endpoint_hbo, 0, 2)./sqrt(size(i_temp_endpoint_hbo, 2)); %calculate standard error   
    for k = 1:size(ses_mindfulness_hbo,1)
        [h,p_hbo,ci,stats_hbo] = ttest(i_temp_mindfulness_hbo(k,:), i_temp_endpoint_hbo(k,:)); %paired ttest between two meditation for each time point in each channel
        ch_temp_hbo_t(k) = stats_hbo.tstat; %t>0 denote mindfulness > intertemporal
        ch_temp_hbo_p(k) = p_hbo;
    end
    Ses_mindfulness_hbo = [Ses_mindfulness_hbo, i_mean_mindfulness_hbo];
    Ses_endpoint_hbo = [Ses_endpoint_hbo, i_mean_endpoint_hbo];
    Ses_mindfulness_hbo_se = [Ses_mindfulness_hbo_se, i_se_mindfulness_hbo];
    Ses_endpoint_hbo_se = [Ses_endpoint_hbo_se, i_se_endpoint_hbo];
    Ses_channel_hbo_t(:, i) = ch_temp_hbo_t;
    Ses_channel_hbo_p(:, i) = ch_temp_hbo_p;
end

clear i_temp_mindfulness_hbo
clear i_temp_endpoint_hbo
clear i_mean_mindfulness_hbo 
clear i_mean_endpoint_hbo 
clear i_se_mindfulness_hbo 
clear i_se_endpoint_hbo
clear temp_mindfulness_hbo 
clear temp_endpoint_hbo 
clear ch_temp_hbo_t 
clear ch_temp_hbo_p 
clear h ci p_hbo stats_hbo 

tTest_Ses_mindfulness_hbo = zeros(5,34);
tTest_Ses_endpoint_hbo = zeros(5,34);
tTest_Ses_contrast_hbo = zeros(5,34);
for i = 1:34
    temp_col_m = Ses_mindfulness_hbo(:, i);
    [h,p_m,ci_m,stats_m] = ttest(temp_col_m); %one-sample ttest
    tTest_Ses_mindfulness_hbo(1, i) = stats_m.tstat;
    tTest_Ses_mindfulness_hbo(2, i) = p_m;
    tTest_Ses_mindfulness_hbo(3, i) = ci_m(1,1);
    tTest_Ses_mindfulness_hbo(4, i) = ci_m(2,1);
    if p_m < 0.001 %alpha level of .001
       tTest_Ses_mindfulness_hbo(5, i) = 1; %accept h1
    else
       tTest_Ses_mindfulness_hbo(5, i) = 0; %accept h0 
    end
    
    temp_col_i = Ses_endpoint_hbo(:, i);
    [h,p_i,ci_i,stats_i] = ttest(temp_col_i); %one-sample ttest
    tTest_Ses_endpoint_hbo(1, i) = stats_i.tstat;
    tTest_Ses_endpoint_hbo(2, i) = p_i;
    tTest_Ses_endpoint_hbo(3, i) = ci_i(1,1);
    tTest_Ses_endpoint_hbo(4, i) = ci_i(2,1);
    if p_i < 0.001 %alpha level of .001
       tTest_Ses_endpoint_hbo(5, i) = 1; %accept h1
    else
       tTest_Ses_endpoint_hbo(5, i) = 0; %accept h0 
    end

    [h,p_c,ci_c,stats_c] = ttest(temp_col_i, temp_col_m); %endpoint - mindfulness; paired ttest
    tTest_Ses_contrast_hbo(1,i) = stats_c.tstat; 
    tTest_Ses_contrast_hbo(2,i) = p_c;
    tTest_Ses_contrast_hbo(3,i) = ci_c(1,1);
    tTest_Ses_contrast_hbo(4,i) = ci_c(2,1); 
    if p_c < 0.001 %alpha level of .001
       tTest_Ses_contrast_hbo(5, i) = 1; %accept h1
    else
       tTest_Ses_contrast_hbo(5, i) = 0; %accept h0 
    end
end

clear temp_col_m temp_col_i
clear h p_m p_i p_c ci_m ci_i ci_c stats_m stats_i stats_c
%% Effect size calculation (Cohen's d)
hbo_mindfulness_ses = zeros(3,34);
hbo_endpoint_ses = zeros(3,34);

for i = 1:34
    temp_hbo_mean = mean(Ses_mindfulness_hbo(:,i));
    temp_hbo_se = std(Ses_mindfulness_hbo(:,i))/sqrt(size(Ses_mindfulness_hbo(:,i),1));
    hbo_mindfulness_ses(1,i) = temp_hbo_mean;
    hbo_mindfulness_ses(2,i) = temp_hbo_se;
    hbo_mindfulness_ses(3,i) = std(Ses_mindfulness_hbo(:,i));
    
    clear temp_hbo_mean temp_hbo_se
    
    temp_hbo_mean = mean(Ses_endpoint_hbo(:,i));
    temp_hbo_se = std(Ses_endpoint_hbo(:,i))/sqrt(size(Ses_endpoint_hbo(:,i),1));
    hbo_endpoint_ses(1,i) = temp_hbo_mean;
    hbo_endpoint_ses(2,i) = temp_hbo_se;
    hbo_endpoint_ses(3,i) = std(Ses_endpoint_hbo(:,i));
    
    clear temp_hbo_mean temp_hbo_se
end

efs_hbo_contrast = zeros(1,34); %overall
n_ses1 = size(Ses_endpoint_hbo,1);
for i = 1:34
    temp_mean_ses = hbo_endpoint_ses(1,i) - hbo_mindfulness_ses(1,i);  %Endpoint - Mindfulness
    temp_poolstd_ses = sqrt((hbo_endpoint_ses(3,i)^2 + hbo_mindfulness_ses(3,i)^2) / 2);
    efs_hbo_contrast(1,i) = temp_mean_ses/temp_poolstd_ses;
end

clear temp_mean_overall temp_mean_ses temp_poolstd_overall

efs_hbo_endpoint = zeros(1,34);
efs_hbo_mindfulness = zeros(1,34); 
for i = 1:34
    efs_hbo_endpoint(1,i) = hbo_endpoint_ses(1,i)/hbo_endpoint_ses(3,i);
    efs_hbo_mindfulness(1,i) = hbo_mindfulness_ses(1,i)/hbo_mindfulness_ses(3,i);
end
%% Visualization for overall activation contrast
cd ..\..\..
cd fNIRS_pics\
%% BA11L & BA46L: Endpoint-focus > Mindfulness
labels_pair1 = {'CH1', 'CH7'};
means_CH1 = [hbo_endpoint_ses(1,1), hbo_mindfulness_ses(1,1)];
means_CH7 = [hbo_endpoint_ses(1,7), hbo_mindfulness_ses(1,7)];
se_CH1 = [hbo_endpoint_ses(2,1), hbo_mindfulness_ses(2,1)];
se_CH7 = [hbo_endpoint_ses(2,7), hbo_mindfulness_ses(2,7)];
means_pair1 = [means_CH1; means_CH7];
se_pair1 = [se_CH1; se_CH7];

Pair1_plot = figure;
set(Pair1_plot, 'Units', 'inches', 'Position', [0, 0, 2, 1.5]);
bar_handle = bar(means_pair1, 'grouped');
hold on;
bar_handle(1).FaceColor = [0.8, 0.1, 0.1]; % Endpoint: red
bar_handle(2).FaceColor = [0.1, 0.1, 0.8]; % Mindfulness: blue
[ngroups, nbars] = size(means_pair1);
group_width = bar_handle(1).BarWidth;
x_positions = zeros(nbars, ngroups);
for i = 1:nbars
    x_positions(i, :) = bar_handle(i).XData + bar_handle(i).XOffset;
    errorbar(x_positions(i,:), means_pair1(:,i), se_pair1(:,i), 'k', 'Color', [0, 0, 0], 'linewidth', 0.3,...
        'linestyle', 'none', 'CapSize', 5);
end
ax = gca;
ax.XColor = 'k'; 
ax.YColor = 'k'; 
ax.TickLength = [0.01, 0.01];
set(gca, 'XTick', 1:ngroups, 'XTickLabel', labels_pair1, 'FontName', 'Times New Roman', 'FontSize', 7);
ylabel('HbO amplitude (uM)', 'FontName', 'Times New Roman', 'FontSize', 7, 'Color', 'k');
title('BA11L & BA46L: Endpoint-focus > Mindfulness', 'FontName', 'Times New Roman', 'Color', 'k', ...
      'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
ylim([-3, 3]);
grid off
box off
legend({'Endpoint-focus', 'Mindfulness'}, 'Location', 'northeast',...
       'FontName', 'Times New Roman', 'FontSize', 7);
hold off;
print(Pair1_plot,'HbO overall contrast Pair1','-dpdf','-r600')
%% BA9L & BA10L: Mindfulness > Endpoint-focus
labels_pair2 = {'CH12', 'CH16'};
means_CH12 = [hbo_endpoint_ses(1,12), hbo_mindfulness_ses(1,12)];
means_CH16 = [hbo_endpoint_ses(1,16), hbo_mindfulness_ses(1,16)];
se_CH12 = [hbo_endpoint_ses(2,12), hbo_mindfulness_ses(2,12)];
se_CH16 = [hbo_endpoint_ses(2,16), hbo_mindfulness_ses(2,16)];
means_pair2 = [means_CH12; means_CH16];
se_pair2 = [se_CH12; se_CH16];

Pair2_plot = figure;
set(Pair2_plot, 'Units', 'inches', 'Position', [0, 0, 2, 1.5]);
bar_handle = bar(means_pair2, 'grouped');
hold on;
bar_handle(1).FaceColor = [0.8, 0.1, 0.1]; % Endpoint: red
bar_handle(2).FaceColor = [0.1, 0.1, 0.8]; % Mindfulness: white
[ngroups, nbars] = size(means_pair2);
group_width = bar_handle(1).BarWidth;
x_positions = zeros(nbars, ngroups);
for i = 1:nbars
    x_positions(i, :) = bar_handle(i).XData + bar_handle(i).XOffset;
    errorbar(x_positions(i,:), means_pair2(:,i), se_pair2(:,i), 'k', 'Color', [0, 0, 0], 'linewidth', 0.3,...
        'linestyle', 'none', 'CapSize', 5);
end
ax = gca;
ax.XColor = 'k'; 
ax.YColor = 'k'; 
ax.TickLength = [0.01, 0.01];
set(gca, 'XTick', 1:ngroups, 'XTickLabel', labels_pair2, 'FontName', 'Times New Roman', 'FontSize', 7);
ylabel('HbO amplitude (uM)', 'FontName', 'Times New Roman', 'FontSize', 7, 'Color', 'k');
title('BA9L & BA10L: Mindfulness > Endpoint-focus', 'FontName', 'Times New Roman', 'Color', 'k', ...
      'FontSize', 7, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
ylim([-3, 3]);
grid off
box off
legend({'Endpoint-focus', 'Mindfulness'}, 'Location', 'northeast',...
       'FontName', 'Times New Roman', 'FontSize', 7);
hold off;
print(Pair2_plot,'HbO overall contrast Pair2','-dpdf','-r600')
%% Epoch-based activation detection and contrast
total_timepoints = 2268; 
points_per_segment = 87;
num_segments = floor(total_timepoints / points_per_segment);
num_channels = size(Ses_endpoint_hbo, 2);
%% Endpoint-focus meditation
segment_means_endpoint = zeros(num_segments, num_channels);
segment_se_endpoint = zeros(num_segments, num_channels);
segment_ttest_endpoint = zeros(num_segments, num_channels); 
segment_pvalue_endpoint = zeros(num_segments, num_channels);
segment_cohend_endpoint = zeros(num_segments, num_channels);
for ch = 1:num_channels
    channel_data = Ses_endpoint_hbo(:, ch); 
    for seg = 1:num_segments
        start_idx = (seg - 1) * points_per_segment + 1;
        end_idx = start_idx + points_per_segment - 1;
        segment_data = channel_data(start_idx:end_idx); 

        segment_means_endpoint(seg, ch) = mean(segment_data); % mean
        segment_se_endpoint(seg, ch) = std(segment_data) ./ sqrt(points_per_segment); % SE

        [~, p, ~, stats] = ttest(segment_data, 0); % one-sample ttest
        segment_ttest_endpoint(seg, ch) = stats.tstat; % for t value
        segment_pvalue_endpoint(seg, ch) = p; % for p value
        
        segment_cohend_endpoint(seg, ch) = mean(segment_data) / std(segment_data); %cohen's d
    end
end

[num_epochs, num_channels] = size(segment_means_endpoint);
time_epochs = (1:num_epochs)';
beta_endpoint = zeros(1, num_channels); % beta
p_endpoint_reg = zeros(1, num_channels); % slope p
r_squared_endpoint = zeros(1, num_channels); % R_square
for ch = 1:num_channels
    activation_values = segment_means_endpoint(:, ch); 
    
    mdl = fitlm(time_epochs, activation_values); % linear regression
    
    beta_endpoint(ch) = mdl.Coefficients.Estimate(2); 
    p_endpoint_reg(ch) = mdl.Coefficients.pValue(2);   
    r_squared_endpoint(ch) = mdl.Rsquared.Ordinary;   
end
%% Mindfulness meditation
segment_means_mindfulness = zeros(num_segments, num_channels);
segment_se_mindfulness = zeros(num_segments, num_channels);
segment_ttest_mindfulness = zeros(num_segments, num_channels); 
segment_pvalue_mindfulness = zeros(num_segments, num_channels);
segment_cohend_mindfulness = zeros(num_segments, num_channels);
for ch = 1:num_channels
    channel_data = Ses_mindfulness_hbo(:, ch); 
    for seg = 1:num_segments
        start_idx = (seg - 1) * points_per_segment + 1;
        end_idx = start_idx + points_per_segment - 1;
        segment_data = channel_data(start_idx:end_idx); 

        segment_means_mindfulness(seg, ch) = mean(segment_data); % mean
        segment_se_mindfulness(seg, ch) = std(segment_data) ./ sqrt(points_per_segment); % SE

        [~, p, ~, stats] = ttest(segment_data, 0); % one-sample ttest
        segment_ttest_mindfulness(seg, ch) = stats.tstat; % for t value
        segment_pvalue_mindfulness(seg, ch) = p; % for p value
        
        segment_cohend_mindfulness(seg, ch) = mean(segment_data) / std(segment_data); %cohen's d
    end
end

[num_epochs, num_channels] = size(segment_means_mindfulness);
time_epochs = (1:num_epochs)';
beta_mindfulness = zeros(1, num_channels); % beta
p_mindfulness_reg = zeros(1, num_channels); % slope p
r_squared_mindfulness = zeros(1, num_channels); % R_square
for ch = 1:num_channels
    activation_values = segment_means_mindfulness(:, ch); 
    
    mdl = fitlm(time_epochs, activation_values); % linear regression
    
    beta_mindfulness(ch) = mdl.Coefficients.Estimate(2); 
    p_mindfulness_reg(ch) = mdl.Coefficients.pValue(2);   
    r_squared_mindfulness(ch) = mdl.Rsquared.Ordinary;   
end
%% Meditation contrast
segment_ttest_comparison = zeros(num_segments, num_channels); 
segment_pvalue_comparison = zeros(num_segments, num_channels); 
segment_cohend_comparison = zeros(num_segments, num_channels); 
for ch = 1:num_channels
    for seg = 1:num_segments
        
        data_endpoint = Ses_endpoint_hbo((seg-1)*points_per_segment+1:seg*points_per_segment, ch);
        data_mindfulness  = Ses_mindfulness_hbo((seg-1)*points_per_segment+1:seg*points_per_segment, ch);

        [~, p, ~, stats] = ttest(data_endpoint, data_mindfulness); % paired ttests (endpoint - mindfulness)

        mean_diff = mean(data_endpoint - data_mindfulness); % mean
        pooled_sd = std(data_endpoint - data_mindfulness); %pooled sd
        cohen_d = mean_diff / pooled_sd;

        segment_ttest_comparison(seg, ch) = stats.tstat; % for t value
        segment_pvalue_comparison(seg, ch) = p; % for p value
        segment_cohend_comparison(seg, ch) = cohen_d; % for Cohen's d
    end
end
%% Visualization for all channels
output_folder = fullfile(pwd, 'Epoch_figs');
invalid_channels = [10, 13, 21, 23, 25, 27, 30, 31, 33, 34];
valid_channels = setdiff(1:num_channels, invalid_channels); 

color_endpoint = [0.8, 0.1, 0.1]; % red (Endpoint)
color_mindfulness = [0.1, 0.1, 0.8];   % blue (Mindfulness)

for ch = valid_channels
    mean_endpoint = segment_means_endpoint(:, ch); 
    se_endpoint = segment_se_endpoint(:, ch);   
    mean_mindfulness = segment_means_mindfulness(:, ch); 
    se_mindfulness = segment_se_mindfulness(:, ch);

    figure('Units', 'inches', 'Position', [0, 0, 4, 3]); 
    hold on;

    errorbar(1:num_segments, mean_endpoint, se_endpoint, '-o', ...
        'Color', color_endpoint, 'MarkerFaceColor', color_endpoint, ...
        'LineWidth', 0.3, 'MarkerSize', 5, 'CapSize', 7);
    errorbar(1:num_segments, mean_mindfulness, se_mindfulness, '-o', ...
        'Color', color_mindfulness, 'MarkerFaceColor', color_mindfulness, ...
        'LineWidth', 0.3, 'MarkerSize', 5, 'CapSize', 7);

    xlabel('Epoch Number', 'FontName', 'Times New Roman', 'FontSize', 7);
    ylabel('HbO Activation', 'FontName', 'Times New Roman', 'FontSize', 7);
    title(['CH ', num2str(ch), ': HbO activation between two meditation practices'], ...
        'FontName', 'Times New Roman', 'FontSize', 7, 'FontWeight', 'bold');

    legend({'Endpoint-focus', 'Mindfulness'}, 'Location', 'northeast', ...
        'FontName', 'Times New Roman', 'FontSize', 7);
    grid off;

    ylim([min([mean_endpoint - se_endpoint; mean_mindfulness - se_mindfulness]) - 0.5, ...
          max([mean_endpoint + se_endpoint; mean_mindfulness + se_mindfulness]) + 0.5]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 7);
    box off;
    hold off;

    save_path = fullfile(output_folder, ['CH', num2str(ch), '.pdf']);
    saveas(gcf, save_path);

    fprintf('Channel %d saved as CH%d.pdf\n', ch, ch);
    close; 
end
disp('All valid channel plots have been saved.');
%% Visualization for selected channels 
output_folder = fullfile(pwd, 'Epoch_figs_selected');
color_endpoint = [0.8, 0.1, 0.1]; % red (Endpoint)
color_mindfulness = [0.1, 0.1, 0.8];   % blue (Mindfulness)
color_inter_fit = [0.5, 0, 0];         % deep red£¨regression line for endpoint)
color_mind_fit = [0, 0, 0.5];          % deep blue£¨regression line for mindfulness£©

epoch_duration = 15; 
time_axis = (1:num_segments)' * epoch_duration - (epoch_duration / 2); % median temporal point for each epoch
channels_to_plot = [7, 11, 18, 22];

for ch = channels_to_plot
    mean_endpoint = segment_means_endpoint(:, ch); 
    se_endpoint = segment_se_endpoint(:, ch);   
    mean_mindfulness = segment_means_mindfulness(:, ch); 
    se_mindfulness = segment_se_mindfulness(:, ch);

    p_inter = polyfit(time_axis, mean_endpoint, 1);
    fit_intertemporal = polyval(p_inter, time_axis);

    p_mind = polyfit(time_axis, mean_mindfulness, 1);
    fit_mindfulness = polyval(p_mind, time_axis);

    figure('Units', 'inches', 'Position', [0, 0, 2, 1.5]); 
    hold on;

    fill([60, 390, 390, 60], [-10, -10, 10, 10], [0.9, 0.9, 0.9], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.7);

    errorbar(time_axis, mean_endpoint, se_endpoint, '-o', ...
        'Color', color_endpoint, 'MarkerFaceColor', color_endpoint, ...
        'LineWidth', 0.3, 'MarkerSize', 3, 'CapSize', 5);
    errorbar(time_axis, mean_mindfulness, se_mindfulness, '-o', ...
        'Color', color_mindfulness, 'MarkerFaceColor', color_mindfulness, ...
        'LineWidth', 0.3, 'MarkerSize', 3, 'CapSize', 5);

    plot(time_axis, fit_intertemporal, '-', 'Color', color_inter_fit, 'LineWidth', 1);
    plot(time_axis, fit_mindfulness, '-', 'Color', color_mind_fit, 'LineWidth', 1);

    plot([time_axis(1), time_axis(end)], [0, 0], '--k', 'LineWidth', 0.5);

    xlabel('Time (s)', 'FontName', 'Times New Roman', 'FontSize', 7);
    ylabel('HbO amplitude (uM)', 'FontName', 'Times New Roman', 'FontSize', 7);
    title(['CH ', num2str(ch), ': HbO activation dynamics'], ...
        'FontName', 'Times New Roman', 'FontSize', 7, 'FontWeight', 'bold');
    ylim([-10, 10]);
    xlim([0, 390]);

    set(gca, 'FontName', 'Times New Roman', 'FontSize', 7);
    box off;
    hold off;

    save_path = fullfile(output_folder, ['CH', num2str(ch), '.pdf']);
    saveas(gcf, save_path);
    fprintf('Channel %d saved as %s\n', ch, save_path);
    close;
end

disp('Selected channel plots with regression lines have been saved to the "Epoch_figs_selected" folder.');
%% Resting-state channel connectivity analysis
num_channels = size(Ses_endpoint_hbo, 2);
r_matrix_endpoint = zeros(num_channels, num_channels);
p_matrix_endpoint = zeros(num_channels, num_channels);
r_matrix_mindfulness = zeros(num_channels, num_channels);
p_matrix_mindfulness = zeros(num_channels, num_channels);
%% Endpoint-focus meditation
for i = 1:num_channels
    for j = i:num_channels
        [r, p] = corr(Ses_endpoint_hbo(:, i), Ses_endpoint_hbo(:, j));
        r_matrix_endpoint(i, j) = r; 
        r_matrix_endpoint(j, i) = r; 
        p_matrix_endpoint(i, j) = p; 
        p_matrix_endpoint(j, i) = p; 
    end
end
%% Mindfulness meditation
for i = 1:num_channels
    for j = i:num_channels
        [r, p] = corr(Ses_mindfulness_hbo(:, i), Ses_mindfulness_hbo(:, j));
        r_matrix_mindfulness(i, j) = r; 
        r_matrix_mindfulness(j, i) = r; 
        p_matrix_mindfulness(i, j) = p; 
        p_matrix_mindfulness(j, i) = p; 
    end
end
%% Fisher-z transformation
z_matrix_endpoint = 0.5 * log((1 + r_matrix_endpoint) ./ (1 - r_matrix_endpoint));
z_matrix_mindfulness = 0.5 * log((1 + r_matrix_mindfulness) ./ (1 - r_matrix_mindfulness));
z_matrix_endpoint(logical(eye(size(z_matrix_endpoint)))) = NaN;
z_matrix_mindfulness(logical(eye(size(z_matrix_mindfulness)))) = NaN;
%% Calculate mean Fisher-z values
mean_z_endpoint = nanmean(z_matrix_endpoint(:)); % Endpoint-focus: 0.54
mean_z_mindfulness = nanmean(z_matrix_mindfulness(:)); % Mindfulness: 0.37
%% Between-meditation analysis: endpoint - mindfulness
z_diff = z_matrix_endpoint - z_matrix_mindfulness;
z_diff(logical(eye(size(z_diff)))) = NaN;
mean_diff = nanmean(z_diff(:));
std_diff = nanstd(z_diff(:));
n_connections = num_channels * (num_channels - 1) / 2; 
se_diff = std_diff / sqrt(n_connections); 
t_stat = mean_diff / se_diff; % t(560)=7.51
df = n_connections - 1; 
p_value = 2 * (1 - tcdf(abs(t_stat), df)); % p<.001
cohens_d = mean_diff / std_diff; % 0.32
%% FC visualization
invalid_channels = [10, 13, 21, 23, 25, 27, 30, 31, 33, 34];
threshold = 0.5; %cut-off for high FC
%% Endpoint-focus meditation
z_matrix_endpoint_tri = tril(z_matrix_endpoint, -1);
z_matrix_endpoint_tri(isnan(z_matrix_endpoint_tri)) = NaN; 
num_channels = size(z_matrix_endpoint, 1);
all_channels = 1:num_channels;
valid_channels = setdiff(all_channels, invalid_channels); % delete pruned channels
z_matrix_filtered = z_matrix_endpoint_tri(valid_channels, valid_channels);

FC_endpoint = figure;
set(FC_endpoint, 'Units', 'inches', 'Position', [0, 0, 4, 4]);
imagesc(z_matrix_filtered);
colormap(autumn);
colorbar;
caxis([-1, 1]);
axis square; 
num_valid_channels = length(valid_channels); 
ax = gca;
ax.XColor = 'k'; 
ax.YColor = 'k'; 
ax.TickLength = [0.01, 0.01];
xticks(1:num_valid_channels);
yticks(1:num_valid_channels);
xticklabels(valid_channels);
yticklabels(valid_channels);
xlabel('Channels', 'FontName', 'Times New Roman', 'FontSize', 7, 'Color', 'k');
ylabel('Channels', 'FontName', 'Times New Roman', 'FontSize', 7, 'Color', 'k');
title('Functional connectivity: Endpoint-focus meditation', 'FontName', 'Times New Roman', 'FontSize', 7,...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');
grid off;
box off;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 7);
print(FC_endpoint,'FC endpoint','-dpdf','-r600')

z_valid_endpoint = z_matrix_filtered + z_matrix_filtered'; 
z_valid_endpoint(abs(z_valid_endpoint) < threshold) = 0;
filename = 'FC_endpoint.edge';
dlmwrite(filename, z_valid_endpoint, 'delimiter', '\t'); %save .edge file for FC endpoint
%% Mindfulness meditation
z_matrix_mindfulness_tri = tril(z_matrix_mindfulness, -1);
z_matrix_mindfulness_tri(isnan(z_matrix_mindfulness_tri)) = NaN; 
num_channels = size(z_matrix_mindfulness, 1);
all_channels = 1:num_channels;
valid_channels = setdiff(all_channels, invalid_channels); % delete pruned channels
z_matrix_filtered = z_matrix_mindfulness_tri(valid_channels, valid_channels);

FC_mindfulness = figure;
set(FC_mindfulness, 'Units', 'inches', 'Position', [0, 0, 4, 4]);
imagesc(z_matrix_filtered);
colormap(autumn);
colorbar;
caxis([-1, 1]);
axis square; 
num_valid_channels = length(valid_channels); 
ax = gca;
ax.XColor = 'k'; 
ax.YColor = 'k'; 
ax.TickLength = [0.01, 0.01];
xticks(1:num_valid_channels);
yticks(1:num_valid_channels);
xticklabels(valid_channels);
yticklabels(valid_channels);
xlabel('Channels', 'FontName', 'Times New Roman', 'FontSize', 7, 'Color', 'k');
ylabel('Channels', 'FontName', 'Times New Roman', 'FontSize', 7, 'Color', 'k');
title('Functional connectivity: Mindfulness meditation', 'FontName', 'Times New Roman', 'FontSize', 7,...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');
grid off;
box off;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 7);
print(FC_mindfulness,'FC mindfulness','-dpdf','-r600')

z_valid_mindfulness = z_matrix_filtered + z_matrix_filtered'; 
z_valid_mindfulness(abs(z_valid_mindfulness) < threshold) = 0;
filename = 'FC_mindfulness.edge';
dlmwrite(filename, z_valid_mindfulness, 'delimiter', '\t'); %save .edge file for FC mindfulness
%% MNI node file generation
[num_data, txt_data, raw_data] = xlsread('FC table.xlsx', 'mni');
mni_data = raw_data;
mni_coords = cell2mat(mni_data(:, 1:3));
channel_labels = cell2mat(mni_data(:, 4));
color = ones(18, 1); %1 for prefrontal lobe's channels
color(19:24) = 2; %2 for temporal lobe's channels
degree_endpoint = sum(abs(z_valid_endpoint) > 0, 2); %calculate the FC counts for each channel
degree_mindfulness = sum(abs(z_valid_mindfulness) > 0, 2); %calculate the FC counts for each channel

endpoint_node = [mni_coords, color, degree_endpoint, channel_labels];
mindfulness_node = [mni_coords, color, degree_mindfulness, channel_labels];

filename_node = 'nodes_endpoint.node';
dlmwrite(filename_node, endpoint_node, 'delimiter', '\t');
filename_node = 'nodes_mindfulness.node';
dlmwrite(filename_node, mindfulness_node, 'delimiter', '\t');
%% Hemodynamic change extraction for each participant
%% For the correlation analysis with delay discounting change
%% Endpoint-focus meditation
[num_timepoints, num_total_columns] = size(ses_endpoint_hbo);
num_subjects = num_total_columns / num_channels;
mean_activation_endpoint = zeros(num_subjects, num_channels); %subj*channel
for subj = 1:num_subjects
    start_col = (subj - 1) * num_channels + 1;
    end_col = subj * num_channels;         
    subj_data = ses_endpoint_hbo(:, start_col:end_col); 
    mean_activation_endpoint(subj, :) = mean(subj_data, 1); % mean extraction
end
%% Mindfulness meditation
[num_timepoints, num_total_columns] = size(ses_mindfulness_hbo);
num_subjects = num_total_columns / num_channels;
mean_activation_mindfulness = zeros(num_subjects, num_channels); %subj*channel
for subj = 1:num_subjects
    start_col = (subj - 1) * num_channels + 1;
    end_col = subj * num_channels;         
    subj_data = ses_mindfulness_hbo(:, start_col:end_col); 
    mean_activation_mindfulness(subj, :) = mean(subj_data, 1); % mean extraction
end