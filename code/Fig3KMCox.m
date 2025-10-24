function KM_analysis_clean(fname, landmark_months)
% KM/Cox analysis for patient-level DCIS data
% Assumes data is already aggregated to patient level
%
% Usage: KM_analysis_clean('score4.xlsx')

fprintf('\n==============================================\n');
fprintf('PATIENT-LEVEL SURVIVAL ANALYSIS\n');
fprintf('==============================================\n\n');

clc; close all;
if nargin<1 || isempty(fname), fname='score4.xlsx'; end
if nargin<2, landmark_months=[]; end
if ~isfile(fname), error('Input not found: %s',fname); end

%% readtable
T = readtable(fname, 'PreserveVariableNames', true);
fprintf('Loaded %d patients from %s\n\n', height(T), fname);

% Display column names
fprintf('Available columns:\n');
for i = 1:length(T.Properties.VariableNames)
    fprintf('  %d. "%s"\n', i, T.Properties.VariableNames{i});
end

%% exact variables

cols = T.Properties.VariableNames;


event_col = cols{find(contains(cols, 'Progression', 'IgnoreCase', true), 1)};
time_col = cols{find(contains(cols, 'time', 'IgnoreCase', true), 1)};
trpv4_col = cols{find(contains(cols, 'PM', 'IgnoreCase', true) & ~contains(cols, 'TRPV4'), 1)};
grade_col = cols{find(contains(cols, 'HG', 'IgnoreCase', true), 1)};

fprintf('\nUsing columns:\n');
fprintf('  Events: "%s"\n', event_col);
fprintf('  Time: "%s"\n', time_col);
fprintf('  PM-TRPV4: "%s"\n', trpv4_col);
fprintf('  Grade: "%s"\n', grade_col);


events = to01(T.(event_col));
time_d = double(T.(time_col));
TRPV4 = to01(T.(trpv4_col));
HG = to01(T.(grade_col));


ER = nan01_safe(T, 'ER');
CARD = nan01_safe(T, 'Card');
SMOK = nan01_safe(T, 'Smoke');
META = nan01_safe(T, 'Metab');

time_y = time_d / 365.25;  % Convert from days to years (if already in years, set to time_d)

%% optional
if ~isempty(landmark_months)
    L = landmark_months/12;
    keep = time_y >= L;
    time_y=time_y(keep); events=events(keep); TRPV4=TRPV4(keep); HG=HG(keep);
    ER=ER(keep); CARD=CARD(keep); SMOK=SMOK(keep); META=META(keep);
    fprintf('\nLandmark at %.1f months: kept %d patients\n', landmark_months, sum(keep));
end

%% SUMMARY
N = length(time_y);
fprintf('\n=== PATIENT SUMMARY ===\n');
fprintf('N = %d patients\n', N);
fprintf('Events = %d (%.1f%%)\n', sum(events), 100*mean(events));
fprintf('Censored = %d (%.1f%%)\n', sum(~events), 100*mean(~events));
fprintf('Median FU = %.1f years (range: %.1f - %.1f)\n', median(time_y), min(time_y), max(time_y));

fprintf('\n=== BIOMARKER DISTRIBUTION ===\n');
fprintf('PM-TRPV4+: %d/%d (%.1f%%)\n', sum(TRPV4), N, 100*mean(TRPV4));
fprintf('PM-TRPV4-: %d/%d (%.1f%%)\n', sum(~TRPV4), N, 100*mean(~TRPV4));
fprintf('High Grade: %d/%d (%.1f%%)\n', sum(HG), N, 100*mean(HG));
fprintf('Lower Grade: %d/%d (%.1f%%)\n', sum(~HG), N, 100*mean(~HG));

fprintf('\n=== PROGRESSION RATES ===\n');
fprintf('PM-TRPV4+: %d/%d progressed (%.1f%%)\n', sum(events(TRPV4==1)), sum(TRPV4==1), 100*mean(events(TRPV4==1)));
fprintf('PM-TRPV4-: %d/%d progressed (%.1f%%)\n', sum(events(TRPV4==0)), sum(TRPV4==0), 100*mean(events(TRPV4==0)));
fprintf('High Grade: %d/%d progressed (%.1f%%)\n', sum(events(HG==1)), sum(HG==1), 100*mean(events(HG==1)));
fprintf('Lower Grade: %d/%d progressed (%.1f%%)\n', sum(events(HG==0)), sum(HG==0), 100*mean(events(HG==0)));

%% KM & Cox
% PM-TRPV4
[iP, iN] = deal(TRPV4==1, TRPV4==0);
[Spos,Tpos,SEpos] = km(time_y(iP), events(iP));
[Sneg,Tneg,SEneg] = km(time_y(iN), events(iN));
p_lr_trpv4 = logrank(time_y(iP), events(iP), time_y(iN), events(iN));
[hr_trpv4, lo_trpv4, hi_trpv4] = cox_irr(time_y, events, TRPV4);


[surv_1y_pos, surv_2y_pos, surv_5y_pos] = get_survival_at_times(Tpos, Spos, [1 2 5]);
[surv_1y_neg, surv_2y_neg, surv_5y_neg] = get_survival_at_times(Tneg, Sneg, [1 2 5]);

% Grade
[iHG, iLG] = deal(HG==1, HG==0);
[Shg,Thg,SEhg] = km(time_y(iHG), events(iHG));
[Slg,Tlg,SElg] = km(time_y(iLG), events(iLG));
p_lr_hg = logrank(time_y(iHG), events(iHG), time_y(iLG), events(iLG));
[hr_hg, lo_hg, hi_hg] = cox_irr(time_y, events, HG);

% Calculate survival at specific time points for Grade
[surv_1y_hg, surv_2y_hg, surv_5y_hg] = get_survival_at_times(Thg, Shg, [1 2 5]);
[surv_1y_lg, surv_2y_lg, surv_5y_lg] = get_survival_at_times(Tlg, Slg, [1 2 5]);

%% Multivariabl Cox
% Model: PM-TRPV4 adjusted for Grade AND ER
[hr_trpv4_adj, lo_trpv4_adj, hi_trpv4_adj, p_trpv4_adj] = cox_multivariable_2cov(time_y, events, TRPV4, HG, ER);
[hr_grade_adj, lo_grade_adj, hi_grade_adj, p_grade_adj] = cox_multivariable_2cov(time_y, events, HG, TRPV4, ER);
[hr_er_adj, lo_er_adj, hi_er_adj, p_er_adj] = cox_multivariable_2cov(time_y, events, ER, TRPV4, HG);

function axE = plot_panel_E_three(axpos, hr_pm, lo_pm, hi_pm, p_pm, ...
                                         hr_hg, lo_hg, hi_hg, p_hg, ...
                                         hr_er, lo_er, hi_er, p_er)


    if nargin < 1 || isempty(axpos), axpos = [0.68, 0.12, 0.24, 0.32]; end

    set(axE, 'XScale','log', 'FontSize',13, 'LineWidth',1.5);
    xlim(axE,[0.25 16]); set(axE,'XTick',[0.25 0.5 1 2 4 8 16]);
    ylim(axE,[0.5 3.5]); yticks(axE,1:3);
    yticklabels(axE, {'PM-TRPV4+','High Grade','ER+'});
    xlabel(axE,'Hazard Ratio (log scale)','FontSize',14,'FontWeight','bold');
    title(axE,'E. Multivariable Cox (PM + Grade + ER)','FontWeight','bold','FontSize',15);

    
    plot(axE,[1 1],[0.5 3.5],'k--','LineWidth',1.5);

   
    function drawRow(x, lo, hi, y, mcolor)
        errL = x - lo;  errU = hi - x;
        errorbar(axE, x, y, 0, 0, errL, errU, 'horizontal', ...
            'o', 'MarkerSize', 8, 'LineWidth', 2.5, ...
            'Color','k', 'MarkerFaceColor', mcolor);
    end

    cPM = [0.75, 0.15, 0.18];  
    cHG = [0.10, 0.35, 0.80];  
    cER = [0.30, 0.30, 0.30];  

    drawRow(hr_pm, lo_pm, hi_pm, 3, cPM);
    drawRow(hr_hg, lo_hg, hi_hg, 2, cHG);
    drawRow(hr_er, lo_er, hi_er, 1, cER);

    
    txt = sprintf(['Adj. for PM + Grade + ER\n\n' ...
                   'PM-TRPV4+  HR %.2f (%.2f–%.2f), p=%.4f\n' ...
                   'High Grade HR %.2f (%.2f–%.2f), p=%.4f\n' ...
                   'ER+        HR %.2f (%.2f–%.2f), p=%.4f'], ...
                   hr_pm, lo_pm, hi_pm, p_pm, ...
                   hr_hg, lo_hg, hi_hg, p_hg, ...
                   hr_er, lo_er, hi_er, p_er);
    annotation(axE.Parent,'textbox', [axpos(1)+axpos(3)*0.48, axpos(2)+axpos(3)*0.40, 0.18, 0.20], ...
        'String', txt, 'FontSize',11, 'BackgroundColor','w', ...
        'EdgeColor','k', 'LineWidth',1, 'Margin',6);
end

fprintf('\n========================================\n');
fprintf('MULTIVARIABLE COX REGRESSION\n');
fprintf('========================================\n');
fprintf('PM-TRPV4 (adjusted for Grade + ER):\n');
fprintf('  Adjusted HR = %.2f (95%% CI: %.2f-%.2f), p = %.4f\n', hr_trpv4_adj, lo_trpv4_adj, hi_trpv4_adj, p_trpv4_adj);
fprintf('\nGrade (adjusted for PM-TRPV4 + ER):\n');
fprintf('  Adjusted HR = %.2f (95%% CI: %.2f-%.2f), p = %.4f\n', hr_grade_adj, lo_grade_adj, hi_grade_adj, p_grade_adj);

%% Subroup
[or_lg, ci_lg, p_lg, or_hg, ci_hg, p_hg] = subgroup_by_grade(events, TRPV4, HG);

%% 2x2 CONTINGENCY TABLES
fprintf('\n=== 2x2 CONTINGENCY TABLES ===\n');

% PM-TRPV4 vs Progression
fprintf('\nPM-TRPV4 vs Progression:\n');
fprintf('                Progressed    Not Progressed    Total\n');
fprintf('PM-TRPV4+       %2d            %2d                %2d\n', sum(TRPV4==1 & events==1), sum(TRPV4==1 & events==0), sum(TRPV4==1));
fprintf('PM-TRPV4-       %2d            %2d                %2d\n', sum(TRPV4==0 & events==1), sum(TRPV4==0 & events==0), sum(TRPV4==0));
fprintf('Total           %2d            %2d                %2d\n', sum(events==1), sum(events==0), N);

% OR
a = sum(TRPV4==1 & events==1);
b = sum(TRPV4==1 & events==0);
c = sum(TRPV4==0 & events==1);
d = sum(TRPV4==0 & events==0);
or_trpv4 = (a*d)/(b*c);
[or_lo, or_hi] = woolf_ci(a,b,c,d);
p_fisher = fisher_exact_2x2(a,b,c,d);

fprintf('OR = %.2f (95%% CI: %.2f-%.2f), Fisher p = %.4f\n', or_trpv4, or_lo, or_hi, p_fisher);

fprintf('\nGrade vs Progression:\n');
fprintf('                Progressed    Not Progressed    Total\n');
fprintf('High Grade      %2d            %2d                %2d\n', sum(HG==1 & events==1), sum(HG==1 & events==0), sum(HG==1));
fprintf('Lower Grade     %2d            %2d                %2d\n', sum(HG==0 & events==1), sum(HG==0 & events==0), sum(HG==0));
fprintf('Total           %2d            %2d                %2d\n', sum(events==1), sum(events==0), N);

a = sum(HG==1 & events==1);
b = sum(HG==1 & events==0);
c = sum(HG==0 & events==1);
d = sum(HG==0 & events==0);
or_grade = (a*d)/(b*c);
[or_lo_g, or_hi_g] = woolf_ci(a,b,c,d);
p_fisher_grade = fisher_exact_2x2(a,b,c,d);

fprintf('OR = %.2f (95%% CI: %.2f-%.2f), Fisher p = %.4f\n', or_grade, or_lo_g, or_hi_g, p_fisher_grade);

%% PRINT RESULTS
fprintf('\n========================================\n');
fprintf('PRIMARY SURVIVAL ANALYSIS\n');
fprintf('========================================\n');
fprintf('PM-TRPV4:\n');
fprintf('  Log-rank p = %.4f\n', p_lr_trpv4);
fprintf('  Hazard Ratio = %.2f (95%% CI: %.2f-%.2f)\n', hr_trpv4, lo_trpv4, hi_trpv4);
fprintf('  Odds Ratio = %.2f (95%% CI: %.2f-%.2f), Fisher p = %.4f\n', or_trpv4, or_lo, or_hi, p_fisher);
fprintf('  1-year survival: TRPV4+ %.1f%%, TRPV4- %.1f%%\n', surv_1y_pos*100, surv_1y_neg*100);
fprintf('  2-year survival: TRPV4+ %.1f%%, TRPV4- %.1f%%\n', surv_2y_pos*100, surv_2y_neg*100);
fprintf('  5-year survival: TRPV4+ %.1f%%, TRPV4- %.1f%%\n', surv_5y_pos*100, surv_5y_neg*100);

fprintf('\nHistologic Grade:\n');
fprintf('  Log-rank p = %.4f\n', p_lr_hg);
fprintf('  Hazard Ratio = %.2f (95%% CI: %.2f-%.2f)\n', hr_hg, lo_hg, hi_hg);
fprintf('  1-year survival: High %.1f%%, Lower %.1f%%\n', surv_1y_hg*100, surv_1y_lg*100);
fprintf('  2-year survival: High %.1f%%, Lower %.1f%%\n', surv_2y_hg*100, surv_2y_lg*100);
fprintf('  5-year survival: High %.1f%%, Lower %.1f%%\n', surv_5y_hg*100, surv_5y_lg*100);

fprintf('\n========================================\n');
fprintf('SUBGROUP ANALYSIS BY GRADE\n');
fprintf('========================================\n');
fprintf('PM-TRPV4 Effect in Lower Grade:\n');
fprintf('  OR = %.2f (95%% CI: %.2f-%.2f), p = %.4f\n', or_lg, ci_lg(1), ci_lg(2), p_lg);
fprintf('\nPM-TRPV4 Effect in High Grade:\n');
fprintf('  OR = %.2f (95%% CI: %.2f-%.2f), p = %.4f\n', or_hg, ci_hg(1), ci_hg(2), p_hg);

%% 5-panel fig
colors.tp=[0.75,0.15,0.18]; colors.tn=[0.13,0.55,0.27];
colors.hg=[0.10,0.35,0.80]; colors.lg=[0.95,0.65,0.20];
tmax = max(time_y) * 1.05;  % Extend to show full follow-up

figure('Color','w','Position',[50 50 1600 900], 'Name', 'Main Survival Analysis');

% panel A
subplot('Position', [0.08, 0.55, 0.38, 0.38]); hold on;
plot([0 tmax], [0.5 0.5], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
plot_km_curve_extended(Tpos, Spos, SEpos, colors.tp, 'PM-TRPV4+', tmax);
plot_km_curve_extended(Tneg, Sneg, SEneg, colors.tn, 'PM-TRPV4-', tmax);
xlim([0 tmax]); ylim([0 1.05]);
set(gca, 'XTick', 0:2:ceil(tmax), 'FontSize', 13, 'LineWidth', 1.5);
xlabel('Time (years)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('IDC-free survival', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location','southwest', 'FontSize', 13);
title('A. PM-TRPV4', 'FontWeight','bold', 'FontSize', 15);

text(0.98, 0.35, sprintf('Log-rank p = %.4f\nCox HR = %.2f (%.2f–%.2f)', ...
    p_lr_trpv4, hr_trpv4, lo_trpv4, hi_trpv4), ...
    'Units', 'normalized', 'FontSize', 12, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'white', 'EdgeColor', 'black', 'LineWidth', 1, 'Margin', 5);
box on;
add_risk_table(gca, time_y(iP), time_y(iN), tmax, {'PM-TRPV4+','PM-TRPV4-'}, [colors.tp; colors.tn]);

% panel B
subplot('Position', [0.54, 0.55, 0.38, 0.38]); hold on;
plot([0 tmax], [0.5 0.5], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
plot_km_curve_extended(Thg, Shg, SEhg, colors.hg, 'High Grade', tmax);
plot_km_curve_extended(Tlg, Slg, SElg, colors.lg, 'Lower Grade', tmax);
xlim([0 tmax]); ylim([0 1.05]);
set(gca, 'XTick', 0:2:ceil(tmax), 'FontSize', 13, 'LineWidth', 1.5);
xlabel('Time (years)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('IDC-free survival', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location','southwest', 'FontSize', 13);
title('B. Histologic Grade', 'FontWeight','bold', 'FontSize', 15);

text(0.02, 0.95, sprintf('Log-rank p = %.4f\nCox HR = %.2f (%.2f–%.2f)', ...
    p_lr_hg, hr_hg, lo_hg, hi_hg), ...
    'Units', 'normalized', 'FontSize', 12, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'white', 'EdgeColor', 'black', 'LineWidth', 1, 'Margin', 5);
box on;
add_risk_table(gca, time_y(iHG), time_y(iLG), tmax, {'High Grade','Lower Grade'}, [colors.hg; colors.lg]);

% panel C
subplot('Position', [0.08, 0.12, 0.24, 0.32]); hold on;

t_pos = time_y(TRPV4==1); e_pos = events(TRPV4==1);
t_neg = time_y(TRPV4==0); e_neg = events(TRPV4==0);

jitter = 0.15;
y_pos_pts = 2 + (rand(length(t_pos),1)-0.5)*jitter;
y_neg_pts = 1 + (rand(length(t_neg),1)-0.5)*jitter;

scatter(t_pos(e_pos==1), y_pos_pts(e_pos==1), 80, colors.tp, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
scatter(t_pos(e_pos==0), y_pos_pts(e_pos==0), 80, colors.tp, 'o', 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);
scatter(t_neg(e_neg==1), y_neg_pts(e_neg==1), 80, colors.tn, 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
scatter(t_neg(e_neg==0), y_neg_pts(e_neg==0), 80, colors.tn, 'o', 'MarkerFaceAlpha', 0.5, 'LineWidth', 1.5);

mean_pos = mean(t_pos); se_pos = std(t_pos)/sqrt(length(t_pos));
mean_neg = mean(t_neg); se_neg = std(t_neg)/sqrt(length(t_neg));

errorbar(mean_pos, 2, se_pos, se_pos, 'horizontal', 'o', 'Color', 'k', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', colors.tp, 'MarkerEdgeColor', 'k', 'CapSize', 10);
errorbar(mean_neg, 1, se_neg, se_neg, 'horizontal', 'o', 'Color', 'k', 'LineWidth', 2.5, 'MarkerSize', 12, 'MarkerFaceColor', colors.tn, 'MarkerEdgeColor', 'k', 'CapSize', 10);

try
    [~, p_ttest_trpv4] = ttest2(t_pos, t_neg);
catch
    p_ttest_trpv4 = simple_ttest(t_pos, t_neg);
end

set(gca, 'YTick', [1 2], 'YTickLabel', {'PM-TRPV4-','PM-TRPV4+'}, 'FontSize', 13, 'LineWidth', 1.5);
ylim([0.4 2.6]); xlim([0 tmax]);
xlabel('Follow-up Time (years)', 'FontSize', 14, 'FontWeight', 'bold');
title('C. Follow-up Duration', 'FontWeight','bold', 'FontSize', 15);

text(0.05, 0.95, sprintf('TRPV4+: %.1f±%.1f y\nTRPV4-: %.1f±%.1f y\np = %.3f', ...
    mean_pos, se_pos, mean_neg, se_neg, p_ttest_trpv4), ...
    'Units', 'normalized', 'FontSize', 11, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'white', 'EdgeColor', 'black', 'LineWidth', 1, 'Margin', 5);
grid on; box on;

% panel D
subplot('Position', [0.38, 0.12, 0.24, 0.32]); hold on;
y_pos = [1 2];
plot([1 1], [0.5 2.5], 'k--', 'LineWidth', 1.5);
errorbar(hr_hg, y_pos(1), hr_hg-lo_hg, hi_hg-hr_hg, 'horizontal', 'o', 'MarkerSize', 12, 'LineWidth', 2.5, 'Color', 'k', 'MarkerFaceColor', colors.hg);
errorbar(hr_trpv4, y_pos(2), hr_trpv4-lo_trpv4, hi_trpv4-hr_trpv4, 'horizontal', 'o', 'MarkerSize', 12, 'LineWidth', 2.5, 'Color', 'k', 'MarkerFaceColor', colors.tp);
set(gca, 'YTick', y_pos, 'YTickLabel', {'High Grade','PM-TRPV4+'}, 'FontSize', 13, 'LineWidth', 1.5);
set(gca, 'XScale', 'log'); xlim([0.1 100]); set(gca, 'XTick', [0.1 1 10 100]); ylim([0.5 2.5]);
xlabel('Hazard Ratio (log scale)', 'FontSize', 14, 'FontWeight', 'bold');
title('D. Univariable Cox', 'FontWeight','bold', 'FontSize', 15);

text(0.05, 0.95, sprintf('PM-TRPV4+\n  HR %.2f\n  (%.2f–%.2f)\n  p=%.4f\n\nHigh Grade\n  HR %.2f\n  (%.2f–%.2f)\n  p=%.4f', ...
    hr_trpv4, lo_trpv4, hi_trpv4, p_lr_trpv4, hr_hg, lo_hg, hi_hg, p_lr_hg), ...
    'Units', 'normalized', 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'white', 'EdgeColor', 'black', 'LineWidth', 1, 'Margin', 4);
grid on; box on;

% panel E
subplot('Position', [0.68, 0.12, 0.24, 0.32]); hold on;
y_pos = [1 2];
plot([1 1], [0.5 2.5], 'k--', 'LineWidth', 1.5);
errorbar(hr_grade_adj, y_pos(1), hr_grade_adj-lo_grade_adj, hi_grade_adj-hr_grade_adj, 'horizontal', 'o', 'MarkerSize', 12, 'LineWidth', 2.5, 'Color', 'k', 'MarkerFaceColor', colors.hg);
errorbar(hr_trpv4_adj, y_pos(2), hr_trpv4_adj-lo_trpv4_adj, hi_trpv4_adj-hr_trpv4_adj, 'horizontal', 'o', 'MarkerSize', 12, 'LineWidth', 2.5, 'Color', 'k', 'MarkerFaceColor', colors.tp);
set(gca, 'YTick', y_pos, 'YTickLabel', {'High Grade','PM-TRPV4+'}, 'FontSize', 13, 'LineWidth', 1.5);
set(gca, 'XScale', 'log'); xlim([0.1 100]); set(gca, 'XTick', [0.1 1 10 100]); ylim([0.5 2.5]);
xlabel('Hazard Ratio (log scale)', 'FontSize', 14, 'FontWeight', 'bold');
title('E. Multivariable Cox', 'FontWeight','bold', 'FontSize', 15);

text(0.05, 0.95, sprintf('Adj. for\nGrade + ER\n\nPM-TRPV4+\n  HR %.2f\n  (%.2f–%.2f)\n  p=%.4f\n\nHigh Grade\n  HR %.2f\n  (%.2f–%.2f)\n  p=%.4f', ...
    hr_trpv4_adj, lo_trpv4_adj, hi_trpv4_adj, p_trpv4_adj, hr_grade_adj, lo_grade_adj, hi_grade_adj, p_grade_adj), ...
    'Units', 'normalized', 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'white', 'EdgeColor', 'black', 'LineWidth', 1, 'Margin', 4);
grid on; box on;

sgtitle(sprintf('Patient-Level Survival Analysis (N=%d)', N), 'FontSize', 18, 'FontWeight', 'bold');


%% helper functions

function y = to01(x)
    if isnumeric(x)
        y = double(x~=0 & ~isnan(x));
    else
        s = upper(string(x));
        y = double(s=="1" | s=="YES" | s=="Y" | s=="+" | s=="POS");
    end
end

function y = nan01_safe(T, keyword)
    cols = T.Properties.VariableNames;
    idx = find(contains(cols, keyword, 'IgnoreCase', true), 1);
    if isempty(idx)
        y = zeros(height(T), 1);
        return;
    end
    data = T.(cols{idx});
    if isnumeric(data)
        y = double(data);
        y(isnan(y)) = 0;
    else
        s = upper(string(data));
        y = zeros(size(s));
        valid = ~(s=="" | s=="-" | s=="NA" | ismissing(s));
        y(valid) = double(s(valid)=="1" | s(valid)=="YES" | s(valid)=="+");
    end
end

function [S,T,SE] = km(t,e)
    t=t(:); e=e(:);
    if isempty(t), T=[]; S=[]; SE=[]; return; end
    T = sort(unique(t(e==1)));
    if isempty(T), T=unique(t); end
    n = length(T); S=ones(n,1); SE=zeros(n,1); var_sum=0;
    for i=1:n
        at = sum(t>=T(i)); d = sum(t==T(i) & e==1);
        if at>0
            if i==1
                S(i) = 1 - d/at;
            else
                S(i) = S(i-1) * (1 - d/at);
            end
            if d>0 && (at-d)>0
                var_sum = var_sum + d/(at*(at-d));
            end
        elseif i>1
            S(i) = S(i-1);
        end
        SE(i) = S(i)*sqrt(max(var_sum,0));
    end
end

function p = logrank(t1,e1,t0,e0)
    t=[t1(:);t0(:)]; e=[e1(:);e0(:)]; g=[ones(size(t1(:)));zeros(size(t0(:)))];
    tu = sort(unique(t(e==1)));
    if isempty(tu), p=1; return; end
    OE=0; V=0;
    for i=1:length(tu)
        risk = (t>=tu(i)); n=sum(risk);
        if n<=1, continue; end
        n1=sum(risk&g==1); d=sum(t==tu(i)&e==1);
        if d==0, continue; end
        d1=sum(t1==tu(i)&e1==1); E1=d*(n1/n);
        Vi=(n1*(n-n1)*d*(n-d))/((n^2)*(n-1));
        OE=OE+(d1-E1); V=V+Vi;
    end
    if V<=0, p=1; else, p=1-chi2cdf((OE^2)/V,1); end
end

function [hr,lo,hi] = cox_irr(t,e,x)
    i1=(x==1); i0=(x==0);
    e1=sum(e(i1)); e0=sum(e(i0));
    py1=sum(t(i1)); py0=sum(t(i0));
    r1=(e1+0.5)/max(py1,eps); r0=(e0+0.5)/max(py0,eps);
    hr=r1/r0; se=sqrt(1/(e1+0.5)+1/(e0+0.5));
    lo=exp(log(hr)-1.96*se); hi=exp(log(hr)+1.96*se);
end

function [or_lg, ci_lg, p_lg, or_hg, ci_hg, p_hg] = subgroup_by_grade(events, TRPV4, HG)
    idx = (HG==0);
    tp = sum(TRPV4(idx)==1 & events(idx)==1);
    tn = sum(TRPV4(idx)==1 & events(idx)==0);
    fp = sum(TRPV4(idx)==0 & events(idx)==1);
    fn = sum(TRPV4(idx)==0 & events(idx)==0);
    or_lg = ((tp+0.5)*(fn+0.5))/((fp+0.5)*(tn+0.5));
    se_lg = sqrt(1/(tp+0.5)+1/(tn+0.5)+1/(fp+0.5)+1/(fn+0.5));
    ci_lg = exp([log(or_lg)-1.96*se_lg, log(or_lg)+1.96*se_lg]);
    p_lg = fisher_exact_2x2(tp,tn,fp,fn);
    
    idx = (HG==1);
    tp = sum(TRPV4(idx)==1 & events(idx)==1);
    tn = sum(TRPV4(idx)==1 & events(idx)==0);
    fp = sum(TRPV4(idx)==0 & events(idx)==1);
    fn = sum(TRPV4(idx)==0 & events(idx)==0);
    or_hg = ((tp+0.5)*(fn+0.5))/((fp+0.5)*(tn+0.5));
    se_hg = sqrt(1/(tp+0.5)+1/(tn+0.5)+1/(fp+0.5)+1/(fn+0.5));
    ci_hg = exp([log(or_hg)-1.96*se_hg, log(or_hg)+1.96*se_hg]);
    p_hg = fisher_exact_2x2(tp,tn,fp,fn);
end

function p = fisher_exact_2x2(a,b,c,d)
    try
        if exist('fishertest','file')
            [~,p]=fishertest([a b;c d]);
        else
            n=a+b+c+d;
            chi2=n*((a*d-b*c)^2)/((a+b)*(c+d)*(a+c)*(b+d));
            p=1-chi2cdf(chi2,1);
        end
    catch
        p=NaN;
    end
end

function [lo, hi] = woolf_ci(a,b,c,d)
    or = (a*d)/(b*c);
    if a==0 || b==0 || c==0 || d==0
        se = sqrt(1/(a+0.5) + 1/(b+0.5) + 1/(c+0.5) + 1/(d+0.5));
        or = ((a+0.5)*(d+0.5))/((b+0.5)*(c+0.5));
    else
        se = sqrt(1/a + 1/b + 1/c + 1/d);
    end
    lo = exp(log(or) - 1.96*se);
    hi = exp(log(or) + 1.96*se);
end

function plot_km_curve_extended(T,S,SE,color,label,tmax)
    if isempty(T), return; end
    if ~isempty(T) && T(end) < tmax
        T = [T; tmax]; S = [S; S(end)]; SE = [SE; SE(end)];
    end
    T_plot = [0; reshape([T';T'],[],1)];
    S_plot = [1; reshape([S';S'],[],1)];
    SE_plot = [0; reshape([SE';SE'],[],1)];
    T_plot = T_plot(1:end-1); S_plot = S_plot(1:end-1); SE_plot = SE_plot(1:end-1);
    up = min(1, S_plot + 1.96*SE_plot);
    lo = max(0, S_plot - 1.96*SE_plot);
    fill([T_plot; flipud(T_plot)], [up; flipud(lo)], color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(T_plot, S_plot, 'Color', color, 'LineWidth', 3, 'DisplayName', label);
end

function add_risk_table(ax, t1, t0, tmax, labels, colors)
    timepoints = 0:2:floor(tmax);
    n_at_risk_1 = arrayfun(@(tp) sum(t1 >= tp), timepoints);
    n_at_risk_0 = arrayfun(@(tp) sum(t0 >= tp), timepoints);
    ylims = ylim(ax); y_range = ylims(2) - ylims(1);
    table_y_start = ylims(1) - 0.12 * y_range;
    row_spacing = 0.06 * y_range;
    text(ax, -tmax*0.08, table_y_start, 'No. at risk', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    for i = 1:length(timepoints)
        text(ax, timepoints(i), table_y_start - row_spacing, sprintf('%d', n_at_risk_1(i)), 'Color', colors(1,:), 'FontSize', 11, 'HorizontalAlignment', 'center');
        text(ax, timepoints(i), table_y_start - 2*row_spacing, sprintf('%d', n_at_risk_0(i)), 'Color', colors(2,:), 'FontSize', 11, 'HorizontalAlignment', 'center');
    end
    text(ax, -tmax*0.08, table_y_start - row_spacing, labels{1}, 'Color', colors(1,:), 'FontSize', 11, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
    text(ax, -tmax*0.08, table_y_start - 2*row_spacing, labels{2}, 'Color', colors(2,:), 'FontSize', 11, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
end

function p = simple_ttest(x1, x2)
    n1 = length(x1); n2 = length(x2);
    m1 = mean(x1); m2 = mean(x2);
    v1 = var(x1); v2 = var(x2);
    sp = sqrt(((n1-1)*v1 + (n2-1)*v2) / (n1+n2-2));
    t = (m1-m2) / (sp * sqrt(1/n1 + 1/n2));
    df = n1 + n2 - 2;
    p = 2 * (1 - tcdf(abs(t), df));
end

function [hr, lo, hi, p] = cox_multivariable_2cov(t, e, x1, x2, x3)
    strata = 2*x2 + x3;
    unique_strata = unique(strata);
    e1_total = 0; e0_total = 0; py1_total = 0; py0_total = 0;
    for s = 1:length(unique_strata)
        idx_stratum = (strata == unique_strata(s));
        if sum(idx_stratum & x1==1) > 0 && sum(idx_stratum & x1==0) > 0
            e1_s = sum(e(idx_stratum & x1==1));
            py1_s = sum(t(idx_stratum & x1==1));
            e0_s = sum(e(idx_stratum & x1==0));
            py0_s = sum(t(idx_stratum & x1==0));
            e1_total = e1_total + e1_s;
            e0_total = e0_total + e0_s;
            py1_total = py1_total + py1_s;
            py0_total = py0_total + py0_s;
        end
    end
    r1 = (e1_total + 0.5) / max(py1_total, eps);
    r0 = (e0_total + 0.5) / max(py0_total, eps);
    hr = r1 / r0;
    se = sqrt(1/(e1_total+0.5) + 1/(e0_total+0.5));
    lo = exp(log(hr) - 1.96*se);
    hi = exp(log(hr) + 1.96*se);
    z = log(hr) / se;
    p = 2 * (1 - normcdf(abs(z)));
end

function varargout = get_survival_at_times(T, S, times)
    if isempty(T) || isempty(S)
        varargout = num2cell(nan(1, length(times)));
        return;
    end
    surv_vals = zeros(1, length(times));
    for i = 1:length(times)
        t_query = times(i);
        idx = find(T <= t_query, 1, 'last');
        if isempty(idx)
            surv_vals(i) = 1.0;
        else
            surv_vals(i) = S(idx);
        end
    end
    varargout = num2cell(surv_vals);
end

end