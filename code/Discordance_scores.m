function trpv4_scores_from_excel(fname)
% TRPV4 scoring analysis 


if nargin<1 || isempty(fname), fname = 'score.xlsx'; end
if ~isfile(fname), error('Cannot find input file: %s', fname); end
fprintf('>>> Loading data from: %s\n', fname);
fprintf('    NOTE: "0.5" in Excel is automatically converted to "0-1" for analysis\n');


CATS   = ["0","0-1","1","2","3"];
GRADE3 = ["LG","IMG","HG"];

SCORE_COLORS = [
    0.90, 0.90, 0.90;  % "0"   - light gray
    0.70, 0.85, 0.95;  % "0-1" - light blue
    0.40, 0.76, 0.65;  % "1"   - teal
    0.99, 0.70, 0.38;  % "2"   - orange
    0.84, 0.19, 0.15;  % "3"   - red
];


opts = detectImportOptions(fname, 'FileType','spreadsheet');
opts.VariableNamingRule = 'preserve';
T0 = readtable(fname, opts);

V  = string(T0.Properties.VariableNames);
% Try both uppercase and lowercase variants
idxA = find(strcmpi(V,'Pathologist A') | strcmpi(V,'pathologist A'),1);
idxB = find(strcmpi(V,'Pathologist B') | strcmpi(V,'pathologist B'),1);
idxC = find(strcmpi(V,'Pathologist C') | strcmpi(V,'pathologist C'),1);
idxP = find(strcmpi(V,'Pathology (Annotated Slides)') | strcmpi(V,'pathology (annotated slides)'),1);

if any(cellfun(@isempty,{idxA,idxB,idxC,idxP}))
    actual_headers = strjoin(T0.Properties.VariableNames, ' | ');
    error(['Missing required columns. Expected (case-insensitive):\n' ...
           '  ''Pathologist A'' | ''Pathologist B'' | ''Pathologist C'' | ''Pathology (Annotated Slides)''\n' ...
           'Found:\n  ''%s'''], actual_headers);
end


L_raw = T0{:,idxA};
S_raw = T0{:,idxB};
A_raw = T0{:,idxC};
P_raw = T0{:,idxP};


if isnumeric(L_raw)
    L = arrayfun(@(x) num2str_custom(x), L_raw, 'UniformOutput', false);
    L = string(L);
else
    L = string(L_raw);
    L(ismissing(L)) = "";  
end
if isnumeric(S_raw)
    S = arrayfun(@(x) num2str_custom(x), S_raw, 'UniformOutput', false);
    S = string(S);
else
    S = string(S_raw);
    S(ismissing(S)) = "";  
end
if isnumeric(A_raw)
    A = arrayfun(@(x) num2str_custom(x), A_raw, 'UniformOutput', false);
    A = string(A);
else
    A = string(A_raw);
    A(ismissing(A)) = "";  
end
P = string(P_raw);
P(ismissing(P)) = "";  

fprintf('Raw data loaded: %d rows\n', numel(L));


fprintf('\nDEBUG - Raw values (first 10 rows):\n');
for i = 1:min(10, numel(L))
    % Handle missing values for display
    L_display = string(L_raw(i)); if ismissing(L_display), L_display = "<missing>"; end
    S_display = string(S_raw(i)); if ismissing(S_display), S_display = "<missing>"; end
    A_display = string(A_raw(i)); if ismissing(A_display), A_display = "<missing>"; end
    fprintf('  Row %d: A_raw="%s" B_raw="%s" C_raw="%s"\n', i, L_display, S_display, A_display);
end

L = normTok(L); S = normTok(S); A = normTok(A);

fprintf('\nDEBUG - Unique values after normalization (0.5→0-1):\n');
fprintf('  Pathologist A: %s\n', strjoin(unique(L(~ismissing(L) & L~="")), ', '));
fprintf('  Pathologist B: %s\n', strjoin(unique(S(~ismissing(S) & S~="")), ', '));
fprintf('  Pathologist C: %s\n', strjoin(unique(A(~ismissing(A) & A~="")), ', '));

[G, keep] = merge_grade_from_pathology(P);
dropped = sum(~keep);
if dropped>0
    fprintf('Dropped %d row(s) that were non-DCIS:\n', dropped);
    dropped_idx = find(~keep);
    for i = 1:min(5, numel(dropped_idx))
        fprintf('  Row %d: "%s"\n', dropped_idx(i), P(dropped_idx(i)));
    end
end
L = L(keep); S = S(keep); A = A(keep); G = G(keep);


isValid = ismember(L,CATS) & ismember(S,CATS) & ismember(A,CATS);
if ~all(isValid)
    invalid_count = sum(~isValid);
    fprintf('Dropping %d row(s) with missing/invalid ratings.\n', invalid_count);
    
    
    invalid_idx = find(~isValid);
    fprintf('\nDEBUG - First 10 invalid rows:\n');
    for i = 1:min(10, numel(invalid_idx))
        idx = invalid_idx(i);
        fprintf('  Row %d: A="%s" B="%s" C="%s" (in CATS: %d %d %d)\n', ...
            idx, L(idx), S(idx), A(idx), ...
            ismember(L(idx),CATS), ismember(S(idx),CATS), ismember(A(idx),CATS));
    end
    
   
    invalid_data = table(L(~isValid), S(~isValid), A(~isValid), G(~isValid), ...
        'VariableNames', {'PathologistA', 'PathologistB', 'PathologistC', 'Grade'});
    writetable(invalid_data, 'INVALID_ROWS_DEBUG.csv');
    fprintf('\n  Exported invalid rows to: INVALID_ROWS_DEBUG.csv\n\n');
end
L=L(isValid); S=S(isValid); A=A(isValid); G=G(isValid);
T = table(L,S,A,G,'VariableNames',{'PathologistA','PathologistB','PathologistC','Grade'});
nItems = height(T);
fprintf('N=%d items after cleaning.\n', nItems);


Dist = catdist(T, CATS);


[wK12,K12] = ckappa(T.PathologistA, T.PathologistB, CATS, true);
[wK13,K13] = ckappa(T.PathologistA, T.PathologistC, CATS, true);
[wK23,K23] = ckappa(T.PathologistB, T.PathologistC, CATS, true);
M = count_rows(T, CATS);
[Fk,Flo,Fhi,Fp] = fleiss_kappa_ci(M);
[Fk_weighted, Fk_weighted_p, Fk_weighted_se] = fleiss_kappa_weighted(M);


fprintf('  Computing bootstrap CI for weighted kappa...');
[Fk_weighted_lo, Fk_weighted_hi] = bootstrap_fleiss_weighted_ci(T, CATS, 2000);
fprintf(' done.\n');

m = 3;
P_i   = (sum(M.^2,2) - m) ./ (m*(m-1));
Pbar  = mean(P_i);
exact_all3 = sum(any(M==3,2));
exact_pct  = exact_all3 / size(M,1);
two_of_three = sum(any(M==2,2) & ~any(M==3,2));
fprintf('\n=== Inter-rater agreement (5-level) ===\n');
fprintf('Cohen kappa (unweighted / linear-weighted):\n');
fprintf('  Pathologist A-B: %.3f / %.3f\n',K12,wK12);
fprintf('  Pathologist A-C: %.3f / %.3f\n',K13,wK13);
fprintf('  Pathologist B-C: %.3f / %.3f\n',K23,wK23);
fprintf('Fleiss kappa (unweighted): %.3f (95%% CI %.3f-%.3f), p=%.3g\n',Fk,Flo,Fhi,Fp);
fprintf('Fleiss kappa (linear-weighted, ordinal): %.3f (95%% CI %.3f-%.3f), p=%.3g\n', ...
    Fk_weighted, Fk_weighted_lo, Fk_weighted_hi, Fk_weighted_p);
fprintf('Observed agreement Pbar: %.3f (%.1f%%)\n',Pbar,100*Pbar);
fprintf('Exact 3/3: %d/%d (%.1f%%), 2/3 majority: %d/%d (%.1f%%)\n', ...
    exact_all3, size(M,1), 100*exact_pct, two_of_three, size(M,1), 100*two_of_three/size(M,1));


C12 = confusion_table(T.PathologistA, T.PathologistB, CATS);
C13 = confusion_table(T.PathologistA, T.PathologistC, CATS);
C23 = confusion_table(T.PathologistB, T.PathologistC, CATS);
figure('Color','w','Position',[100 100 1200 400]);
subplot(1,3,1); heat5(C12,CATS,'Pathologist A','Pathologist B');
subplot(1,3,2); heat5(C13,CATS,'Pathologist A','Pathologist C');
subplot(1,3,3); heat5(C23,CATS,'Pathologist B','Pathologist C');
sgtitle('Rater vs rater (5x5)', 'FontSize', 18, 'FontWeight', 'bold');


all_scores = [T.PathologistA; T.PathologistB; T.PathologistC];
pooled_counts = arrayfun(@(c) sum(all_scores==c), CATS);


figure('Color','w','Position',[100 100 600 500]);
b1 = bar(100*pooled_counts/sum(pooled_counts), 0.7);
b1.FaceColor = 'flat';
b1.CData = SCORE_COLORS;
set(gca, 'XTickLabel', CATS, 'FontSize', 14, 'LineWidth', 1.5);
xlabel('PM TRPV4 score', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Proportion of ratings (%%)', 'FontSize', 16, 'FontWeight', 'bold');
title('Overall distribution (percent)', 'FontSize', 18, 'FontWeight', 'bold');
box on; grid on; grid minor;
ylim([0 max(100*pooled_counts/sum(pooled_counts))*1.15]);


[~,pct_merged] = grouped_by_grade(T,CATS,GRADE3);
figure('Color','w','Position',[100 100 750 500]);
hb = bar(pct_merged','grouped'); 
for ii=1:numel(hb)
    hb(ii).BarWidth = 0.7;
    hb(ii).FaceColor = SCORE_COLORS(ii,:);
end
set(gca, 'XTickLabel', GRADE3, 'FontSize', 14, 'LineWidth', 1.5);
xlabel('Merged grade', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Proportion of ratings (%%)', 'FontSize', 16, 'FontWeight', 'bold');
legend(CATS, 'Location', 'northeastoutside', 'FontSize', 13);
title('Score distributions by grade', 'FontSize', 18, 'FontWeight', 'bold');
box on; grid on;
ylim([0 max(pct_merged(:))*1.15]);


pm_bin = dichotomize_pm(T, CATS);
[bA, bB, bC] = deal(pm_from_rater(T.PathologistA), pm_from_rater(T.PathologistB), pm_from_rater(T.PathologistC));
[K12b,~] = kappa_binary(bA,bB); [K13b,~] = kappa_binary(bA,bC); [K23b,~] = kappa_binary(bB,bC);
Mb3 = [bA bB bC]; Fb = fleiss_kappa_binary(Mb3);
fprintf('\n=== Binary PM+/PM- ===\n');
fprintf('Cohen kappa (binary): A-B %.3f | A-C %.3f | B-C %.3f\n',K12b,K13b,K23b);
fprintf('Fleiss kappa (binary): %.3f\n',Fb);
[grpRates, grpCI, grpN] = pm_rate_by_grade(pm_bin,G,GRADE3);
fprintf('PM+ by grade (Wilson 95%% CI):\n');
for i=1:numel(GRADE3)
    fprintf('  %s: %.1f%% (%.1f-%.1f%%), n=%d\n', GRADE3(i), 100*grpRates(i), 100*grpCI(i,1), 100*grpCI(i,2), grpN(i));
end


figure('Color','w','Position',[100 100 600 500]);
hold on; x = 1:3;


for i=1:3
    idx = G==GRADE3(i);
    pm_vals = double(pm_bin(idx)) * 100;
    n_pts = numel(pm_vals);
    if n_pts > 0
        jitter = (rand(n_pts,1) - 0.5) * 0.2;
        scatter(x(i) + jitter, pm_vals, 50, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 0.4);
    end
end


lower_error = 100*grpRates - 100*grpCI(:,1);
upper_error = 100*grpCI(:,2) - 100*grpRates;

errorbar(x, 100*grpRates, lower_error, upper_error, ...
    'o', 'MarkerSize', 14, 'LineWidth', 3, 'Color', [0.85 0.2 0.15], ...
    'MarkerFaceColor', [0.85 0.2 0.15], 'CapSize', 14);

set(gca, 'XTick', x, 'XTickLabel', GRADE3, 'FontSize', 14, 'LineWidth', 1.5);
xlim([0.4 3.6]);
ylim([-5 110]);
ylabel('PM-positive rate (%%)', 'FontSize', 16, 'FontWeight', 'bold'); 
title('PM-TRPV4 positivity by grade (95%% Wilson CI)', 'FontSize', 18, 'FontWeight', 'bold'); 
grid on; box on;

[z_trend,p_trend] = cochran_armitage(pm_bin,G,GRADE3);
[rho,p_spear]     = spearman_scores_vs_grade(T);
fprintf('Cochran-Armitage trend (LG->IMG->HG): z=%.2f, p=%.3g\n', z_trend, p_trend);
fprintf('Spearman rho (ordinal score vs grade): rho=%.3f, p=%.3g\n', rho, p_spear);


writetable(T,'trpv4_scores_with_grades_cleaned.csv');
writetable(Dist,'trpv4_score_distribution.csv');
writematrix(C12,'confusion_PathologistA_vs_PathologistB.csv');
writematrix(C13,'confusion_PathologistA_vs_PathologistC.csv');
writematrix(C23,'confusion_PathologistB_vs_PathologistC.csv');
Tbin = table(G, pm_bin, 'VariableNames',{'Grade','PM_Positive'});
writetable(Tbin,'pm_binary_by_item.csv');
fprintf('\nSaved CSV files\n');

% helper
function s = normTok(s)
    s = string(s);
    % Replace missing/NaN with empty string
    s(ismissing(s)) = "";
    % Replace various dash types
    s = replace(s, char(8211), "-");  % en dash
    s = replace(s, char(8212), "-");  % em dash
    s = replace(s, char(8722), "-");  % minus sign
    s = replace(s, "–", "-");         % en dash (UTF-8)
    s = replace(s, "—", "-");         % em dash (UTF-8)
    % Trim whitespace
    s = strtrim(s);
    
    % *** KEY FIX: Convert 0.5 (Excel-friendly) to 0-1 ***
    s = replace(s, "0.5", "0-1");
    
    % Handle Excel formula issues: if "0-1" was interpreted as -1
    s(s == "-1") = "0-1";
    
    % Excel might convert 0-1 to 0.1 or other formats
    s = replace(s, "0.1", "0-1");
    s = replace(s, "0 -1", "0-1");
    s = replace(s, "0- 1", "0-1");
    s = replace(s, "0 - 1", "0-1");
    
    % Handle if Excel stored as "Jan-00" or date format
    s = replace(s, "Jan-00", "0-1");
    s = replace(s, "01-Jan", "0-1");
    
    % Remove any "+" suffixes
    s = replace(s, "0-1+", "0-1");
    s = replace(s, "1+", "1");
    s = replace(s, "2+", "2");
    s = replace(s, "3+", "3");
    
    
    s = lower(s);
end

function str = num2str_custom(x)
    % Custom numeric to string converter
    if isnan(x) || ismissing(x)
        str = "";
    elseif x == -1
        str = "0-1";  
    else
        str = num2str(x);
    end
end

function [G, keep] = merge_grade_from_pathology(ptext)
    G = strings(size(ptext)); keep = false(size(ptext));
    for i=1:numel(ptext)
        t = upper(string(ptext(i)));
        hasDCIS = contains(t,"DCIS");
        hasLG   = contains(t,"LG") | contains(t,"LOW");
        hasIMG  = contains(t,"IMG") | contains(t,"INTERMEDIATE");
        hasHG   = contains(t,"HG") | contains(t,"HIGH");
        
        if hasDCIS
            if hasHG
                G(i) = "HG"; keep(i)=true;
            elseif hasIMG
                G(i) = "IMG"; keep(i)=true;
            elseif hasLG
                G(i) = "LG";  keep(i)=true;
            else
                G(i) = "IMG"; keep(i)=true;
            end
        else
            G(i) = "";    keep(i)=false;
        end
    end
end

function Dist = catdist(T, CATS)
    f = @(s) arrayfun(@(c) sum(s==c), CATS);
    Dist = table(CATS', f(T.PathologistA)', f(T.PathologistB)', f(T.PathologistC)', ...
        'VariableNames',{'Category','PathologistA','PathologistB','PathologistC'});
end

function M = count_rows(T, CATS)
    n = height(T); k = numel(CATS); M = zeros(n,k);
    for i=1:n
        v = [T.PathologistA(i), T.PathologistB(i), T.PathologistC(i)];
        for t = v
            idx = find(CATS==t,1);
            if ~isempty(idx), M(i,idx) = M(i,idx)+1; end
        end
    end
end

function [wK,K] = ckappa(a,b,CATS,weighted)
    if nargin<4, weighted=false; end
    cats = CATS; m = numel(cats); M = zeros(m);
    for i=1:m
        for j=1:m
            M(i,j) = sum(a==cats(i) & b==cats(j));
        end
    end
    N = sum(M,'all');
    if N==0, K=NaN; wK=NaN; return; end
    po = trace(M)/N;
    r = sum(M,2)/N; c = sum(M,1,'omitnan')/N;
    pe = r'*c';
    if abs(1-pe) < eps
        K = NaN;
    else
        K = (po-pe)/(1-pe);
    end

    if ~weighted, wK = K; return; end
    W = zeros(m);
    for ii=1:m
        for jj=1:m
            W(ii,jj) = 1 - abs(ii-jj)/(m-1);
        end
    end
    Pobs = M/N; Pexp = (r)*(c);
    pw = sum(W.*Pobs,'all'); pew = sum(W.*Pexp,'all');
    if abs(1-pew) < eps
        wK = NaN;
    else
        wK = (pw - pew) / (1 - pew);
    end
end

function [k, k_lo, k_hi, pval] = fleiss_kappa_ci(M)
    n = size(M,1);
    if n == 0
       k=NaN; k_lo=NaN; k_hi=NaN; pval=NaN;
       return;
    end
    m = sum(M(1,:));
    if any(sum(M,2) ~= m)
        error('Each row must sum to the same number of raters.');
    end

    p = sum(M,1) / (n*m);
    P_i = (sum(M.^2, 2) - m) ./ (max(eps, m * (m-1)));
    Pbar = mean(P_i);
    Pexp = sum(p.^2);

    denominator = 1 - Pexp;
    if abs(denominator) < eps
        if abs(Pbar - Pexp) < eps
           k = 1;
        else
           k = NaN;
        end
        se = NaN; k_lo = NaN; k_hi = NaN; pval = NaN;
    else
        k = (Pbar - Pexp) / denominator;
        s = sum(p .* (1 - p));
        if s <= 0
            varK_corrected = NaN;
        else
             varK_corrected = (2 / (n * m * (m - 1))) * s / (denominator^2);
        end

        if isnan(varK_corrected) || varK_corrected < 0
            se = NaN;
        else
            se = sqrt(varK_corrected);
        end

        if isnan(se) || se == 0
            k_lo = k; k_hi = k;
            pval = NaN;
        else
            z_alpha = 1.96;
            k_lo = k - z_alpha * se;
            k_hi = k + z_alpha * se;
            zstat = k / se;
            pval = erfc(abs(zstat) / sqrt(2));
        end
    end
end

function [k_weighted, pval, se] = fleiss_kappa_weighted(M)
    
    
    n = size(M, 1);  
    k = size(M, 2);  
    m = sum(M(1, :));  
    
    if any(sum(M, 2) ~= m)
        error('All subjects must have the same number of raters');
    end
    
   
    weights = zeros(k, k);
    for i = 1:k
        for j = 1:k
            weights(i, j) = 1 - abs(i - j) / (k - 1);
        end
    end
    
    
    p = sum(M, 1) / (n * m);
    
    % Observed weighted agreement
    P_i_weighted = zeros(n, 1);
    for subj = 1:n
        for c1 = 1:k
            for c2 = 1:k
                P_i_weighted(subj) = P_i_weighted(subj) + weights(c1, c2) * M(subj, c1) * M(subj, c2);
            end
        end
        P_i_weighted(subj) = (P_i_weighted(subj) - m) / (m * (m - 1));
    end
    P_obs_weighted = mean(P_i_weighted);
    
   
    P_exp_weighted = 0;
    for c1 = 1:k
        for c2 = 1:k
            P_exp_weighted = P_exp_weighted + weights(c1, c2) * p(c1) * p(c2);
        end
    end
    
   
    denominator = 1 - P_exp_weighted;
    if abs(denominator) < eps
        if abs(P_obs_weighted - P_exp_weighted) < eps
            k_weighted = 1;
        else
            k_weighted = NaN;
        end
        se = NaN;
        pval = NaN;
        return;
    else
        k_weighted = (P_obs_weighted - P_exp_weighted) / denominator;
    end
    
    
    var_P_i = var(P_i_weighted);
    
    sum_term = 0;
    for c = 1:k
        sum_term = sum_term + p(c) * (1 - p(c));
    end
    
    var_kappa = (2 / (n * m * (m - 1))) * sum_term / (denominator^2);
    
  
    if var_P_i > 0
        var_kappa_improved = var_P_i / (n * denominator^2);
        se = sqrt(var_kappa_improved);
    else
        se = sqrt(max(0, var_kappa));
    end
    
    if isnan(se) || se <= 0
        pval = NaN;
    else
        z_stat = k_weighted / se;
        pval = 2 * (1 - normcdf(abs(z_stat)));  
    end
end

function C = confusion_table(a,b,CATS)
    m = numel(CATS); C = zeros(m);
    for i=1:m
        for j=1:m
            C(i,j) = sum(a==CATS(i) & b==CATS(j));
        end
    end
end

function heat5(C,CATS,xlab,ylab)
    imagesc(C); axis square; 
    cb = colorbar;
    cb.FontSize = 12;
    set(gca,'XTick',1:numel(CATS),'XTickLabel',CATS, ...
            'YTick',1:numel(CATS),'YTickLabel',CATS,'FontSize',13, 'LineWidth', 1.5);
    xlabel(xlab, 'FontSize', 15, 'FontWeight', 'bold'); 
    ylabel(ylab, 'FontSize', 15, 'FontWeight', 'bold'); 
    title('Confusion (counts)', 'FontSize', 16, 'FontWeight', 'bold');
end

function [counts_merged,pct_merged] = grouped_by_grade(T,CATS,GR)
    Kc = numel(CATS); Gc = numel(GR); counts_merged = zeros(Kc,Gc);
    for gix = 1:Gc
        idx = T.Grade==GR(gix);
        scores_g = [T.PathologistA(idx); T.PathologistB(idx); T.PathologistC(idx)];
        counts_merged(:,gix) = arrayfun(@(c) sum(scores_g==c), CATS);
    end
    col_sums = sum(counts_merged,1);
    pct_merged = 100 * counts_merged ./ max(col_sums,eps);
end

function y = pm_from_rater(col)
    positive_scores = ["1", "2", "3"];
    y = ismember(col, positive_scores);
end

function pm = dichotomize_pm(T, CATS)
    n = height(T); pm = false(n,1);
    for i=1:n
        v = [T.PathologistA(i), T.PathologistB(i), T.PathologistC(i)];
        if any(ismember(v, ["1","2","3"]))
            pm(i) = true;
        elseif all(v == "0-1")
             pm(i) = true;
        else
             pm(i) = false;
        end
    end
end

function [kappa, tab] = kappa_binary(a,b)
    a = logical(a); b = logical(b);
    tab = confusionmat(a, b);
    N = sum(tab,'all');
    if N == 0, kappa = NaN; return; end
    po = (tab(1,1)+tab(2,2))/N;
    p_a_true = sum(a)/numel(a); p_a_false = 1-p_a_true;
    p_b_true = sum(b)/numel(b); p_b_false = 1-p_b_true;
    pe = (p_a_true * p_b_true) + (p_a_false * p_b_false);
    if abs(1-pe) < eps
        kappa = NaN;
    else
        kappa = (po-pe)/(1-pe);
    end
end

function Fb = fleiss_kappa_binary(M3)
    N = size(M3,1);
    if N == 0, Fb = NaN; return; end
    Mb = [sum(~M3,2) sum(M3,2)];
    m = 3;
    if any(sum(Mb,2) ~= m)
        error('Binary rating counts do not sum to raters.');
    end
    p = sum(Mb,1)/(N*m);
    P_i = (sum(Mb.^2,2)-m)./(max(eps, m*(m-1)));
    Pbar = mean(P_i);
    Pexp = sum(p.^2);
    denominator = 1 - Pexp;
     if abs(denominator) < eps
        if abs(Pbar - Pexp) < eps
            Fb = 1;
        else
            Fb = NaN;
        end
    else
        Fb = (Pbar - Pexp)/denominator;
    end
end

function [rates, ci, Ns] = pm_rate_by_grade(pm_bin,G,GR)
    rates=zeros(numel(GR),1); ci=zeros(numel(GR),2); Ns=zeros(numel(GR),1);
    for i=1:numel(GR)
        idx = G==GR(i);
        Ns(i)=sum(idx);
        k=sum(pm_bin(idx));
        n=Ns(i);
        if n==0
            rates(i)=NaN; ci(i,:)=[NaN NaN];
        else
            rates(i)=k/n; [lo,hi]=wilson_ci(k,n,0.05); ci(i,:)=[lo hi];
        end
    end
end

function [lo,hi] = wilson_ci(k,n,alpha)
    if nargin<3, alpha=0.05; end
    if n == 0, lo=NaN; hi=NaN; return; end
    z = norminv_local(1-alpha/2);
    p = k/n; denom = 1 + z^2/n;
    center = (p + z^2/(2*n))/denom;
    term_under_sqrt = (p*(1-p)/n + z^2/(4*n^2));
    if term_under_sqrt < 0, term_under_sqrt = 0; end
    half   = (z / denom) * sqrt(term_under_sqrt);
    lo = max(0, center - half); hi = min(1, center + half);
end

function [z,p] = cochran_armitage(pm_bin,G,GR)
    x=(0:numel(GR)-1)';
    n=zeros(numel(GR),1); y=n;
    for i=1:numel(GR)
        idx=G==GR(i);
        n(i)=sum(idx);
        y(i)=sum(pm_bin(idx));
    end
    N=sum(n);
    if N==0, z=NaN; p=NaN; return; end
    p_hat=sum(y)/N;
    if p_hat == 0 || p_hat == 1
        z=0; p=1; return;
    end
    x_bar = sum(x.*n)/N;
    w=x-x_bar;
    numerator=sum(w.*(y - n*p_hat));
    denominator_sq = p_hat*(1-p_hat) * sum(w.^2 .* n);
    if denominator_sq <= 0
        z = 0; p = 1;
    else
        den=sqrt(denominator_sq);
        z = numerator/den;
        p = erfc(abs(z)/sqrt(2));
    end
end

function val = norminv_local(p)
    if p<=0 || p>=1, error('p must be in (0,1)'); end
    a = [ -3.969683028665376e+01,  2.209460984245205e+02,...
          -2.759285104469687e+02,  1.383577518672690e+02,...
          -3.066479806614716e+01,  2.506628277459239e+00 ];
    b = [ -5.447609879822406e+01,  1.615858368580409e+02,...
          -1.556989798598866e+02,  6.680131188771972e+01,...
          -1.328068155288572e+01 ];
    c = [ -7.784894002430293e-03, -3.223964580411365e-01,...
          -2.400758277161838e+00, -2.549732539343734e+00,...
           4.374664141464968e+00,  2.938163982698783e+00 ];
    d = [ 7.784695709041462e-03,  3.224671290700398e-01,...
          2.445134137142996e+00,  3.754408661907416e+00 ];
    plow = 0.02425; phigh = 1 - plow;
    if p < plow
        q = sqrt(-2*log(p));
        val = polyval(c,q)/ (polyval(d,q)*q+1);
    elseif phigh < p
        q = sqrt(-2*log(1-p));
        val = -polyval(c,q)/ (polyval(d,q)*q+1);
    else
        q = p - 0.5; r = q*q;
        val = q * polyval(a,r) / (polyval(b,r)*r+1);
    end
end

function [rho, pval] = spearman_scores_vs_grade(T)
    score_map = containers.Map(["0", "0-1", "1", "2", "3"], [0, 0.5, 1, 2, 3]);
    grade_map = containers.Map(["LG", "IMG", "HG"], [1, 2, 3]);

    n = height(T);
    numeric_scores = zeros(n, 3);

    try
        numeric_scores(:, 1) = cell2mat(values(score_map, cellstr(T.PathologistA)));
        numeric_scores(:, 2) = cell2mat(values(score_map, cellstr(T.PathologistB)));
        numeric_scores(:, 3) = cell2mat(values(score_map, cellstr(T.PathologistC)));
    catch ME
        error('Error converting scores: %s', ME.message);
    end

    avg_numeric_score = mean(numeric_scores, 2);

    try
        numeric_grade = cell2mat(values(grade_map, cellstr(T.Grade)));
    catch ME
        error('Error converting grades: %s', ME.message);
    end

    if numel(unique(avg_numeric_score(~isnan(avg_numeric_score)))) < 2 || numel(unique(numeric_grade(~isnan(numeric_grade)))) < 2
        rho = NaN; pval = NaN;
    else
       [rho, pval] = corr(avg_numeric_score, numeric_grade, 'Type', 'Spearman', 'rows', 'complete');
    end
end

function [ci_lo, ci_hi] = bootstrap_fleiss_weighted_ci(T, CATS, n_boot)
    
    if nargin < 3, n_boot = 2000; end
    
    n_items = height(T);
    kappa_boot = zeros(n_boot, 1);
    
    M_orig = count_rows(T, CATS);
    kappa_orig = fleiss_kappa_weighted(M_orig);
   
    for b = 1:n_boot
        
        boot_idx = randi(n_items, n_items, 1);
        T_boot = T(boot_idx, :);
        
        
        M_boot = count_rows(T_boot, CATS);
        kappa_boot(b) = fleiss_kappa_weighted(M_boot);
    end
    
   
    ci_lo = prctile(kappa_boot, 2.5);
    ci_hi = prctile(kappa_boot, 97.5);
end

end 