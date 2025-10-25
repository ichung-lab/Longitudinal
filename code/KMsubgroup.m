function Figure4_subgroup_analysis(fname)
% Supple Fig 2S: Subgroup analysis

if nargin<1 || isempty(fname), fname='score4.xlsx'; end
if ~isfile(fname), error('Input not found: %s',fname); end

%% read columns
T = readtable(fname,'PreserveVariableNames',true);
cols = T.Properties.VariableNames;

event_col = pickcol(cols, ["Progression","Event","Status"]);
time_col  = pickcol(cols, ["time","followup","follow_up","fu","Time"]);
pm_col    = pickcol(cols, ["PM","PM_TRPV4","TRPV4_PM","PM-TRPV4","PMTRPV4"]);
grade_col = pickcol(cols, ["HG","HighGrade","Grade","HistologicGrade"]);

% Extract + coerce
events = to01(T.(event_col));
traw   = double(T.(time_col));
PM     = to01(T.(pm_col));
HG     = to01(T.(grade_col));

% Auto units: days -> years if values look big
if max(traw) > 50, time_y = traw/365.25; else, time_y = traw; end
N = numel(time_y);

%% splitting subgroups
idx_lg = (HG==0);
idx_hg = (HG==1);

%LG+IMG
tl = time_y(idx_lg); el = events(idx_lg); pml = PM(idx_lg);
[iP_lg,iN_lg] = deal(pml==1, pml==0);
[Spos_lg,Tpos_lg,SEpos_lg] = km(tl(iP_lg), el(iP_lg));
[Sneg_lg,Tneg_lg,SEneg_lg] = km(tl(iN_lg), el(iN_lg));
p_lr_lg = logrank(tl(iP_lg), el(iP_lg), tl(iN_lg), el(iN_lg));
[hr_lg, lo_lg, hi_lg, p_cx_lg] = cox_uni(tl, el, pml);

% HG
th = time_y(idx_hg); eh = events(idx_hg); pmh = PM(idx_hg);
[iP_hg,iN_hg] = deal(pmh==1, pmh==0);
[Spos_hg,Tpos_hg,SEpos_hg] = km(th(iP_hg), eh(iP_hg));
[Sneg_hg,Tneg_hg,SEneg_hg] = km(th(iN_hg), eh(iN_hg));
p_lr_hg = logrank(th(iP_hg), eh(iP_hg), th(iN_hg), eh(iN_hg));
[hr_hg, lo_hg, hi_hg, p_cx_hg] = cox_uni(th, eh, pmh);

% Cos: grade/PM interaction
[p_interaction, ~, ~] = cox_interaction(time_y, events, PM, HG);

%% 2x2
prog_tp_lg = sum(PM==1 & HG==0 & events==1);
tot_tp_lg  = sum(PM==1 & HG==0);
prog_tn_lg = sum(PM==0 & HG==0 & events==1);
tot_tn_lg  = sum(PM==0 & HG==0);

prog_tp_hg = sum(PM==1 & HG==1 & events==1);
tot_tp_hg  = sum(PM==1 & HG==1);
prog_tn_hg = sum(PM==0 & HG==1 & events==1);
tot_tn_hg  = sum(PM==0 & HG==1);

rate = 100*[prog_tp_lg/max(tot_tp_lg,1), prog_tn_lg/max(tot_tn_lg,1), ...
            prog_tp_hg/max(tot_tp_hg,1), prog_tn_hg/max(tot_tn_hg,1)];

%% Fig
colors.tp=[0.75,0.15,0.18]; colors.tn=[0.13,0.55,0.27];
colors.lg=[0.95,0.65,0.20]; colors.hg=[0.10,0.35,0.80];
tmax = max(time_y)*1.05;

figure('Color','w','Position',[50 50 1600 900],'Name','Figure 4 - Subgroup Analysis');

posA=[0.07,0.55,0.40,0.38]; posB=[0.55,0.55,0.40,0.38];
posC=[0.07,0.10,0.40,0.35]; posD=[0.55,0.10,0.40,0.35];

% A. lower g
subplot('Position',posA); hold on;
plot([0 tmax],[0.5 0.5],'k--','LineWidth',1,'HandleVisibility','off');
plot_km_curve_extended(Tpos_lg,Spos_lg,SEpos_lg,colors.tp,'PM-TRPV4+',tmax);
plot_km_curve_extended(Tneg_lg,Sneg_lg,SEneg_lg,colors.tn,'PM-TRPV4-',tmax);
xlim([0 tmax]); ylim([0 1.05]); set(gca,'XTick',0:2:ceil(tmax),'FontSize',13,'LineWidth',1.5);
xlabel('Time (years)','FontSize',14,'FontWeight','bold'); ylabel('IDC-free survival','FontSize',14,'FontWeight','bold');
legend('Location','southwest','FontSize',13); title('A. Lower Grade DCIS','FontWeight','bold','FontSize',15);
text(0.98,0.35,sprintf(['n = %d\nLog-rank p = %.4f\nCox HR = %.2f (%.2f–%.2f)\nCox p = %.4f'], ...
    sum(idx_lg), p_lr_lg, hr_lg, lo_lg, hi_lg, p_cx_lg), ...
    'Units','normalized','FontSize',12,'HorizontalAlignment','right','VerticalAlignment','top', ...
    'BackgroundColor','white','EdgeColor','black','LineWidth',1,'Margin',5);
box on; add_risk_table(gca, tl(iP_lg), tl(iN_lg), tmax, {'PM-TRPV4+','PM-TRPV4-'}, [colors.tp; colors.tn], 13);

% B. HG
subplot('Position',posB); hold on;
plot([0 tmax],[0.5 0.5],'k--','LineWidth',1,'HandleVisibility','off');
plot_km_curve_extended(Tpos_hg,Spos_hg,SEpos_hg,colors.tp,'PM-TRPV4+',tmax);
plot_km_curve_extended(Tneg_hg,Sneg_hg,SEneg_hg,colors.tn,'PM-TRPV4-',tmax);
xlim([0 tmax]); ylim([0 1.05]); set(gca,'XTick',0:2:ceil(tmax),'FontSize',13,'LineWidth',1.5);
xlabel('Time (years)','FontSize',14,'FontWeight','bold'); ylabel('IDC-free survival','FontSize',14,'FontWeight','bold');
legend('Location','southwest','FontSize',13); title('B. High Grade DCIS','FontWeight','bold','FontSize',15);
text(0.98,0.35,sprintf(['n = %d\nLog-rank p = %.4f\nCox HR = %.2f (%.2f–%.2f)\nCox p = %.4f'], ...
    sum(idx_hg), p_lr_hg, hr_hg, lo_hg, hi_hg, p_cx_hg), ...
    'Units','normalized','FontSize',12,'HorizontalAlignment','right','VerticalAlignment','top', ...
    'BackgroundColor','white','EdgeColor','black','LineWidth',1,'Margin',5);
box on; add_risk_table(gca, th(iP_hg), th(iN_hg), tmax, {'PM-TRPV4+','PM-TRPV4-'}, [colors.tp; colors.tn], 13);

% C. Forest 
subplot('Position',posC); hold on;
ypos=[1 2]; xline(1,'k--','LineWidth',1.5);
errorbar(hr_hg, ypos(1), hr_hg-lo_hg, hi_hg-hr_hg, 'horizontal','o','MarkerSize',12,'LineWidth',2.5,'Color','k','MarkerFaceColor',colors.hg);
errorbar(hr_lg, ypos(2), hr_lg-lo_lg, hi_lg-hr_lg, 'horizontal','o','MarkerSize',12,'LineWidth',2.5,'Color','k','MarkerFaceColor',colors.lg);
set(gca,'YTick',ypos,'YTickLabel',{'High Grade','Lower Grade'},'FontSize',13,'LineWidth',1.5,'XScale','log');
xlim([0.1 100]); set(gca,'XTick',[0.1 1 10 100]); ylim([0.5 2.5]);
xlabel('Hazard Ratio (log scale)','FontSize',14,'FontWeight','bold'); title('C. PM-TRPV4 Effect by Grade','FontWeight','bold','FontSize',15);
txt = sprintf(['Lower Grade  HR %.2f (%.2f–%.2f), Cox p=%.4f\n' ...
               'High Grade   HR %.2f (%.2f–%.2f), Cox p=%.4f\n' ...
               'Interaction p=%.4f'], hr_lg,lo_lg,hi_lg,p_cx_lg, hr_hg,lo_hg,hi_hg,p_cx_hg, p_interaction);
annotation('textbox',[0.095 0.16 0.28 0.16],'String',txt,'FontSize',11,'EdgeColor','k','BackgroundColor','w');

grid on; box on;

% D. Progression-rate 
subplot('Position',posD); hold on;
xpos=[1 2 4 5]; vals=rate;
bar_colors=[colors.tp; colors.tn; colors.tp; colors.tn];
for i=1:4, bar(xpos(i), vals(i), 0.7, 'FaceColor',bar_colors(i,:), 'EdgeColor','k','LineWidth',1.5); end
counts=[prog_tp_lg, prog_tn_lg, prog_tp_hg, prog_tn_hg];
totals=[tot_tp_lg, tot_tn_lg, tot_tp_hg, tot_tn_hg];
for i=1:4
    text(xpos(i), vals(i)+5, sprintf('%d/%d',counts(i),totals(i)), 'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
    text(xpos(i), max(vals(i)*0.5,3), sprintf('%.1f%%',vals(i)), 'HorizontalAlignment','center','FontSize',12,'FontWeight','bold','Color','w');
end
text(1.5,-8,'Lower Grade','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
text(4.5,-8,'High Grade','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
set(gca,'XTick',xpos,'XTickLabel',{'TRPV4+','TRPV4-','TRPV4+','TRPV4-'},'FontSize',12,'LineWidth',1.5);
ylim([0 max([vals 80])+15]); ylabel('Progression Rate (%)','FontSize',14,'FontWeight','bold');
title('D. Progression Rates by PM-TRPV4 and Grade','FontWeight','bold','FontSize',15);
grid on; box on;

sgtitle(sprintf('PM-TRPV4 Effect within Grade Subgroups (N=%d)',N), 'FontSize',18,'FontWeight','bold');

%% helpers
function name = pickcol(vnames, keys)
    for k=1:numel(keys)
        j = find(strcmpi(vnames,keys{k}) | contains(vnames,keys{k},'IgnoreCase',true),1);
        if ~isempty(j), name=vnames{j}; return; end
    end
    error('Missing column matching: %s', strjoin(string(keys),', '));
end

function y = to01(x)
    if isnumeric(x)
        y = double(x~=0 & ~isnan(x));
    else
        s = upper(string(x));
        y = double(s=="1" | s=="YES" | s=="Y" | s=="+" | s=="POS" | s=="TRUE");
    end
end

function [S,T,SE] = km(t,e)
    t=t(:); e=e(:);
    if isempty(t), T=[]; S=[]; SE=[]; return; end
    T = sort(unique(t(e==1))); if isempty(T), T = sort(unique(t)); end
    n=numel(T); S=ones(n,1); SE=zeros(n,1); v=0;
    for i=1:n
        ti=T(i); at=sum(t>=ti); d=sum(t==ti & e==1);
        if i==1, S(i) = 1 - d/max(at,1); else, S(i)=S(i-1)*(1 - d/max(at,1)); end
        if d>0 && at>d, v=v + d/(at*(at-d)); end
        SE(i) = S(i)*sqrt(max(v,0));
    end
end

function p = logrank(t1,e1,t0,e0)
    t=[t1(:);t0(:)]; e=[e1(:);e0(:)]; g=[ones(numel(t1),1);zeros(numel(t0),1)];
    tu=sort(unique(t(e==1))); OE=0; V=0;
    for i=1:numel(tu)
        ti=tu(i); risk=(t>=ti); n=sum(risk); if n<=1, continue; end
        n1=sum(risk & g==1); d=sum(t==ti & e==1); if d==0, continue; end
        d1=sum(t1==ti & e1==1); E1=d*(n1/n); Vi=(n1*(n-n1)*d*(n-d))/((n^2)*(n-1));
        OE=OE+(d1-E1); V=V+Vi;
    end
    p = (V<=0)*1 + (V>0)*(1-chi2cdf((OE^2)/V,1));
end

function [hr,lo,hi,p] = cox_uni(t,e,x)
    try
        [b,~,stats] = coxphfit(x, t, 'Censoring',1-e, 'Ties','efron');
        hr=exp(b); se=stats.se; lo=exp(b-1.96*se); hi=exp(b+1.96*se);
        z=b./se; p=2*(1-normcdf(abs(z)));
        hr=hr(1); lo=lo(1); hi=hi(1); p=p(1);
    catch
        % Fallback: Haldane–Anscombe rate ratio
        i1=(x==1); i0=(x==0);
        e1=sum(e(i1)); e0=sum(e(i0)); py1=sum(t(i1)); py0=sum(t(i0));
        r1=(e1+0.5)/max(py1,eps); r0=(e0+0.5)/max(py0,eps);
        hr=r1/r0; se=sqrt(1/(e1+0.5)+1/(e0+0.5));
        lo=exp(log(hr)-1.96*se); hi=exp(log(hr)+1.96*se);
        z=log(hr)/se; p=2*(1-normcdf(abs(z)));
    end
end

function [p_int,b,se] = cox_interaction(t,e,pm,hg)
    X=[pm(:), hg(:), pm(:).*hg(:)];
    mask=all(~isnan(X),2) & ~isnan(t(:)) & ~isnan(e(:));
    try
        [b,~,stats] = coxphfit(X(mask,:), t(mask), 'Censoring',1-e(mask), 'Ties','efron');
        se=stats.se(:); z=b(3)/se(3); p_int=2*(1-normcdf(abs(z)));
    catch
        p_int=NaN; b=[NaN;NaN;NaN]; se=[NaN;NaN;NaN];
    end
end

function plot_km_curve_extended(T,S,SE,color,label,tmax)
    if isempty(T), return; end
    if T(end) < tmax, T=[T; tmax]; S=[S; S(end)]; SE=[SE; SE(end)]; end
    Tpl=[0; reshape([T';T'],[],1)]; Tpl=Tpl(1:end-1);
    Spl=[1; reshape([S';S'],[],1)]; Spl=Spl(1:end-1);
    SEpl=[0; reshape([SE';SE'],[],1)]; SEpl=SEpl(1:end-1);
    up=min(1, Spl+1.96*SEpl); lo=max(0, Spl-1.96*SEpl);
    fill([Tpl; flipud(Tpl)], [up; flipud(lo)], color, 'FaceAlpha',0.2, 'EdgeColor','none','HandleVisibility','off');
    plot(Tpl,Spl,'Color',color,'LineWidth',3,'DisplayName',label);
end

function add_risk_table(ax,t1,t0,tmax,labels,colors,fs)
    if nargin<7, fs=13; end
    tps=0:2:floor(tmax);
    n1=arrayfun(@(tp) sum(t1>=tp), tps);
    n0=arrayfun(@(tp) sum(t0>=tp), tps);
    yl=ylim(ax); yr=yl(2)-yl(1); y0=yl(1)-0.12*yr; dy=0.06*yr;
    text(ax, -tmax*0.08, y0, 'No. at risk','FontSize',fs,'FontWeight','bold','HorizontalAlignment','right');
    for i=1:numel(tps)
        text(ax, tps(i), y0-dy,  sprintf('%d',n1(i)),'Color',colors(1,:),'FontSize',fs,'HorizontalAlignment','center');
        text(ax, tps(i), y0-2*dy,sprintf('%d',n0(i)),'Color',colors(2,:),'FontSize',fs,'HorizontalAlignment','center');
    end
    text(ax, -tmax*0.08, y0-dy,  labels{1}, 'Color',colors(1,:),'FontSize',fs,'HorizontalAlignment','right','FontWeight','bold');
    text(ax, -tmax*0.08, y0-2*dy,labels{2}, 'Color',colors(2,:),'FontSize',fs,'HorizontalAlignment','right','FontWeight','bold');
end

end
