function laminate_viz(res)
% LAMINATE_VIZ  Generate all 21 engineering visualizations for CLT laminate analysis.
% Called automatically by laminate_analysis(). Can also be called standalone:
%   res = laminate_analysis();  laminate_viz(res);

viz_01_laminate_stack(res);
viz_02_fiber_orientation(res);
viz_03_Q_heatmap(res);
viz_04_Qbar_heatmaps(res);
viz_05_A_heatmap(res);
viz_06_B_heatmap(res);
viz_07_D_heatmap(res);
viz_08_ABD_heatmap(res);
viz_09_stress_thickness(res);
viz_10_strain_thickness(res);
viz_11_TsaiWu_bar(res);
viz_12_TsaiHill_bar(res);
viz_13_MaxStress_bar(res);
viz_14_delamination(res);
viz_15_buckling_SF(res);
viz_16_deformation(res);
viz_17_buckling_mode(res);
viz_18_polar_stiffness(res);
viz_19_failure_envelope(res);
viz_20_ply_contribution(res);
viz_21_dashboard(res);
end

% =========================================================================
% HELPER: map angle value to a consistent RGB color
% =========================================================================
function col = ang_color(ang, all_angles)
uniq = unique(all_angles);
cmap = lines(numel(uniq));
idx  = find(uniq == ang, 1);
col  = cmap(idx, :);
end

% =========================================================================
% HELPER: traffic-light color from 0-1 value
% =========================================================================
function col = risk_color(v)
if     v > 0.7; col = [0.90, 0.15, 0.15];
elseif v > 0.4; col = [1.00, 0.70, 0.00];
else;           col = [0.10, 0.72, 0.33];
end
end

% =========================================================================
% VIZ 1 — LAMINATE STACK DIAGRAM
% =========================================================================
function viz_01_laminate_stack(res)
figure('Name','1 - Laminate Stack','NumberTitle','off', ...
       'Color','w','Position',[30 50 420 520]);
hold on; axis off;
n = res.n;  t = res.mat.t;  h = res.h;

% Draw plies from top (ply n) down to bottom (ply 1)
for k = n:-1:1
    row  = n - k;           % 0 = top row
    yb   = (k-1)*t / h;
    yt   = k*t / h;
    col  = ang_color(res.theta_deg(k), res.theta_deg);
    fill([0.1 0.9 0.9 0.1], [yb yb yt yt], col, ...
         'EdgeColor','k','LineWidth',0.8,'FaceAlpha',0.80);
    text(0.50, (yb+yt)/2, sprintf('Ply %2d    %+g\xB0', k, res.theta_deg(k)), ...
         'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color','k');
end

% Labels
text(0.50, 1.06, 'TOP',    'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','Color',[0.5,0,0]);
text(0.50,-0.05, 'BOTTOM', 'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','Color',[0.5,0,0]);
line([0.1 0.9],[0,0],'Color','k','LineWidth',2);
line([0.1 0.9],[1,1],'Color','k','LineWidth',2);

% Legend for angles
uniq = unique(res.theta_deg);
for ui = 1:numel(uniq)
    col = ang_color(uniq(ui), res.theta_deg);
    patch(NaN,NaN,col,'DisplayName',sprintf('%+g\xB0',uniq(ui)));
end
legend('Location','eastoutside','FontSize',9,'Box','off');

title('Laminate Stack Diagram','FontSize',12,'FontWeight','bold');
ylim([-0.12, 1.12]);
end

% =========================================================================
% VIZ 2 — FIBER ORIENTATION (top-view per ply)
% =========================================================================
function viz_02_fiber_orientation(res)
n = res.n;
nc = ceil(sqrt(n));  nr = ceil(n/nc);
figure('Name','2 - Fiber Orientations','NumberTitle','off', ...
       'Color','w','Position',[460 50 700 500]);
for k = 1:n
    subplot(nr, nc, k);  hold on; axis equal; axis off;
    ang = res.theta_deg(k);
    c = cosd(ang);  s = sind(ang);
    col = ang_color(ang, res.theta_deg);
    rectangle('Position',[-0.55,-0.55,1.1,1.1],'EdgeColor','k','LineWidth',1.2);
    for t_off = -0.4:0.13:0.4
        x1 = -0.5*c - t_off*s;   y1 = -0.5*s + t_off*c;
        x2 =  0.5*c - t_off*s;   y2 =  0.5*s + t_off*c;
        line([x1,x2],[y1,y2],'Color',col,'LineWidth',1.8);
    end
    title(sprintf('Ply %d: %+g\xB0',k,ang),'FontSize',8,'FontWeight','bold');
    xlim([-0.65,0.65]); ylim([-0.65,0.65]);
end
sgtitle('Fiber Orientation — Top View','FontSize',12,'FontWeight','bold');
end

% =========================================================================
% VIZ 3 — Q MATRIX HEATMAP
% =========================================================================
function viz_03_Q_heatmap(res)
figure('Name','3 - Q Matrix','NumberTitle','off','Color','w','Position',[30 620 400 360]);
Q = res.Q;
imagesc(Q/1e9); colorbar; colormap(gca, parula);
set(gca,'XTick',1:3,'YTick',1:3, ...
    'XTickLabel',{'1','2','6'},'YTickLabel',{'1','2','6'},'FontSize',10);
for i=1:3; for j=1:3
    text(j,i,sprintf('%.2f',Q(i,j)/1e9), ...
         'HorizontalAlignment','center','FontSize',9,'Color','w','FontWeight','bold');
end; end
title('Reduced Stiffness [Q] (GPa)','FontSize',11,'FontWeight','bold');
xlabel('Index'); ylabel('Index'); axis square;
end

% =========================================================================
% VIZ 4 — Qbar HEATMAPS (one per unique angle)
% =========================================================================
function viz_04_Qbar_heatmaps(res)
uniq = unique(res.theta_deg, 'stable');
nu   = numel(uniq);
figure('Name','4 - Qbar Heatmaps','NumberTitle','off','Color','w', ...
       'Position',[440 620 220*nu 340]);
for ui = 1:nu
    ang = uniq(ui);
    ki  = find(res.theta_deg == ang, 1);
    Qb  = res.Qbar(:,:,ki);
    subplot(1, nu, ui);
    imagesc(Qb/1e9); colorbar; colormap(gca, parula);
    set(gca,'XTick',1:3,'YTick',1:3, ...
        'XTickLabel',{'1','2','6'},'YTickLabel',{'1','2','6'},'FontSize',8);
    for i=1:3; for j=1:3
        text(j,i,sprintf('%.1f',Qb(i,j)/1e9), ...
             'HorizontalAlignment','center','FontSize',7,'Color','w','FontWeight','bold');
    end; end
    title(sprintf('\\bar{Q}(%+g\xB0) [GPa]',ang),'FontSize',9,'FontWeight','bold');
    axis square;
end
sgtitle('Transformed Stiffness \bar{Q}','FontSize',12,'FontWeight','bold');
end

% =========================================================================
% VIZ 5 — A MATRIX HEATMAP
% =========================================================================
function viz_05_A_heatmap(res)
figure('Name','5 - A Matrix','NumberTitle','off','Color','w','Position',[30 50 400 360]);
A = res.A;
lbl = {'A_{11}','A_{12}','A_{16}';'A_{21}','A_{22}','A_{26}';'A_{61}','A_{62}','A_{66}'};
imagesc(A/1e6); colorbar; colormap(gca, hot);
for i=1:3; for j=1:3
    text(j,i,sprintf('%s\n%.0f',lbl{i,j},A(i,j)/1e6), ...
         'HorizontalAlignment','center','FontSize',8,'Color','w','FontWeight','bold');
end; end
set(gca,'XTick',1:3,'YTick',1:3,'XTickLabel',{'1','2','6'},'YTickLabel',{'1','2','6'},'FontSize',10);
title('Extensional Stiffness [A] (MN/m)','FontSize',11,'FontWeight','bold');
axis square;
end

% =========================================================================
% VIZ 6 — B MATRIX HEATMAP
% =========================================================================
function viz_06_B_heatmap(res)
figure('Name','6 - B Matrix','NumberTitle','off','Color','w','Position',[450 50 400 360]);
B = res.B;
clim_val = max(abs(B(:)));
if clim_val < 1e-10; clim_val = 1; end
imagesc(B, [-clim_val, clim_val]); colorbar; colormap(gca, coolwarm_cmap());
lbl = {'B_{11}','B_{12}','B_{16}';'B_{21}','B_{22}','B_{26}';'B_{61}','B_{62}','B_{66}'};
for i=1:3; for j=1:3
    text(j,i,sprintf('%s\n%.2e',lbl{i,j},B(i,j)), ...
         'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
end; end
set(gca,'XTick',1:3,'YTick',1:3,'XTickLabel',{'1','2','6'},'YTickLabel',{'1','2','6'},'FontSize',10);
title('Extension-Bending Coupling [B] (N)','FontSize',11,'FontWeight','bold');
if max(abs(B(:))) < 1e-8
    text(2, 1.5, 'B \approx 0   (Symmetric Laminate)', ...
         'HorizontalAlignment','center','FontSize',10,'Color',[0,0.6,0],'FontWeight','bold');
end
axis square;
end

% =========================================================================
% VIZ 7 — D MATRIX HEATMAP
% =========================================================================
function viz_07_D_heatmap(res)
figure('Name','7 - D Matrix','NumberTitle','off','Color','w','Position',[880 50 400 360]);
D = res.D;
lbl = {'D_{11}','D_{12}','D_{16}';'D_{21}','D_{22}','D_{26}';'D_{61}','D_{62}','D_{66}'};
imagesc(D); colorbar; colormap(gca, parula);
for i=1:3; for j=1:3
    text(j,i,sprintf('%s\n%.2f',lbl{i,j},D(i,j)), ...
         'HorizontalAlignment','center','FontSize',8,'Color','w','FontWeight','bold');
end; end
set(gca,'XTick',1:3,'YTick',1:3,'XTickLabel',{'1','2','6'},'YTickLabel',{'1','2','6'},'FontSize',10);
title('Bending Stiffness [D] (N\cdotm)','FontSize',11,'FontWeight','bold');
axis square;
end

% =========================================================================
% VIZ 8 — FULL ABD HEATMAP
% =========================================================================
function viz_08_ABD_heatmap(res)
figure('Name','8 - ABD Matrix','NumberTitle','off','Color','w','Position',[30 470 520 440]);
ABD = res.ABD;
imagesc(log10(abs(ABD)+1)); colorbar; colormap(gca, parula);
hold on;
line([3.5,3.5],[0.5,6.5],'Color','r','LineWidth',2.5);
line([0.5,6.5],[3.5,3.5],'Color','r','LineWidth',2.5);
row_lbl = {'1','2','6','1','2','6'};
col_lbl = {'1','2','6','1','2','6'};
prefix  = {'A','A','A','B','B','B';'A','A','A','B','B','B';'A','A','A','B','B','B'; ...
           'B','B','B','D','D','D';'B','B','B','D','D','D';'B','B','B','D','D','D'};
for i=1:6; for j=1:6
    text(j,i,sprintf('%s_{%s%s}',prefix{i,j},row_lbl{i},col_lbl{j}), ...
         'HorizontalAlignment','center','FontSize',6,'Color','w');
end; end
set(gca,'XTick',1:6,'YTick',1:6,'FontSize',8);
text(1.5,3.8,'A','FontSize',16,'Color','r','FontWeight','bold','HorizontalAlignment','center');
text(4.5,3.8,'B','FontSize',16,'Color','r','FontWeight','bold','HorizontalAlignment','center');
text(4.5,1.5,'B','FontSize',16,'Color','r','FontWeight','bold','HorizontalAlignment','center');
text(4.5,5.0,'D','FontSize',16,'Color','r','FontWeight','bold','HorizontalAlignment','center');
title('Full Laminate Stiffness [ABD]  (log scale)','FontSize',11,'FontWeight','bold');
end

% =========================================================================
% VIZ 9 — GLOBAL STRESS THROUGH THICKNESS
% =========================================================================
function viz_09_stress_thickness(res)
figure('Name','9 - Stress Through Thickness','NumberTitle','off','Color','w', ...
       'Position',[560 470 780 440]);
z_mm = res.z_pts * 1e3;
sig  = res.sig_g / 1e6;   % MPa
comp_lbl = {'\sigma_x (MPa)','\sigma_y (MPa)','\tau_{xy} (MPa)'};
col_list  = {[0.1,0.3,0.9],[0.8,0.1,0.1],[0.1,0.6,0.2]};

for c = 1:3
    subplot(1,3,c); hold on; grid on; box on;
    plot(sig(c,:), z_mm, '-o', 'Color',col_list{c},'LineWidth',1.8,'MarkerSize',4,'MarkerFaceColor',col_list{c});
    xline(0,'k--','LineWidth',0.8);
    [~,imax]=max(sig(c,:));  [~,imin]=min(sig(c,:));
    plot(sig(c,imax),z_mm(imax),'r^','MarkerSize',9,'MarkerFaceColor','r','DisplayName','Max tensile');
    plot(sig(c,imin),z_mm(imin),'bv','MarkerSize',9,'MarkerFaceColor','b','DisplayName','Max compressive');
    for zb = res.z*1e3
        yline(zb,'k:','Alpha',0.35,'LineWidth',0.6);
    end
    xlabel(comp_lbl{c},'FontSize',9);
    ylabel('z (mm)','FontSize',9);
    title(comp_lbl{c},'FontSize',10,'FontWeight','bold');
    legend('show','Location','best','FontSize',7);
end
sgtitle('Global Stress Through Laminate Thickness','FontSize',12,'FontWeight','bold');
end

% =========================================================================
% VIZ 10 — GLOBAL STRAIN THROUGH THICKNESS
% =========================================================================
function viz_10_strain_thickness(res)
figure('Name','10 - Strain Through Thickness','NumberTitle','off','Color','w', ...
       'Position',[30 50 780 440]);
z_mm = res.z_pts * 1e3;
eps  = res.eps_g * 1e6;   % microstrain
comp_lbl = {'\epsilon_x (\mu\epsilon)','\epsilon_y (\mu\epsilon)','\gamma_{xy} (\mu\epsilon)'};
col_list  = {[0.1,0.3,0.9],[0.8,0.1,0.1],[0.1,0.6,0.2]};

for c = 1:3
    subplot(1,3,c); hold on; grid on; box on;
    plot(eps(c,:), z_mm, '-o', 'Color',col_list{c},'LineWidth',1.8,'MarkerSize',4,'MarkerFaceColor',col_list{c});
    xline(0,'k--','LineWidth',0.8);
    for zb = res.z*1e3
        yline(zb,'k:','Alpha',0.35,'LineWidth',0.6);
    end
    xlabel(comp_lbl{c},'FontSize',9);
    ylabel('z (mm)','FontSize',9);
    title(comp_lbl{c},'FontSize',10,'FontWeight','bold');
end
sgtitle('Global Strain Through Laminate Thickness','FontSize',12,'FontWeight','bold');
end

% =========================================================================
% VIZ 11 — TSAI-WU FAILURE INDEX BAR CHART
% =========================================================================
function viz_11_TsaiWu_bar(res)
figure('Name','11 - Tsai-Wu Failure Index','NumberTitle','off','Color','w', ...
       'Position',[830 50 580 400]);
n=res.n; FI=res.FI_TW;
cols = arrayfun(@(v) risk_color(v), FI, 'UniformOutput', false);
cols = vertcat(cols{:});
b = bar(1:n, FI, 'FaceColor','flat'); b.CData = cols;
hold on; grid on;
yline(1,'r--','LineWidth',2.0,'Label','Failure  FI=1', ...
      'LabelHorizontalAlignment','right','FontSize',9);
[mx,ic] = max(FI);
plot(ic, mx, 'k*','MarkerSize',13,'LineWidth',2);
text(ic, mx+0.04, sprintf('Critical Ply %d\nFI = %.4f',ic,mx), ...
     'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
xlabel('Ply Number','FontSize',11); ylabel('Tsai-Wu Index','FontSize',11);
title('Tsai-Wu Failure Index per Ply','FontSize',12,'FontWeight','bold');
xticks(1:n);
xticklabels(arrayfun(@(k)sprintf('%d\n(%+g\xB0)',k,res.theta_deg(k)),(1:n)','UniformOutput',false));
ylim([0, max(max(FI)*1.3, 1.3)]);
patch(NaN,NaN,[0.10,0.72,0.33],'DisplayName','Safe');
patch(NaN,NaN,[1.00,0.70,0.00],'DisplayName','Near limit');
patch(NaN,NaN,[0.90,0.15,0.15],'DisplayName','Failed');
legend('Location','northeast','FontSize',9);
end

% =========================================================================
% VIZ 12 — TSAI-HILL FAILURE INDEX BAR CHART
% =========================================================================
function viz_12_TsaiHill_bar(res)
figure('Name','12 - Tsai-Hill Failure Index','NumberTitle','off','Color','w', ...
       'Position',[30 520 580 400]);
n=res.n; FI=res.FI_TH;
cols = arrayfun(@(v) risk_color(v), FI, 'UniformOutput', false);
cols = vertcat(cols{:});
b = bar(1:n, FI, 'FaceColor','flat'); b.CData = cols;
hold on; grid on;
yline(1,'r--','LineWidth',2.0,'Label','Failure  FI=1', ...
      'LabelHorizontalAlignment','right','FontSize',9);
[mx,ic] = max(FI);
plot(ic, mx, 'k*','MarkerSize',13,'LineWidth',2);
text(ic, mx+0.04, sprintf('Critical Ply %d\nFI = %.4f',ic,mx), ...
     'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
xlabel('Ply Number','FontSize',11); ylabel('Tsai-Hill Index','FontSize',11);
title('Tsai-Hill Failure Index per Ply','FontSize',12,'FontWeight','bold');
xticks(1:n);
xticklabels(arrayfun(@(k)sprintf('%d\n(%+g\xB0)',k,res.theta_deg(k)),(1:n)','UniformOutput',false));
ylim([0, max(max(FI)*1.3, 1.3)]);
end

% =========================================================================
% VIZ 13 — MAXIMUM STRESS CRITERION
% =========================================================================
function viz_13_MaxStress_bar(res)
figure('Name','13 - Max Stress Criterion','NumberTitle','off','Color','w', ...
       'Position',[630 520 580 400]);
n=res.n; FI=res.FI_MS;
cols = arrayfun(@(v) risk_color(v), FI, 'UniformOutput', false);
cols = vertcat(cols{:});
b = bar(1:n, FI, 'FaceColor','flat'); b.CData = cols;
hold on; grid on;
yline(1,'r--','LineWidth',2.0,'Label','Failure','LabelHorizontalAlignment','right','FontSize',9);
[mx,ic] = max(FI);
plot(ic, mx, 'k*','MarkerSize',13,'LineWidth',2);
text(ic, mx+0.04, sprintf('Critical Ply %d\nUtil = %.4f',ic,mx), ...
     'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
xlabel('Ply Number','FontSize',11); ylabel('Utilization Ratio','FontSize',11);
title('Maximum Stress Criterion — Utilization per Ply','FontSize',12,'FontWeight','bold');
xticks(1:n);
xticklabels(arrayfun(@(k)sprintf('%d\n(%+g\xB0)',k,res.theta_deg(k)),(1:n)','UniformOutput',false));
ylim([0, max(max(FI)*1.3, 1.3)]);
end

% =========================================================================
% VIZ 14 — DELAMINATION RISK
% =========================================================================
function viz_14_delamination(res)
figure('Name','14 - Delamination Risk','NumberTitle','off','Color','w', ...
       'Position',[30 50 640 400]);
risk  = res.delam_risk;
n_int = numel(risk);
cols  = arrayfun(@(v) risk_color(v), risk, 'UniformOutput', false);
cols  = vertcat(cols{:});
b = bar(1:n_int, risk, 'FaceColor','flat'); b.CData = cols;
hold on; grid on;
yline(0.7,'r:','LineWidth',1.5,'Label','High risk threshold','LabelHorizontalAlignment','right');
yline(0.4,'y:','LineWidth',1.5,'Label','Medium risk threshold','LabelHorizontalAlignment','right');
[~,ic] = max(risk);
plot(ic, risk(ic)+0.04, 'kv','MarkerSize',10,'MarkerFaceColor','k');
text(ic, risk(ic)+0.08, sprintf('Highest Risk\nInterface %d|%d',ic,ic+1), ...
     'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
patch(NaN,NaN,[0.10,0.72,0.33],'DisplayName','Low (\Delta\theta \leq 36\circ)');
patch(NaN,NaN,[1.00,0.70,0.00],'DisplayName','Medium');
patch(NaN,NaN,[0.90,0.15,0.15],'DisplayName','High (\Delta\theta > 63\circ)');
legend('Location','northeast','FontSize',9);
xlabel('Interface (between Ply k and k+1)','FontSize',11);
ylabel('Normalized Delamination Risk','FontSize',11);
title('Delamination Risk at Each Ply Interface','FontSize',12,'FontWeight','bold');
lbl = arrayfun(@(k)sprintf('%d|%d\n%+g\xB0|%+g\xB0',k,k+1,res.theta_deg(k),res.theta_deg(k+1)), ...
               (1:n_int)','UniformOutput',false);
xticks(1:n_int); xticklabels(lbl);
ylim([0, 1.25]);
end

% =========================================================================
% VIZ 15 — BUCKLING SAFETY FACTOR
% =========================================================================
function viz_15_buckling_SF(res)
figure('Name','15 - Buckling Safety Factor','NumberTitle','off','Color','w', ...
       'Position',[690 50 640 400]);
Nx_app = abs(res.NMvec(1));
N_cr   = res.N_cr;
SF     = res.SF_buckle;

subplot(1,2,1); hold on; grid on; box on;
vals = [Nx_app, N_cr];
b = bar(vals,'FaceColor','flat');
if SF >= 1
    b.CData = [[0.90,0.15,0.15]; [0.10,0.72,0.33]];
else
    b.CData = [[0.90,0.15,0.15]; [0.90,0.15,0.15]];
end
set(gca,'XTick',1:2,'XTickLabel',{'N_{x} applied','N_{cr}'},'FontSize',10);
ylabel('Load Intensity (N/m)','FontSize',10);
title('Applied vs Critical Buckling Load','FontSize',10,'FontWeight','bold');
for i=1:2
    text(i, vals(i)*1.03, sprintf('%.0f N/m',vals(i)), ...
         'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end

subplot(1,2,2); hold on; axis equal; axis off;
% Semi-circular gauge
th_g = linspace(0, pi, 200);
r_o = 1.0;  r_i = 0.62;
fill([r_i*cos(th_g), r_o*cos(fliplr(th_g))], ...
     [r_i*sin(th_g), r_o*sin(fliplr(th_g))], [0.88,0.88,0.88]);
% Colored arc up to SF position
SF_frac = min(SF/4, 1);
th_fill = linspace(0, pi*SF_frac, 100);
fill_col = [0.10,0.72,0.33]*(SF>=1) + [0.90,0.15,0.15]*(SF<1);
fill([r_i*cos(th_fill), r_o*cos(fliplr(th_fill))], ...
     [r_i*sin(th_fill), r_o*sin(fliplr(th_fill))], fill_col, 'FaceAlpha',0.85);
% Needle
th_needle = pi * min(SF/4, 1);
plot([0, 0.9*cos(th_needle)],[0, 0.9*sin(th_needle)],'k-','LineWidth',3);
% Labels
text(0,  0.10, sprintf('SF = %.2f', SF),'HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
if SF >= 1
    text(0, -0.25, 'SAFE',   'HorizontalAlignment','center','FontSize',13,'FontWeight','bold','Color',[0,0.55,0]);
else
    text(0, -0.25, 'UNSAFE', 'HorizontalAlignment','center','FontSize',13,'FontWeight','bold','Color',[0.85,0,0]);
end
text(-1.05, 0, '0',   'FontSize',9,'HorizontalAlignment','center');
text( 1.05, 0, '4+',  'FontSize',9,'HorizontalAlignment','center');
text( 0,    1.12,'SF','FontSize',9,'HorizontalAlignment','center');
title(sprintf('Safety Factor Gauge  (mode m=%d)', res.m_cr),'FontSize',10,'FontWeight','bold');
sgtitle('Buckling Analysis','FontSize',13,'FontWeight','bold');
end

% =========================================================================
% VIZ 16 — LAMINATE DEFORMATION
% =========================================================================
function viz_16_deformation(res)
figure('Name','16 - Laminate Deformation','NumberTitle','off','Color','w', ...
       'Position',[30 500 760 420]);
kap = res.kap;
a = res.a_pan;  b = res.b_pan;
[X,Y] = meshgrid(linspace(0,a,35), linspace(0,b,35));
W = 0.5*kap(1)*X.^2 + 0.5*kap(2)*Y.^2 + kap(3)*X.*Y;

scale = 0;
if max(abs(W(:))) > 1e-15
    scale = 0.04 * max(a,b) / max(abs(W(:)));
end
W_plot = W * scale;

subplot(1,2,1);
surf(X*1e3, Y*1e3, zeros(size(X)), 'FaceColor',[0.65,0.85,1], ...
     'EdgeColor','k','FaceAlpha',0.5,'LineWidth',0.4);
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('w');
title('Original (Flat) Shape','FontSize',10,'FontWeight','bold');
grid on; view(30,25);

subplot(1,2,2);
surf(X*1e3, Y*1e3, W_plot, 'FaceColor','interp','EdgeAlpha',0.25,'FaceAlpha',0.90);
colormap(gca, cool); colorbar;
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('w (exaggerated)');
if scale > 0
    title(sprintf('Deformed Shape  (×%.0f exaggeration)', scale),'FontSize',10,'FontWeight','bold');
else
    title('Deformed Shape  (no curvature applied)','FontSize',10,'FontWeight','bold');
end
grid on; view(30,25);
sgtitle('Laminate Bending Deformation under Applied Moments','FontSize',12,'FontWeight','bold');
end

% =========================================================================
% VIZ 17 — BUCKLING MODE SHAPE
% =========================================================================
function viz_17_buckling_mode(res)
figure('Name','17 - Buckling Mode Shape','NumberTitle','off','Color','w', ...
       'Position',[810 500 520 420]);
m = res.m_cr;  n_mode = 1;
a = res.a_pan; b = res.b_pan;
[X,Y] = meshgrid(linspace(0,a,55), linspace(0,b,55));
W = sin(m*pi*X/a) .* sin(n_mode*pi*Y/b);
surf(X*1e3, Y*1e3, W, 'FaceColor','interp','EdgeAlpha',0.12,'FaceAlpha',0.92);
colormap(cool); colorbar;
xlabel('x (mm)','FontSize',11); ylabel('y (mm)','FontSize',11);
zlabel('w (normalized)','FontSize',11);
title(sprintf('First Buckling Mode Shape  (m=%d, n=%d)\nN_{cr} = %.0f N/m', ...
      m, n_mode, res.N_cr),'FontSize',12,'FontWeight','bold');
grid on; view(42,28);
end

% =========================================================================
% VIZ 18 — POLAR STIFFNESS PLOT
% =========================================================================
function viz_18_polar_stiffness(res)
figure('Name','18 - Polar Stiffness','NumberTitle','off','Color','w', ...
       'Position',[30 50 520 480]);
a_comp = inv(res.A);
a11 = a_comp(1,1); a12 = a_comp(1,2); a22 = a_comp(2,2);
a66 = a_comp(3,3); a16 = a_comp(1,3); a26 = a_comp(2,3);
h = res.h;

ang = deg2rad(0:1:360);
Ex_theta = zeros(size(ang));
for i = 1:numel(ang)
    c = cos(ang(i));  s = sin(ang(i));
    a11_rot = a11*c^4 + (2*a12+a66)*c^2*s^2 + a22*s^4 + 2*a16*c^3*s + 2*a26*c*s^3;
    Ex_theta(i) = 1 / (h * a11_rot + eps);
end
Ex_GPa = Ex_theta / 1e9;

ax = polaraxes;
polarplot(ax, ang, Ex_GPa, 'b-', 'LineWidth', 2.5);
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaTick = 0:30:330;
title('Laminate Stiffness E(\theta)  [GPa]','FontSize',12,'FontWeight','bold');
% Annotate max/min
[Emax, im] = max(Ex_GPa);  [Emin, imn] = min(Ex_GPa);
hold(ax,'on');
polarplot(ax, ang(im),  Emax, 'r^','MarkerSize',10,'MarkerFaceColor','r','DisplayName',sprintf('E_{max}=%.1f GPa',Emax));
polarplot(ax, ang(imn), Emin, 'bv','MarkerSize',10,'MarkerFaceColor','b','DisplayName',sprintf('E_{min}=%.1f GPa',Emin));
legend(ax,'Location','southoutside','FontSize',9);
end

% =========================================================================
% VIZ 19 — FAILURE ENVELOPE (Tsai-Wu, sigma1 vs sigma2, tau12=0)
% =========================================================================
function viz_19_failure_envelope(res)
figure('Name','19 - Failure Envelope','NumberTitle','off','Color','w', ...
       'Position',[570 50 540 460]);
mat = res.mat;
F1  = 1/mat.Xt - 1/mat.Xc;
F2  = 1/mat.Yt - 1/mat.Yc;
F11 = 1/(mat.Xt*mat.Xc);
F22 = 1/(mat.Yt*mat.Yc);
F12 = -0.5*sqrt(F11*F22);

% Solve quadratic F22*s2^2 + (F2+2F12*s1)*s2 + (F11*s1^2+F1*s1-1) = 0
s1_vec = linspace(-mat.Xc, mat.Xt, 600) / 1e6;   % MPa
s2_hi  = nan(size(s1_vec));
s2_lo  = nan(size(s1_vec));
for i = 1:numel(s1_vec)
    s1 = s1_vec(i)*1e6;
    A_q = F22;
    B_q = F2 + 2*F12*s1;
    C_q = F11*s1^2 + F1*s1 - 1;
    disc = B_q^2 - 4*A_q*C_q;
    if disc >= 0
        s2_hi(i) = (-B_q + sqrt(disc))/(2*A_q) / 1e6;
        s2_lo(i) = (-B_q - sqrt(disc))/(2*A_q) / 1e6;
    end
end
ok = ~isnan(s2_hi);
hold on; grid on;
fill([s1_vec(ok), fliplr(s1_vec(ok))], [s2_hi(ok), fliplr(s2_lo(ok))], ...
     [0.80,1.00,0.80],'EdgeColor',[0,0.6,0],'LineWidth',2.0, ...
     'FaceAlpha',0.35,'DisplayName','Safe region');
plot(s1_vec(ok), s2_hi(ok), 'g-','LineWidth',2,'DisplayName','Failure boundary');
plot(s1_vec(ok), s2_lo(ok), 'g-','LineWidth',2,'HandleVisibility','off');

% Current operating point (critical ply)
[~,ic] = max(res.FI_TW);
zm = (res.z(ic)+res.z(ic+1))/2;
ep = res.eps0 + zm*res.kap;
sg = res.Qbar(:,:,ic)*ep;
c = cosd(res.theta_deg(ic)); s_k = sind(res.theta_deg(ic));
T = [c^2,s_k^2,2*c*s_k; s_k^2,c^2,-2*c*s_k; -c*s_k,c*s_k,c^2-s_k^2];
sl = T*sg;
plot(sl(1)/1e6, sl(2)/1e6, 'r*','MarkerSize',14,'LineWidth',2, ...
     'DisplayName',sprintf('Operating point (Ply %d)',ic));

xline(0,'k:','LineWidth',0.8); yline(0,'k:','LineWidth',0.8);
xlabel('\sigma_1 (MPa)','FontSize',11);
ylabel('\sigma_2 (MPa)','FontSize',11);
title('Tsai-Wu Failure Envelope  (\tau_{12}=0 plane)','FontSize',12,'FontWeight','bold');
legend('Location','northeast','FontSize',9);
end

% =========================================================================
% VIZ 20 — PLY CONTRIBUTION ANALYSIS
% =========================================================================
function viz_20_ply_contribution(res)
figure('Name','20 - Ply Contribution','NumberTitle','off','Color','w', ...
       'Position',[30 50 800 540]);
n    = res.n;
data = [res.ply_Ex_pct; res.ply_Ey_pct; res.ply_Gxy_pct; res.ply_D11_pct];
ttls = {'Contribution to E_x (%)','Contribution to E_y (%)', ...
        'Contribution to G_{xy} (%)','Contribution to D_{11} (%)'};
cmap_p = parula(n);
for p = 1:4
    subplot(2,2,p);
    bh = bar(1:n, data(p,:),'FaceColor','flat');
    for k=1:n; bh.CData(k,:) = cmap_p(k,:); end
    hold on; grid on; box on;
    [~,imax] = max(data(p,:));
    text(imax, data(p,imax)+1.5, sprintf('%.1f%%',data(p,imax)), ...
         'HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Color','r');
    xlabel('Ply','FontSize',9); ylabel('%','FontSize',9);
    title(ttls{p},'FontSize',10,'FontWeight','bold');
    xticks(1:n);
    xticklabels(arrayfun(@(k)sprintf('%d\n(%+g\xB0)',k,res.theta_deg(k)),(1:n)','UniformOutput',false));
end
sgtitle('Ply Contribution to Laminate Stiffness','FontSize',13,'FontWeight','bold');
end

% =========================================================================
% VIZ 21 — ENGINEERING DASHBOARD (summary figure)
% =========================================================================
function viz_21_dashboard(res)
figure('Name','21 - Engineering Dashboard','NumberTitle','off','Color','w', ...
       'Position',[50 50 1300 820]);

[mTW, iTW] = max(res.FI_TW);
[~,  iDR ] = max(res.delam_risk);

% ---- (1) Laminate stack mini ----
subplot(3,4,1); hold on; axis off;
n=res.n; h=res.h;
for k=1:n
    yb=(k-1)/n; yt=k/n;
    fill([0,1,1,0],[yb,yb,yt,yt],ang_color(res.theta_deg(k),res.theta_deg), ...
         'EdgeColor','k','LineWidth',0.7,'FaceAlpha',0.8);
    text(0.5,(yb+yt)/2,sprintf('%+g\xB0',res.theta_deg(k)),...
         'HorizontalAlignment','center','FontSize',7,'FontWeight','bold');
end
text(0.5,1.06,'TOP','HorizontalAlignment','center','FontSize',7,'Color',[0.5,0,0]);
text(0.5,-0.06,'BTM','HorizontalAlignment','center','FontSize',7,'Color',[0.5,0,0]);
title('Layup','FontSize',9,'FontWeight','bold');

% ---- (2) Engineering constants ----
subplot(3,4,2); axis off;
text(0.5,0.96,'ENGINEERING CONSTANTS','HorizontalAlignment','center','FontWeight','bold','FontSize',9,'Units','normalized');
vals2 = {sprintf('E_x   = %.2f GPa',res.Ex/1e9); ...
         sprintf('E_y   = %.2f GPa',res.Ey/1e9); ...
         sprintf('G_{xy} = %.2f GPa',res.Gxy/1e9); ...
         sprintf('\\nu_{xy} = %.4f',res.nuxy); ...
         sprintf('h     = %.3f mm',res.h*1e3); ...
         sprintf('\\rho A  = %.3f kg/m^2',res.rho_area)};
for vi=1:numel(vals2)
    text(0.05,0.83-(vi-1)*0.13,vals2{vi},'FontSize',9,'Units','normalized');
end

% ---- (3) Failure summary ----
subplot(3,4,3); axis off;
text(0.5,0.96,'FAILURE SUMMARY','HorizontalAlignment','center','FontWeight','bold','FontSize',9,'Units','normalized');
col_tw = risk_color(mTW);
text(0.05,0.82,sprintf('Max Tsai-Wu:   %.4f',mTW),'FontSize',9,'Color',col_tw,'FontWeight','bold','Units','normalized');
text(0.05,0.70,sprintf('Critical Ply:  %d (%+g\xB0)',iTW,res.theta_deg(iTW)),'FontSize',9,'Units','normalized');
text(0.05,0.58,sprintf('Max Tsai-Hill: %.4f',max(res.FI_TH)),'FontSize',9,'Units','normalized');
text(0.05,0.46,sprintf('Max MS util:   %.4f',max(res.FI_MS)),'FontSize',9,'Units','normalized');
if mTW < 1
    text(0.5,0.28,'NO FAILURE','HorizontalAlignment','center','FontSize',12,'FontWeight','bold','Color',[0,0.55,0],'Units','normalized');
else
    text(0.5,0.28,'FAILURE','HorizontalAlignment','center','FontSize',12,'FontWeight','bold','Color',[0.85,0,0],'Units','normalized');
end

% ---- (4) Buckling summary ----
subplot(3,4,4); axis off;
text(0.5,0.96,'BUCKLING','HorizontalAlignment','center','FontWeight','bold','FontSize',9,'Units','normalized');
col_sf = risk_color(1 - min(res.SF_buckle/3, 1));
text(0.05,0.82,sprintf('N_{cr} = %.0f N/m',res.N_cr),'FontSize',9,'Units','normalized');
text(0.05,0.70,sprintf('N_x   = %.0f N/m',abs(res.NMvec(1))),'FontSize',9,'Units','normalized');
text(0.05,0.58,sprintf('SF    = %.2f',res.SF_buckle),'FontSize',9,'Color',col_sf,'FontWeight','bold','Units','normalized');
text(0.05,0.46,sprintf('Mode  m = %d half-waves',res.m_cr),'FontSize',9,'Units','normalized');
[~,idr]=max(res.delam_risk);
text(0.05,0.30,sprintf('Max delam. risk: Interface %d|%d',idr,idr+1),'FontSize',9,'Units','normalized');

% ---- (5) ABD heatmap mini ----
subplot(3,4,5);
imagesc(log10(abs(res.ABD)+1)); colormap(gca,parula); axis off; axis square;
hold on;
line([3.5,3.5],[0.5,6.5],'Color','r','LineWidth',1.5);
line([0.5,6.5],[3.5,3.5],'Color','r','LineWidth',1.5);
title('[ABD] (log scale)','FontSize',9,'FontWeight','bold');

% ---- (6) Tsai-Wu bar mini ----
subplot(3,4,6);
bar(1:res.n, res.FI_TW,'FaceColor',[0.2,0.5,0.9]);
hold on; grid on;
yline(1,'r--','LineWidth',1.5);
xlabel('Ply','FontSize',8); ylabel('FI','FontSize',8);
title('Tsai-Wu Index','FontSize',9,'FontWeight','bold');
ylim([0, max(max(res.FI_TW)*1.3,1.3)]);

% ---- (7) Stress through thickness mini ----
subplot(3,4,7); hold on; grid on;
plot(res.sig_g(1,:)/1e6, res.z_pts*1e3,'b-o','MarkerSize',3,'LineWidth',1.3);
xline(0,'k:');
for zb=res.z*1e3; yline(zb,'k:','Alpha',0.3,'LineWidth',0.5); end
xlabel('\sigma_x (MPa)','FontSize',8); ylabel('z (mm)','FontSize',8);
title('\sigma_x Through Thickness','FontSize',9,'FontWeight','bold');

% ---- (8) Polar stiffness mini ----
subplot(3,4,8);
a_comp = inv(res.A);
a11=a_comp(1,1); a12=a_comp(1,2); a22=a_comp(2,2);
a66=a_comp(3,3); a16=a_comp(1,3); a26=a_comp(2,3);
ang_p = deg2rad(0:2:360);
Ex_p  = zeros(size(ang_p));
for i=1:numel(ang_p)
    c=cos(ang_p(i)); s=sin(ang_p(i));
    a11r = a11*c^4+(2*a12+a66)*c^2*s^2+a22*s^4+2*a16*c^3*s+2*a26*c*s^3;
    Ex_p(i) = 1/(res.h*a11r+eps);
end
polarplot(ang_p, Ex_p/1e9,'b-','LineWidth',1.8);
title('E(\theta) [GPa]','FontSize',9,'FontWeight','bold');

% ---- (9) Buckling bar mini ----
subplot(3,4,9);
bar_d = [abs(res.NMvec(1)), res.N_cr];
bh = bar(bar_d,'FaceColor','flat');
bh.CData = [[0.90,0.15,0.15];[0.10,0.72,0.33]];
set(gca,'XTick',1:2,'XTickLabel',{'N_{app}','N_{cr}'},'FontSize',8);
ylabel('N (N/m)','FontSize',8); grid on;
title(sprintf('Buckling  SF=%.2f',res.SF_buckle),'FontSize',9,'FontWeight','bold');

% ---- (10) Delamination mini ----
subplot(3,4,10);
risk = res.delam_risk;
cols_d = arrayfun(@(v) risk_color(v), risk, 'UniformOutput', false);
cols_d = vertcat(cols_d{:});
bh2 = bar(1:numel(risk), risk,'FaceColor','flat'); bh2.CData=cols_d;
grid on; ylim([0,1.25]);
xlabel('Interface','FontSize',8); ylabel('Risk','FontSize',8);
title('Delamination Risk','FontSize',9,'FontWeight','bold');

% ---- (11) Ply contribution D11 mini ----
subplot(3,4,11);
bar(1:res.n, res.ply_D11_pct,'FaceColor',[0.3,0.7,0.5]);
grid on;
xlabel('Ply','FontSize',8); ylabel('%','FontSize',8);
title('D_{11} Ply Contribution (%)','FontSize',9,'FontWeight','bold');

% ---- (12) Text summary box ----
subplot(3,4,12); axis off;
text(0.5,1.00,'LAMINATE SUMMARY','HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Units','normalized');
sumlines = {sprintf('Plies: %d',res.n); ...
            sprintf('Layup: [%s]',num2str(res.theta_deg)); ...
            sprintf('h = %.3f mm',res.h*1e3); ...
            sprintf('\\rho A = %.3f kg/m^2',res.rho_area); ...
            ' '; ...
            sprintf('E_x = %.2f GPa',res.Ex/1e9); ...
            sprintf('E_y = %.2f GPa',res.Ey/1e9); ...
            sprintf('G_{xy} = %.2f GPa',res.Gxy/1e9); ...
            sprintf('\\nu_{xy} = %.4f',res.nuxy); ...
            ' '; ...
            sprintf('Max \\sigma_x = %.1f MPa',max(abs(res.sig_g(1,:)))/1e6); ...
            sprintf('Max TW index = %.4f',mTW); ...
            sprintf('Critical Ply = %d',iTW); ...
            sprintf('N_{cr} = %.0f N/m',res.N_cr); ...
            sprintf('SF = %.2f',res.SF_buckle); ...
            sprintf('Max delam. iface = %d|%d',iDR,iDR+1)};
for vi=1:numel(sumlines)
    text(0.03, 0.92-(vi-1)*0.055, sumlines{vi},'FontSize',7.5,'Units','normalized');
end

sgtitle('COMPOSITE LAMINATE — ENGINEERING DASHBOARD', ...
        'FontSize',14,'FontWeight','bold','Color',[0.08,0.08,0.45]);
end

% =========================================================================
% HELPER: blue-white-red diverging colormap
% =========================================================================
function cmap = coolwarm_cmap()
n = 64;
r = [linspace(0.17,1,n/2), linspace(1,0.70,n/2)]';
g = [linspace(0.51,1,n/2), linspace(1,0.09,n/2)]';
b = [linspace(0.72,1,n/2), linspace(1,0.09,n/2)]';
cmap = [r, g, b];
end
