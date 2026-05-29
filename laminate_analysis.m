function res = laminate_analysis()
% LAMINATE_ANALYSIS  Interactive Phase 1 Classical Laminate Theory analysis.
% Run from the MATLAB command window:  res = laminate_analysis();

clc;
fprintf('=======================================================\n');
fprintf('       COMPOSITE LAMINATE ANALYSIS  -  PHASE 1\n');
fprintf('=======================================================\n\n');

% -----------------------------------------------------------------------
%  SECTION A: MATERIAL PROPERTIES
% -----------------------------------------------------------------------
fprintf('--- MATERIAL PROPERTIES ---\n');
fprintf('Enter elastic moduli and strengths for a single ply.\n\n');

mat.E1 = prompt_scalar( ...
    'Longitudinal modulus E1 [GPa]  (e.g. 181): ', 1e9);

mat.E2 = prompt_scalar( ...
    'Transverse modulus   E2 [GPa]  (e.g. 10.3): ', 1e9);

mat.G12 = prompt_scalar( ...
    'Shear modulus       G12 [GPa]  (e.g. 7.17): ', 1e9);

mat.nu12 = prompt_scalar( ...
    'Major Poisson ratio nu12 [ ]   (e.g. 0.28): ', 1);

mat.t = prompt_scalar( ...
    'Ply thickness         t  [mm]  (e.g. 0.125): ', 1e-3);

fprintf('\n Strengths:\n');
mat.Xt = prompt_scalar( ...
    '  Longitudinal tensile  Xt [MPa]  (e.g. 1500): ', 1e6);

mat.Xc = prompt_scalar( ...
    '  Longitudinal compress Xc [MPa]  (e.g. 1500): ', 1e6);

mat.Yt = prompt_scalar( ...
    '  Transverse tensile    Yt [MPa]  (e.g.   40): ', 1e6);

mat.Yc = prompt_scalar( ...
    '  Transverse compress   Yc [MPa]  (e.g.  246): ', 1e6);

mat.S = prompt_scalar( ...
    '  In-plane shear        S  [MPa]  (e.g.   68): ', 1e6);

mat.rho = prompt_scalar( ...
    'Ply density          rho [kg/m^3] (e.g. 1600): ', 1);

mat.nu21 = mat.nu12 * mat.E2 / mat.E1;

% -----------------------------------------------------------------------
%  SECTION B: LAYUP DEFINITION
% -----------------------------------------------------------------------
fprintf('\n--- LAYUP DEFINITION ---\n');
fprintf('Enter ply angles in degrees, from PLY 1 (bottom) to PLY N (top).\n');
fprintf('Format : space-separated list inside square brackets.\n');
fprintf('Example : [0 45 -45 90 90 -45 45 0]\n\n');

theta_deg = prompt_vector('Ply angles [deg]: ');
n  = numel(theta_deg);
th = deg2rad(theta_deg);
h  = n * mat.t;
z  = linspace(-h/2, h/2, n+1);   % ply boundary z-coordinates [m]

fprintf('  -> %d plies detected.  Total thickness = %.3f mm\n', n, h*1e3);

% -----------------------------------------------------------------------
%  SECTION C: APPLIED LOADS
% -----------------------------------------------------------------------
fprintf('\n--- APPLIED LOADS ---\n');
fprintf('In-plane stress resultants [N/m] and moment resultants [N.m/m].\n\n');

Nx  = prompt_scalar('Nx   [N/m]     (e.g. 1000, tension positive): ', 1);
Ny  = prompt_scalar('Ny   [N/m]     (e.g.    0): ', 1);
Nxy = prompt_scalar('Nxy  [N/m]     (e.g.    0): ', 1);
Mx  = prompt_scalar('Mx   [N.m/m]   (e.g.   10, sagging positive): ', 1);
My  = prompt_scalar('My   [N.m/m]   (e.g.    0): ', 1);
Mxy = prompt_scalar('Mxy  [N.m/m]   (e.g.    0): ', 1);
NMvec = [Nx; Ny; Nxy; Mx; My; Mxy];

% -----------------------------------------------------------------------
%  SECTION D: PANEL DIMENSIONS  (for buckling)
% -----------------------------------------------------------------------
fprintf('\n--- PANEL DIMENSIONS (for buckling) ---\n');
fprintf('Simply-supported rectangular plate.\n\n');
a_pan = prompt_scalar('Panel length  a [mm]  (x-direction, e.g. 500): ', 1e-3);
b_pan = prompt_scalar('Panel width   b [mm]  (y-direction, e.g. 500): ', 1e-3);

% -----------------------------------------------------------------------
%  COMPUTATION
% -----------------------------------------------------------------------
fprintf('\n=======================================================\n');
fprintf('  Computing...\n');
fprintf('=======================================================\n\n');

%% 1. REDUCED STIFFNESS MATRIX  [Q]
d = 1 - mat.nu12 * mat.nu21;
Q = [mat.E1/d,          mat.nu12*mat.E2/d, 0;
     mat.nu12*mat.E2/d, mat.E2/d,          0;
     0,                 0,                 mat.G12];

%% 2. TRANSFORMED STIFFNESS  Qbar(k)
Q11=Q(1,1); Q12=Q(1,2); Q22=Q(2,2); Q66=Q(3,3);
Qbar = zeros(3,3,n);
for k = 1:n
    c=cos(th(k)); s=sin(th(k));
    c2=c^2; s2=s^2; c4=c^4; s4=s^4; c2s2=c2*s2;
    Qbar(1,1,k) = Q11*c4 + 2*(Q12+2*Q66)*c2s2 + Q22*s4;
    Qbar(1,2,k) = (Q11+Q22-4*Q66)*c2s2 + Q12*(c4+s4);
    Qbar(2,2,k) = Q11*s4 + 2*(Q12+2*Q66)*c2s2 + Q22*c4;
    Qbar(1,3,k) = (Q11-Q12-2*Q66)*c^3*s  - (Q22-Q12-2*Q66)*c*s^3;
    Qbar(2,3,k) = (Q11-Q12-2*Q66)*c*s^3  - (Q22-Q12-2*Q66)*c^3*s;
    Qbar(3,3,k) = (Q11+Q22-2*Q12-2*Q66)*c2s2 + Q66*(c4+s4);
    Qbar(2,1,k)=Qbar(1,2,k); Qbar(3,1,k)=Qbar(1,3,k); Qbar(3,2,k)=Qbar(2,3,k);
end

%% 3. ABD MATRIX
A=zeros(3); B=zeros(3); D=zeros(3);
for k = 1:n
    dz  = z(k+1)  - z(k);
    dz2 = z(k+1)^2 - z(k)^2;
    dz3 = z(k+1)^3 - z(k)^3;
    A = A + Qbar(:,:,k) * dz;
    B = B + Qbar(:,:,k) * dz2/2;
    D = D + Qbar(:,:,k) * dz3/3;
end
ABD = [A, B; B, D];

%% 4. LAMINATE ENGINEERING CONSTANTS
a_comp = inv(A);
Ex   = 1 / (h * a_comp(1,1));
Ey   = 1 / (h * a_comp(2,2));
Gxy  = 1 / (h * a_comp(3,3));
nuxy = -a_comp(1,2) / a_comp(1,1);

%% 5. MIDPLANE STRAINS AND CURVATURES
eps_kap = ABD \ NMvec;
eps0 = eps_kap(1:3);
kap  = eps_kap(4:6);

%% 6. STRESS AND STRAIN THROUGH THICKNESS
npts = 2*n;
z_s  = zeros(1,npts); eps_g=zeros(3,npts); sig_g=zeros(3,npts); sig_l=zeros(3,npts);
idx  = 1;
for k = 1:n
    for zi = [z(k), z(k+1)]
        ep = eps0 + zi*kap;
        sg = Qbar(:,:,k)*ep;
        c=cos(th(k)); s=sin(th(k));
        T=[c^2,s^2,2*c*s; s^2,c^2,-2*c*s; -c*s,c*s,c^2-s^2];
        z_s(idx)=zi; eps_g(:,idx)=ep; sig_g(:,idx)=sg; sig_l(:,idx)=T*sg;
        idx=idx+1;
    end
end
[z_u,ia] = unique(z_s,'stable');
eps_g_u=eps_g(:,ia); sig_g_u=sig_g(:,ia); sig_l_u=sig_l(:,ia);

%% 7. FAILURE INDICES  (ply midplane)
FI_TW=zeros(1,n); FI_TH=zeros(1,n); FI_MS=zeros(1,n);
F1=1/mat.Xt-1/mat.Xc; F2=1/mat.Yt-1/mat.Yc;
F11=1/(mat.Xt*mat.Xc); F22=1/(mat.Yt*mat.Yc);
F66=1/mat.S^2; F12=-0.5*sqrt(F11*F22);
for k=1:n
    zm=(z(k)+z(k+1))/2;
    ep=eps0+zm*kap; sg=Qbar(:,:,k)*ep;
    c=cos(th(k)); s=sin(th(k));
    T=[c^2,s^2,2*c*s; s^2,c^2,-2*c*s; -c*s,c*s,c^2-s^2];
    sl=T*sg; s1=sl(1); s2=sl(2); t12=sl(3);
    FI_TW(k)=F1*s1+F2*s2+F11*s1^2+F22*s2^2+F66*t12^2+2*F12*s1*s2;
    X=mat.Xt*(s1>=0)+mat.Xc*(s1<0); Y=mat.Yt*(s2>=0)+mat.Yc*(s2<0);
    FI_TH(k)=(s1/X)^2-s1*s2/X^2+(s2/Y)^2+(t12/mat.S)^2;
    FI_MS(k)=max([abs(s1)/(mat.Xt*(s1>=0)+mat.Xc*(s1<0)), ...
                  abs(s2)/(mat.Yt*(s2>=0)+mat.Yc*(s2<0)), abs(t12)/mat.S]);
end

%% 8. DELAMINATION RISK
delam_risk=zeros(1,n-1);
for k=1:n-1
    dt=abs(theta_deg(k+1)-theta_deg(k));
    delam_risk(k)=min(dt,180-dt)/90;
end

%% 9. BUCKLING
D11=D(1,1); D12=D(1,2); D22=D(2,2); D66=D(3,3);
N_cr_all=zeros(1,10);
for m=1:10
    R=m*b_pan/a_pan;
    N_cr_all(m)=(pi/b_pan)^2*(D11*R^2+2*(D12+2*D66)+D22/R^2);
end
N_cr=min(N_cr_all); m_cr=find(N_cr_all==N_cr,1);
SF_buckle=N_cr/max(abs(Nx),1e-10);

%% 10. PLY CONTRIBUTIONS
ply_Ex_pct=zeros(1,n); ply_Ey_pct=zeros(1,n);
ply_Gxy_pct=zeros(1,n); ply_D11_pct=zeros(1,n);
for k=1:n
    dz=z(k+1)-z(k); dz3=(z(k+1)^3-z(k)^3)/3;
    ply_Ex_pct(k)=Qbar(1,1,k)*dz; ply_Ey_pct(k)=Qbar(2,2,k)*dz;
    ply_Gxy_pct(k)=Qbar(3,3,k)*dz; ply_D11_pct(k)=Qbar(1,1,k)*dz3;
end
pct=@(v)100*v/(sum(abs(v))+eps);
ply_Ex_pct=pct(ply_Ex_pct); ply_Ey_pct=pct(ply_Ey_pct);
ply_Gxy_pct=pct(ply_Gxy_pct); ply_D11_pct=pct(ply_D11_pct);

% -----------------------------------------------------------------------
%  RESULTS OUTPUT
% -----------------------------------------------------------------------
is_sym = max(abs(theta_deg - flip(theta_deg))) < 0.1;
[mTW, iTW] = max(FI_TW);
[~,  iDR ] = max(delam_risk);

fprintf('=======================================================\n');
fprintf('  RESULTS\n');
fprintf('=======================================================\n\n');

fprintf('LAMINATE GEOMETRY\n');
fprintf('  Layup          : [%s]\n', num2str(theta_deg));
fprintf('  No. of plies   : %d\n', n);
fprintf('  Symmetry       : %s\n', ternary(is_sym,'Symmetric','Unsymmetric'));
fprintf('  Total thickness: %.4f mm\n', h*1e3);
fprintf('  Weight/area    : %.4f kg/m^2\n\n', mat.rho*h);

fprintf('ENGINEERING CONSTANTS (in-plane)\n');
fprintf('  Ex   = %10.3f GPa\n', Ex/1e9);
fprintf('  Ey   = %10.3f GPa\n', Ey/1e9);
fprintf('  Gxy  = %10.3f GPa\n', Gxy/1e9);
fprintf('  nuxy = %10.4f\n\n',   nuxy);

fprintf('MIDPLANE STRAINS & CURVATURES\n');
fprintf('  eps0x  = %10.4e\n', eps0(1));
fprintf('  eps0y  = %10.4e\n', eps0(2));
fprintf('  gam0xy = %10.4e\n', eps0(3));
fprintf('  kx     = %10.4e  1/m\n', kap(1));
fprintf('  ky     = %10.4e  1/m\n', kap(2));
fprintf('  kxy    = %10.4e  1/m\n\n', kap(3));

fprintf('Q MATRIX  [GPa]\n');
print_matrix(Q/1e9);

fprintf('A MATRIX  [MN/m]\n');
print_matrix(A/1e6);

fprintf('B MATRIX  [N]\n');
print_matrix(B);
if max(abs(B(:))) < 1e-8
    fprintf('  -> B ≈ 0  (symmetric laminate, no extension-bending coupling)\n');
end
fprintf('\n');

fprintf('D MATRIX  [N.m]\n');
print_matrix(D);

fprintf('FAILURE INDICES  (evaluated at ply midplane)\n');
fprintf('  %-6s  %-8s  %-12s  %-12s  %-12s  %-12s\n', ...
        'Ply','Angle','Tsai-Wu','Tsai-Hill','Max Stress','Status');
fprintf('  %s\n', repmat('-',1,66));
for k=1:n
    if FI_TW(k)>=1
        status='FAILED';
    elseif FI_TW(k)>=0.7
        status='Near limit';
    else
        status='OK';
    end
    fprintf('  %-6d  %-8.1f  %-12.4f  %-12.4f  %-12.4f  %s\n', ...
            k, theta_deg(k), FI_TW(k), FI_TH(k), FI_MS(k), status);
end
fprintf('  -> Critical ply (max Tsai-Wu): Ply %d (%+g deg)  FI = %.4f\n\n', ...
        iTW, theta_deg(iTW), mTW);

fprintf('DELAMINATION RISK\n');
fprintf('  %-12s  %-14s  %-12s  %s\n','Interface','Delta-theta','Risk','Level');
fprintf('  %s\n', repmat('-',1,52));
for k=1:n-1
    dt=abs(theta_deg(k+1)-theta_deg(k));
    dt=min(dt,180-dt);
    if delam_risk(k)>0.7; lvl='HIGH';
    elseif delam_risk(k)>0.4; lvl='MEDIUM';
    else; lvl='LOW'; end
    fprintf('  %d | %d         %-14.1f  %-12.3f  %s\n', ...
            k,k+1,dt,delam_risk(k),lvl);
end
fprintf('  -> Highest risk: Interface %d|%d\n\n', iDR, iDR+1);

fprintf('BUCKLING  (simply-supported, uniaxial Nx)\n');
fprintf('  Critical load N_cr = %.2f N/m  (m=%d half-waves)\n', N_cr, m_cr);
fprintf('  Applied    Nx      = %.2f N/m\n', Nx);
fprintf('  Safety factor SF   = %.3f', SF_buckle);
if SF_buckle >= 1
    fprintf('  -> SAFE\n\n');
else
    fprintf('  -> UNSAFE (buckling predicted)\n\n');
end

fprintf('PLY CONTRIBUTIONS TO STIFFNESS  [%%]\n');
fprintf('  %-6s  %-8s  %-10s  %-10s  %-10s  %-10s\n', ...
        'Ply','Angle','Ex','Ey','Gxy','D11');
fprintf('  %s\n', repmat('-',1,60));
for k=1:n
    fprintf('  %-6d  %-8.1f  %-10.2f  %-10.2f  %-10.2f  %-10.2f\n', ...
            k, theta_deg(k), ply_Ex_pct(k), ply_Ey_pct(k), ply_Gxy_pct(k), ply_D11_pct(k));
end
fprintf('\n');

fprintf('GLOBAL STRESSES THROUGH THICKNESS  [MPa]\n');
fprintf('  %-8s  %-12s  %-12s  %-12s  %-6s\n','z [mm]','sigma_x','sigma_y','tau_xy','Ply');
fprintf('  %s\n', repmat('-',1,56));
for i=1:numel(z_u)
    ply_k=max(1,min(n,sum(z_u(i)>=z)-1));  % which ply this z belongs to
    if z_u(i)==z(end); ply_k=n; end
    fprintf('  %-8.3f  %-12.3f  %-12.3f  %-12.3f  %d\n', ...
            z_u(i)*1e3, sig_g_u(1,i)/1e6, sig_g_u(2,i)/1e6, sig_g_u(3,i)/1e6, ply_k);
end
fprintf('\n');

fprintf('=======================================================\n');
fprintf('  ANALYSIS COMPLETE\n');
fprintf('=======================================================\n');

% -----------------------------------------------------------------------
%  PACK RESULTS  (return struct for optional post-processing)
% -----------------------------------------------------------------------
res.n=n; res.theta_deg=theta_deg; res.mat=mat; res.h=h; res.z=z;
res.Q=Q; res.Qbar=Qbar; res.A=A; res.B=B; res.D=D; res.ABD=ABD;
res.Ex=Ex; res.Ey=Ey; res.Gxy=Gxy; res.nuxy=nuxy;
res.eps0=eps0; res.kap=kap;
res.z_pts=z_u; res.sig_g=sig_g_u; res.eps_g=eps_g_u; res.sig_l=sig_l_u;
res.FI_TW=FI_TW; res.FI_TH=FI_TH; res.FI_MS=FI_MS;
res.delam_risk=delam_risk; res.N_cr=N_cr; res.SF_buckle=SF_buckle; res.m_cr=m_cr;
res.NMvec=NMvec; res.a_pan=a_pan; res.b_pan=b_pan;
res.ply_Ex_pct=ply_Ex_pct; res.ply_Ey_pct=ply_Ey_pct;
res.ply_Gxy_pct=ply_Gxy_pct; res.ply_D11_pct=ply_D11_pct;
res.rho_area=mat.rho*h;
end

% =========================================================================
%  INPUT HELPERS
% =========================================================================
function val = prompt_scalar(msg, scale)
% Reads a single number and multiplies by scale to convert to SI units.
while true
    raw = input(['  ' msg]);
    if isnumeric(raw) && isscalar(raw)
        val = raw * scale;
        return;
    end
    fprintf('  Invalid input. Enter a single number.\n');
end
end

function vec = prompt_vector(msg)
% Reads a row vector entered as [a b c ...] or a b c ...
while true
    raw = input(['  ' msg]);
    if isnumeric(raw) && isvector(raw) && numel(raw) >= 1
        vec = raw(:)';
        return;
    end
    fprintf('  Invalid input. Enter values inside square brackets, e.g. [0 45 -45 90].\n');
end
end

% =========================================================================
%  OUTPUT HELPERS
% =========================================================================
function print_matrix(M)
for i=1:size(M,1)
    fprintf('  [');
    fprintf('  %10.4f', M(i,:));
    fprintf('  ]\n');
end
fprintf('\n');
end

function s = ternary(cond, s_true, s_false)
if cond; s=s_true; else; s=s_false; end
end
