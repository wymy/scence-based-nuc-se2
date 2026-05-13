function sbnuc_mf_se2_lie_ablation_demo()
% SBNUC_MF_SE2_LIE_ABLATION_DEMO (fixed: non-empty Fig2 + uitable compatible)
% ------------------------------------------------------------
% 输出：
%   图1 Success map（Δt-Δθ 网格热力图，3模式并排）
%   图2 收敛曲线（自动选择一个“可成功且差异明显”的代表点：Energy & NCC vs Iter）
%   表1 全局统计表（SuccessRate/MeanIters/MeanTime/MeanEnergy/MeanNCC）
%
% 关键修复：
%   1) Fig2 为空：E_hist/N_hist 可能为空 -> 强制写入初值/末值，且代表点自动选择 + 二次复跑
%   2) uitable 报错/显示空：table 的 string 列不兼容 -> 转成 cellstr(char) 给 uitable
% ------------------------------------------------------------

%% ========== 用户参数区 ==========
bin_dir     = 'E:\scene based\第二次采\外场晃动3000帧\打快门30分钟后\';
target_name = '';   % 留空自动取最新 .bin
W = 1024; H = 1280;

start_frame = 1561;
window_len  = 30;

% ======= Ablation 扰动网格（更“好看”的中等扰动） =======
grid_rot_n   = 5;
grid_trans_n = 5;

rot_deg_max  = 2;     % 原来 5 -> 2
trans_max_px = 4;     % 原来 8 -> 4

use_signed_trans = false;   % false: [0..max], true: [-max..max]
frames_used      = 12;      % 消融抽样帧数
rng_seed         = 0;

% ======= 粗/细配准（初始化） =======
pyr_levels_coarse = 2;
imreg_affine_pyr  = 3;

% ======= SE(2) refine（论文式参数） =======
pyr_levels_refine = 3;
iters_per_lvl     = 10;
sample_stride     = 3;
huber_delta       = 1.5;
lambda_mag        = 0.30;
lambda_dir        = 0.20;
stop_eps          = 1e-3;
step_clip         = 5;
damping           = 1e-6;

% ======= Success 判据（可按你的数据调） =======
ncc_success_th   = 0.70;
energy_drop_min  = 1e-6;
min_valid_pix    = 1200;

% ======= 输出控制 =======
do_postfusion_edgestrength = false;   % 默认关（避免跑太久）
save_outputs = false;                 % true 会保存 png/csv
output_dir   = fullfile(pwd, 'ablation_outputs');
% ===================================

rng(rng_seed);

%% A) 读取窗口帧
[newfilename, multi_frames] = read_bin_window(bin_dir, target_name, W, H, start_frame, window_len);
N = size(multi_frames,3);
center_idx = ceil(N/2);
Iref = multi_frames(:,:,center_idx);

fprintf('Loaded: %s\n', newfilename);
fprintf('Window: [%d..%d], N=%d, center=%d\n', start_frame, start_frame+N-1, N, center_idx);

%% B) 选择用于消融的帧（均匀抽样）
idx_all = setdiff(1:N, center_idx);
if frames_used >= numel(idx_all)
    idx_use = idx_all(:)';
else
    idx_use = round(linspace(1, numel(idx_all), frames_used));
    idx_use = idx_all(idx_use);
end
fprintf('Ablation frames used: %d\n', numel(idx_use));

%% C) 扰动网格
if use_signed_trans
    trans_grid = linspace(-trans_max_px, trans_max_px, grid_trans_n);
else
    trans_grid = linspace(0, trans_max_px, grid_trans_n);
end
rot_grid = linspace(0, rot_deg_max, grid_rot_n);

fprintf('Grid trans: [%s], rot(deg): [%s]\n', num2str(trans_grid), num2str(rot_grid));
fprintf('Modes: NoRefine / Additive / Lie\n\n');

%% D) 扫描：每个帧 -> init -> 加扰动 -> 3模式 refine（第一遍：只统计，不记录曲线）
modes = {'NoRefine','Additive','Lie'};
S = struct();
for mi=1:numel(modes)
    S.(modes{mi}) = init_stats_container(numel(idx_use), numel(rot_grid), numel(trans_grid));
end

for fi = 1:numel(idx_use)
    k = idx_use(fi);
    Imov = multi_frames(:,:,k);

    % (1) 初始化：粗+affine
    Gm2r_init = estimate_init_affine_m2r(Imov, Iref, pyr_levels_coarse, imreg_affine_pyr);
    H_r2m_init = inv(Gm2r_init); % ref->mov

    for ri = 1:numel(rot_grid)
        for ti = 1:numel(trans_grid)
            dth = deg2rad(rot_grid(ri));
            dt  = trans_grid(ti);

            % (2) 人为加扰动（ref->mov）
            H0 = se2_exp(dth, dt, 0.0) * H_r2m_init;

            for mi=1:numel(modes)
                mode = modes{mi};

                tStart = tic;
                switch mode
                    case 'NoRefine'
                        Hhat = H0;
                        it_used = 0;
                        [E0, ncc0, valid0] = eval_energy_ncc(Imov, Iref, H0, ...
                            lambda_mag, lambda_dir, huber_delta, sample_stride, min_valid_pix);
                        E1 = E0; ncc1 = ncc0; valid1 = valid0;

                    case 'Additive'
                        [Hhat, it_used, E0, E1, ncc0, ncc1, valid1] = ...
                            refine_se2_additive(Imov, Iref, H0, ...
                                pyr_levels_refine, iters_per_lvl, sample_stride, ...
                                lambda_mag, lambda_dir, huber_delta, ...
                                stop_eps, step_clip, damping, min_valid_pix, false);

                    case 'Lie'
                        [Hhat, it_used, E0, E1, ncc0, ncc1, valid1] = ...
                            refine_se2_lie(Imov, Iref, H0, ...
                                pyr_levels_refine, iters_per_lvl, sample_stride, ...
                                lambda_mag, lambda_dir, huber_delta, ...
                                stop_eps, step_clip, damping, min_valid_pix, false);
                end
                tCost = toc(tStart);

                % success 判据：有效域足够 + NCC够高 + 能量下降
                succ = valid1 && isfinite(ncc1) && (ncc1 >= ncc_success_th) && ...
                       isfinite(E0) && isfinite(E1) && ((E0 - E1) >= energy_drop_min);

                S.(mode) = update_stats(S.(mode), fi, ri, ti, succ, it_used, tCost, E1, ncc1);
            end
        end
    end

    fprintf('Frame %d/%d done (k=%d)\n', fi, numel(idx_use), k);
end

%% E) 打印消融统计（命令行）
fprintf('\n=== Ablation sweep (perturb init) ===\n');
fprintf('Ablation frames used: %d\n', numel(idx_use));
fprintf('Grid: rot=%d x trans=%d, modes=%d\n\n', numel(rot_grid), numel(trans_grid), numel(modes));

fprintf('=== Ablation Summary (global averages over grid) ===\n');
Rsum = struct();
for mi=1:numel(modes)
    mode = modes{mi};
    R = summarize_stats(S.(mode));
    Rsum.(mode) = R;
    fprintf('[%s] SuccessRate=%.3f | MeanIters=%.2f | MeanTime=%.3fs | MeanEnergy=%.6g | MeanNCC=%.4f\n', ...
        mode, R.success_rate, R.mean_iters, R.mean_time, R.mean_energy, R.mean_ncc);
end

%% ========== 输出 1：图1 Success Map 热力图 ==========
fig1 = figure('Name','Fig1 Success Map (Δt-Δθ grid)', 'Color','w', 'Position',[80 80 1400 420]);
for mi=1:numel(modes)
    mode = modes{mi};
    succ_map = squeeze(mean(S.(mode).succ, 1));  % [nR x nT] 对帧平均
    subplot(1,3,mi);
    imagesc(trans_grid, rot_grid, succ_map);
    axis image; axis xy;
    xlabel('\Delta t (px)');
    ylabel('\Delta \theta (deg)');
    title(sprintf('%s Success Map', mode));
    colorbar;
    clim([0 1]);
end
sgtitle(sprintf('Success Rate Heatmaps | rot\\_max=%.2fdeg, trans\\_max=%.2fpx', rot_deg_max, trans_max_px));

%% ========== 输出 2：表1 全局统计表 ==========
T = table( ...
    cellstr(string(modes(:))), ...
    [Rsum.NoRefine.success_rate; Rsum.Additive.success_rate; Rsum.Lie.success_rate], ...
    [Rsum.NoRefine.mean_iters;   Rsum.Additive.mean_iters;   Rsum.Lie.mean_iters], ...
    [Rsum.NoRefine.mean_time;    Rsum.Additive.mean_time;    Rsum.Lie.mean_time], ...
    [Rsum.NoRefine.mean_energy;  Rsum.Additive.mean_energy;  Rsum.Lie.mean_energy], ...
    [Rsum.NoRefine.mean_ncc;     Rsum.Additive.mean_ncc;     Rsum.Lie.mean_ncc], ...
    'VariableNames', {'Mode','SuccessRate','MeanIters','MeanTime_s','MeanEnergy','MeanNCC'});

disp(' ');
disp('=== Table1 Global Stats (moderate perturb grid) ===');
disp(T);

figT = figure('Name','Table1 Global Stats', 'Color','w', 'Position',[200 200 820 220]);

% uitable 兼容：把 string/char 都变成 cellstr + numeric
Tc = table2cell(T);
Tc(:,1) = cellstr(T.Mode);  % 确保第一列是 char cell
uitable(figT, 'Data', Tc, 'ColumnName', T.Properties.VariableNames, ...
    'Units','normalized', 'Position',[0 0 1 1]);

%% ========== 输出 3：图2 收敛曲线（自动选代表点 + 二次复跑，保证不空） ==========
[rep, okrep] = pick_representative_point(S, rot_grid, trans_grid);
if ~okrep
    warning('No suitable representative point found for BOTH Additive & Lie success. Fig2 will try Lie-only point.');
end

% 选择一个帧，使该点下（优先：Additive&Lie都成功；否则：Lie成功）
[rep_fi, rep_ok_frame] = pick_representative_frame(S, rep.ri, rep.ti, okrep);
if ~rep_ok_frame
    warning('No representative frame found. Skip Fig2.');
else
    k_rep = idx_use(rep_fi);
    Imov = multi_frames(:,:,k_rep);

    % 重新计算 init（与 sweep 一致）
    Gm2r_init = estimate_init_affine_m2r(Imov, Iref, pyr_levels_coarse, imreg_affine_pyr);
    H_r2m_init = inv(Gm2r_init); % ref->mov

    dth = deg2rad(rot_grid(rep.ri));
    dt  = trans_grid(rep.ti);
    H0  = se2_exp(dth, dt, 0.0) * H_r2m_init;

    % 复跑：记录历史
    [H_add,~,E0A,E1A,N0A,N1A,validA,E_hist_add,N_hist_add] = ...
        refine_se2_additive(Imov, Iref, H0, ...
            pyr_levels_refine, iters_per_lvl, sample_stride, ...
            lambda_mag, lambda_dir, huber_delta, ...
            stop_eps, step_clip, damping, min_valid_pix, true); %#ok<ASGLU>

    [H_lie,~,E0L,E1L,N0L,N1L,validL,E_hist_lie,N_hist_lie] = ...
        refine_se2_lie(Imov, Iref, H0, ...
            pyr_levels_refine, iters_per_lvl, sample_stride, ...
            lambda_mag, lambda_dir, huber_delta, ...
            stop_eps, step_clip, damping, min_valid_pix, true); %#ok<ASGLU>

    % 安全保护：至少有2点
    if numel(E_hist_add) < 2, E_hist_add = [E0A; E1A]; end
    if numel(E_hist_lie) < 2, E_hist_lie = [E0L; E1L]; end
    if numel(N_hist_add) < 2, N_hist_add = [N0A; N1A]; end
    if numel(N_hist_lie) < 2, N_hist_lie = [N0L; N1L]; end

    fig2 = figure('Name','Fig2 Convergence Curves (Representative point)', 'Color','w', 'Position',[120 120 1100 420]);

    % 能量归一化
    Eadd = E_hist_add / max(E_hist_add(1), eps);
    Elie = E_hist_lie / max(E_hist_lie(1), eps);

    subplot(1,2,1);
    plot(1:numel(Eadd), Eadd, '-o'); hold on;
    plot(1:numel(Elie), Elie, '-o'); grid on;
    xlabel('Iteration');
    ylabel('Normalized Energy (E/E_0)');
    title('Energy convergence');
    legend('Additive','Lie','Location','northeast');

    subplot(1,2,2);
    plot(1:numel(N_hist_add), N_hist_add, '-o'); hold on;
    plot(1:numel(N_hist_lie), N_hist_lie, '-o'); grid on;
    xlabel('Iteration');
    ylabel('NCC');
    title('NCC convergence');
    legend('Additive','Lie','Location','southeast');

    sgtitle(sprintf('Representative point: rot=%.2fdeg, trans=%.2fpx (frame k=%d)', ...
        rot_grid(rep.ri), trans_grid(rep.ti), k_rep));

    % 可选：代表点 warp 对比（init / Additive / Lie）
    Iw0  = warp_by_Gm2r(Imov, inv(H0));
    IwA  = warp_by_Gm2r(Imov, inv(H_add));
    IwL  = warp_by_Gm2r(Imov, inv(H_lie));

    fig3 = figure('Name','Representative warp compare', 'Color','w', 'Position',[120 580 1400 520]);
    subplot(2,3,1); imshow(robust_vis(Iref)); title('Ref (center)');
    subplot(2,3,2); imshow(robust_vis(Imov)); title(sprintf('Mov (k=%d)',k_rep));
    subplot(2,3,3); imshow(robust_vis(Iw0));  title('Warp (perturbed init)');

    subplot(2,3,4); imshow(robust_vis(Iref-IwA)); title('Residual (ref - AdditiveWarp)');
    subplot(2,3,5); imshow(robust_vis(Iref-IwL)); title('Residual (ref - LieWarp)');
    subplot(2,3,6); imshow(robust_vis(IwL));  title('Warp after Lie');
end

%% 保存输出（可选）
if save_outputs
    if ~exist(output_dir,'dir'), mkdir(output_dir); end
    exportgraphics(fig1, fullfile(output_dir, 'Fig1_SuccessMap.png'));
    exportgraphics(figT, fullfile(output_dir, 'Table1_GlobalStats.png'));
    writetable(T, fullfile(output_dir, 'Table1_GlobalStats.csv'));

    if exist('fig2','var') && isvalid(fig2), exportgraphics(fig2, fullfile(output_dir, 'Fig2_Convergence.png')); end
    if exist('fig3','var') && isvalid(fig3), exportgraphics(fig3, fullfile(output_dir, 'Fig3_WarpCompare.png')); end

    fprintf('Saved outputs to: %s\n', output_dir);
end

fprintf('\nDone.\n');
end

%% ============================================================
%% ======================== 代表点选择 =========================
%% ============================================================
function [rep, ok] = pick_representative_point(S, rot_grid, trans_grid)
% 目标：找一个“Lie 成功率高且与 Additive 差异大”的点，并要求两者都有成功（用于画对比曲线）
succA = squeeze(mean(S.Additive.succ, 1)); % [nR x nT]
succL = squeeze(mean(S.Lie.succ, 1));      % [nR x nT]

mask_both = (succA > 0) & (succL > 0);
ok = any(mask_both(:));

if ok
    score = (succL - succA);
    score(~mask_both) = -inf;
    [~, idx] = max(score(:));
else
    % fallback：只要 Lie 有成功
    mask_lie = (succL > 0);
    if any(mask_lie(:))
        score = succL;
        score(~mask_lie) = -inf;
        [~, idx] = max(score(:));
    else
        % 全部失败：硬选右上角
        idx = numel(rot_grid) * numel(trans_grid);
    end
end

[ri, ti] = ind2sub(size(succL), idx);
rep.ri = ri;
rep.ti = ti;
end

function [fi_pick, ok] = pick_representative_frame(S, ri, ti, require_both)
% 优先：该点下 Additive&Lie 都成功的帧；否则：Lie 成功的帧
succA_f = squeeze(S.Additive.succ(:,ri,ti));
succL_f = squeeze(S.Lie.succ(:,ri,ti));

if require_both
    idx = find(succA_f & succL_f, 1, 'first');
    if ~isempty(idx)
        fi_pick = idx; ok = true; return;
    end
end

idx = find(succL_f, 1, 'first');
if ~isempty(idx)
    fi_pick = idx; ok = true; return;
end

% 再退一步：随便选第一帧（可能失败，但不会崩）
fi_pick = 1;
ok = false;
end

%% ============================================================
%% ======================== 统计容器 ===========================
%% ============================================================
function T = init_stats_container(nF, nR, nT)
T.succ  = false(nF,nR,nT);
T.iters = zeros(nF,nR,nT);
T.time  = zeros(nF,nR,nT);
T.energy= nan(nF,nR,nT);
T.ncc   = nan(nF,nR,nT);
end

function T = update_stats(T, fi, ri, ti, succ, it_used, tCost, E1, ncc1)
T.succ(fi,ri,ti)   = succ;
T.iters(fi,ri,ti)  = it_used;
T.time(fi,ri,ti)   = tCost;
T.energy(fi,ri,ti) = E1;
T.ncc(fi,ri,ti)    = ncc1;
end

function R = summarize_stats(T)
succ = T.succ(:);
R.success_rate = mean(succ);

if any(succ)
    R.mean_iters  = mean(T.iters(succ));
    R.mean_time   = mean(T.time(succ));
    R.mean_energy = mean(T.energy(succ),'omitnan');
    R.mean_ncc    = mean(T.ncc(succ),'omitnan');
else
    R.mean_iters  = NaN;
    R.mean_time   = NaN;
    R.mean_energy = NaN;
    R.mean_ncc    = NaN;
end
end

%% ============================================================
%% =============== 初始化：粗配准 + affine ======================
%% ============================================================
function Gm2r_col = estimate_init_affine_m2r(Imov_raw, Iref_raw, pyr_levels_coarse, imreg_affine_pyr)
Iref = mat2gray(Iref_raw);
Imov = mat2gray(Imov_raw);

% (1) 粗：imregcorr rigid 或 phasecorr 平移兜底
tform0 = affine2d(eye(3));
if exist('imregcorr','file')==2
    try
        tform0 = imregcorr(Imov, Iref, 'rigid');
    catch
        [dx,dy] = estimate_translation_phasecorr(Imov, Iref, pyr_levels_coarse);
        tform0 = affine2d([1 0 0; 0 1 0; dx dy 1]);
    end
else
    [dx,dy] = estimate_translation_phasecorr(Imov, Iref, pyr_levels_coarse);
    tform0 = affine2d([1 0 0; 0 1 0; dx dy 1]);
end

% (2) 细：imregtform affine
tform_aff = tform0;
if exist('imregtform','file')==2
    try
        [opt, met] = imregconfig('monomodal');
        opt.MaximumIterations = 80;
        opt.MinimumStepLength = 1e-4;
        opt.RelaxationFactor  = 0.7;
        tform_aff = imregtform(Imov, Iref, 'affine', opt, met, ...
            'InitialTransformation', tform0, 'PyramidLevels', imreg_affine_pyr);
    catch
        tform_aff = tform0;
    end
end

% MATLAB row-convention -> column convention
Gm2r_col = (tform_aff.T)'; % p_ref = G * p_mov
Gm2r_col(3,:) = [0 0 1];
Gm2r_col(:,3) = [Gm2r_col(1,3); Gm2r_col(2,3); 1];
end

%% ============================================================
%% ================= Additive refine (对照组) ===================
%% ============================================================
function [H_r2m, it_used, E0, E1, ncc0, ncc1, valid1, E_hist, ncc_hist] = refine_se2_additive(Imov_raw, Iref_raw, H0, ...
    pyr_levels, iters_per_lvl, sample_stride, ...
    lambda_mag, lambda_dir, huber_delta, ...
    stop_eps, step_clip, damping, min_valid_pix, record_hist)

if nargin < 15, record_hist = false; end
E_hist = []; ncc_hist = [];

H_r2m = H0;

[E0, ncc0, valid0] = eval_energy_ncc(Imov_raw, Iref_raw, H_r2m, lambda_mag, lambda_dir, huber_delta, sample_stride, min_valid_pix);
E1 = E0; ncc1 = ncc0; valid1 = valid0;
it_used = 0;

% ★关键：保证历史不空（即使 0 次迭代也能画）
if record_hist
    E_hist = E0;
    ncc_hist = ncc0;
end

for lvl = pyr_levels:-1:1
    s = 1/(2^(lvl-1));
    Imov = imresize(double(Imov_raw), s);
    Iref = imresize(double(Iref_raw), s);

    [h,w] = size(Iref);
    cx = (w+1)/2; cy = (h+1)/2;

    [Hc, Tc, iTc] = to_centered(H_r2m, cx, cy);

    for it=1:iters_per_lvl
        [d, ok, Ecur, ncccur] = one_gn_step(Imov, Iref, Hc, ...
            lambda_mag, lambda_dir, huber_delta, sample_stride, damping, min_valid_pix);
        if ~ok, break; end

        it_used = it_used + 1;

        nd = norm(d);
        if nd > step_clip
            d = d * (step_clip / nd);
        end

        % Additive：在 centered 坐标内直接加角度/平移
        R = Hc(1:2,1:2);
        t = Hc(1:2,3);

        th_now = atan2(R(2,1), R(1,1));
        th_new = th_now + d(1);

        Rn = [cos(th_new) -sin(th_new); sin(th_new) cos(th_new)];
        tn = t + [d(2); d(3)];
        Hc = [Rn tn; 0 0 1];

        if record_hist
            E_hist(end+1,1) = Ecur; %#ok<AGROW>
            ncc_hist(end+1,1)= ncccur; %#ok<AGROW>
        end

        if norm(d) < stop_eps
            break;
        end
    end

    H_r2m = Tc * Hc * iTc;
end

[E1, ncc1, valid1] = eval_energy_ncc(Imov_raw, Iref_raw, H_r2m, lambda_mag, lambda_dir, huber_delta, sample_stride, min_valid_pix);

% ★关键：末值也写入
if record_hist
    E_hist(end+1,1) = E1;
    ncc_hist(end+1,1)= ncc1;
end
end

%% ============================================================
%% =================== Lie refine (论文组) ======================
%% ============================================================
function [H_r2m, it_used, E0, E1, ncc0, ncc1, valid1, E_hist, ncc_hist] = refine_se2_lie(Imov_raw, Iref_raw, H0, ...
    pyr_levels, iters_per_lvl, sample_stride, ...
    lambda_mag, lambda_dir, huber_delta, ...
    stop_eps, step_clip, damping, min_valid_pix, record_hist)

if nargin < 15, record_hist = false; end
E_hist = []; ncc_hist = [];

H_r2m = H0;

[E0, ncc0, valid0] = eval_energy_ncc(Imov_raw, Iref_raw, H_r2m, lambda_mag, lambda_dir, huber_delta, sample_stride, min_valid_pix);
E1 = E0; ncc1 = ncc0; valid1 = valid0;
it_used = 0;

% ★关键：保证历史不空
if record_hist
    E_hist = E0;
    ncc_hist = ncc0;
end

for lvl = pyr_levels:-1:1
    s = 1/(2^(lvl-1));
    Imov = imresize(double(Imov_raw), s);
    Iref = imresize(double(Iref_raw), s);

    [h,w] = size(Iref);
    cx = (w+1)/2; cy = (h+1)/2;

    [Hc, Tc, iTc] = to_centered(H_r2m, cx, cy);

    for it=1:iters_per_lvl
        [d, ok, Ecur, ncccur] = one_gn_step(Imov, Iref, Hc, ...
            lambda_mag, lambda_dir, huber_delta, sample_stride, damping, min_valid_pix);
        if ~ok, break; end

        it_used = it_used + 1;

        nd = norm(d);
        if nd > step_clip
            d = d * (step_clip / nd);
        end

        % Lie 群复合：H <- Exp(d) * H
        Expd = se2_exp(d(1), d(2), d(3));
        Hc = Expd * Hc;

        if record_hist
            E_hist(end+1,1) = Ecur; %#ok<AGROW>
            ncc_hist(end+1,1)= ncccur; %#ok<AGROW>
        end

        if norm(d) < stop_eps
            break;
        end
    end

    H_r2m = Tc * Hc * iTc;
end

[E1, ncc1, valid1] = eval_energy_ncc(Imov_raw, Iref_raw, H_r2m, lambda_mag, lambda_dir, huber_delta, sample_stride, min_valid_pix);

% ★关键：末值也写入
if record_hist
    E_hist(end+1,1) = E1;
    ncc_hist(end+1,1)= ncc1;
end
end

%% ============================================================
%% ===================== 单步 GN/IRLS ==========================
%% ============================================================
function [d, ok, E, ncc] = one_gn_step(Imov, Iref, Hc, ...
    lambda_mag, lambda_dir, huber_delta, sample_stride, damping, min_valid_pix)

[h,w] = size(Iref);

% ref 梯度
[Ir_x, Ir_y] = grad_xy(Iref);
Ir_mag = hypot(Ir_x, Ir_y);
Ir_ux  = Ir_x ./ (Ir_mag + 1e-9);
Ir_uy  = Ir_y ./ (Ir_mag + 1e-9);

% mov 一阶 + Hessian
[Im_x, Im_y, Im_xx, Im_xy, Im_yy] = grad_hessian(Imov);

% 采样网格（ref）
xs = 1:sample_stride:w;
ys = 1:sample_stride:h;
[X,Y] = meshgrid(xs,ys);
Xv = X(:); Yv = Y(:);

cx = (w+1)/2; cy = (h+1)/2;
x0 = Xv - cx;
y0 = Yv - cy;

R = Hc(1:2,1:2);
t = Hc(1:2,3);

xm0 = R(1,1)*x0 + R(1,2)*y0 + t(1);
ym0 = R(2,1)*x0 + R(2,2)*y0 + t(2);

xm = xm0 + cx;
ym = ym0 + cy;

valid = (xm>=2) & (xm<=w-1) & (ym>=2) & (ym<=h-1);
if nnz(valid) < min_valid_pix
    ok = false; d = [0;0;0]; E = inf; ncc = -inf; return;
end

xv = Xv(valid); yv = Yv(valid);
x0v = x0(valid); y0v = y0(valid);
xm0v = xm0(valid); ym0v = ym0(valid);
xmv = xm(valid); ymv = ym(valid);

% 强度残差
Iw = interp2(Imov, xmv, ymv, 'linear', NaN);
Ir = Iref(sub2ind([h,w], yv, xv));
rI = Iw - Ir;

% moving 梯度 + Hessian 采样
gx_m = interp2(Im_x,  xmv, ymv, 'linear', NaN);
gy_m = interp2(Im_y,  xmv, ymv, 'linear', NaN);
gxx  = interp2(Im_xx, xmv, ymv, 'linear', NaN);
gxy  = interp2(Im_xy, xmv, ymv, 'linear', NaN);
gyy  = interp2(Im_yy, xmv, ymv, 'linear', NaN);

% grad wrt ref: g_w = R^T g_m
Rt = R';
gwx = Rt(1,1)*gx_m + Rt(1,2)*gy_m;
gwy = Rt(2,1)*gx_m + Rt(2,2)*gy_m;
mag_w = hypot(gwx, gwy);

idx = sub2ind([h,w], yv, xv);
mag_r = Ir_mag(idx);
urx   = Ir_ux(idx);
ury   = Ir_uy(idx);

rMag = mag_w - mag_r;

uwx = gwx ./ (mag_w + 1e-9);
uwy = gwy ./ (mag_w + 1e-9);
dotur = uwx.*urx + uwy.*ury;
rDir = 1 - dotur;

r = [rI;
     sqrt(lambda_mag)*rMag;
     sqrt(lambda_dir)*rDir];

% Huber 权重
Wv = huber_weights(r, huber_delta);

% ===== 解析雅可比 =====
Jmat = [0 -1; 1 0];

dxm0_dth = -ym0v;
dym0_dth =  xm0v;

dIw_dth  = gx_m .* dxm0_dth + gy_m .* dym0_dth;
dIw_dtx  = gx_m;
dIw_dty  = gy_m;

dRt_dth = -Rt * Jmat;

dgm_x_th = gxx .* dxm0_dth + gxy .* dym0_dth;
dgm_y_th = gxy .* dxm0_dth + gyy .* dym0_dth;

dgm_x_tx = gxx;  dgm_y_tx = gxy;
dgm_x_ty = gxy;  dgm_y_ty = gyy;

dg_wx_th = dRt_dth(1,1).*gx_m + dRt_dth(1,2).*gy_m + Rt(1,1).*dgm_x_th + Rt(1,2).*dgm_y_th;
dg_wy_th = dRt_dth(2,1).*gx_m + dRt_dth(2,2).*gy_m + Rt(2,1).*dgm_x_th + Rt(2,2).*dgm_y_th;

dg_wx_tx = Rt(1,1).*dgm_x_tx + Rt(1,2).*dgm_y_tx;
dg_wy_tx = Rt(2,1).*dgm_x_tx + Rt(2,2).*dgm_y_tx;

dg_wx_ty = Rt(1,1).*dgm_x_ty + Rt(1,2).*dgm_y_ty;
dg_wy_ty = Rt(2,1).*dgm_x_ty + Rt(2,2).*dgm_y_ty;

dMag_dth = uwx .* dg_wx_th + uwy .* dg_wy_th;
dMag_dtx = uwx .* dg_wx_tx + uwy .* dg_wy_tx;
dMag_dty = uwx .* dg_wx_ty + uwy .* dg_wy_ty;

invMag = 1 ./ (mag_w + 1e-9);
P11 = 1 - uwx.*uwx;  P12 = -uwx.*uwy;
P21 = P12;           P22 = 1 - uwy.*uwy;

dux_th = (P11.*dg_wx_th + P12.*dg_wy_th) .* invMag;
duy_th = (P21.*dg_wx_th + P22.*dg_wy_th) .* invMag;

dux_tx = (P11.*dg_wx_tx + P12.*dg_wy_tx) .* invMag;
duy_tx = (P21.*dg_wx_tx + P22.*dg_wy_tx) .* invMag;

dux_ty = (P11.*dg_wx_ty + P12.*dg_wy_ty) .* invMag;
duy_ty = (P21.*dg_wx_ty + P22.*dg_wy_ty) .* invMag;

dDir_dth = -(urx .* dux_th + ury .* duy_th);
dDir_dtx = -(urx .* dux_tx + ury .* duy_tx);
dDir_dty = -(urx .* dux_ty + ury .* duy_ty);

JI = [dIw_dth, dIw_dtx, dIw_dty];
JM = sqrt(lambda_mag) * [dMag_dth, dMag_dtx, dMag_dty];
JD = sqrt(lambda_dir) * [dDir_dth, dDir_dtx, dDir_dty];
J  = [JI; JM; JD];

JW = J .* Wv;
A  = J' * JW + damping*eye(3);
b  = -(J' * (Wv .* r));
d  = A \ b;

E = sum(Wv .* (r.^2));
ncc = ncc_on_samples(Ir, Iw);

ok = all(isfinite(d)) && isfinite(E) && isfinite(ncc);
end

%% ============================================================
%% ==================== 能量 / NCC 评估 ========================
%% ============================================================
function [E, ncc, valid_ok] = eval_energy_ncc(Imov_raw, Iref_raw, H_r2m, ...
    lambda_mag, lambda_dir, huber_delta, sample_stride, min_valid_pix)

Imov = double(Imov_raw);
Iref = double(Iref_raw);

[h,w] = size(Iref);
cx = (w+1)/2; cy = (h+1)/2;

[Hc, ~, ~] = to_centered(H_r2m, cx, cy);

[Ir_x, Ir_y] = grad_xy(Iref);
Ir_mag = hypot(Ir_x, Ir_y);
Ir_ux  = Ir_x ./ (Ir_mag + 1e-9);
Ir_uy  = Ir_y ./ (Ir_mag + 1e-9);

[Im_x, Im_y] = grad_xy(Imov);

xs = 1:sample_stride:w;
ys = 1:sample_stride:h;
[X,Y] = meshgrid(xs,ys);
Xv = X(:); Yv = Y(:);

x0 = Xv - cx;
y0 = Yv - cy;

R = Hc(1:2,1:2);
t = Hc(1:2,3);

xm0 = R(1,1)*x0 + R(1,2)*y0 + t(1);
ym0 = R(2,1)*x0 + R(2,2)*y0 + t(2);

xm = xm0 + cx;
ym = ym0 + cy;

valid = (xm>=2) & (xm<=w-1) & (ym>=2) & (ym<=h-1);
if nnz(valid) < min_valid_pix
    E = inf; ncc = -inf; valid_ok = false; return;
end

xv  = Xv(valid); yv = Yv(valid);
xmv = xm(valid); ymv = ym(valid);

Iw = interp2(Imov, xmv, ymv, 'linear', NaN);
Ir = Iref(sub2ind([h,w], yv, xv));
rI = Iw - Ir;

gx_m = interp2(Im_x, xmv, ymv, 'linear', NaN);
gy_m = interp2(Im_y, xmv, ymv, 'linear', NaN);

Rt = R';
gwx = Rt(1,1)*gx_m + Rt(1,2)*gy_m;
gwy = Rt(2,1)*gx_m + Rt(2,2)*gy_m;
mag_w = hypot(gwx, gwy);

idx = sub2ind([h,w], yv, xv);
mag_r = Ir_mag(idx);
urx = Ir_ux(idx); ury = Ir_uy(idx);

rMag = mag_w - mag_r;

uwx = gwx ./ (mag_w + 1e-9);
uwy = gwy ./ (mag_w + 1e-9);
dotur = uwx.*urx + uwy.*ury;
rDir = 1 - dotur;

r = [rI;
     sqrt(lambda_mag)*rMag;
     sqrt(lambda_dir)*rDir];

Wv = huber_weights(r, huber_delta);
E = sum(Wv .* (r.^2));
ncc = ncc_on_samples(Ir, Iw);
valid_ok = true;
end

function v = ncc_on_samples(a,b)
a = double(a(:)); b = double(b(:));
m = isfinite(a) & isfinite(b);
a = a(m); b = b(m);
if numel(a) < 200
    v = -inf; return;
end
a = a - mean(a); b = b - mean(b);
v = (a'*b) / (sqrt(sum(a.^2))*sqrt(sum(b.^2)) + 1e-12);
end

%% ============================================================
%% ======================= warp 工具 ===========================
%% ============================================================
function Iw = warp_by_Gm2r(Imov_raw, Gm2r_col)
[h,w] = size(Imov_raw);
Rout = imref2d([h,w]);
tform = affine2d(Gm2r_col'); % row convention
fillv = median(Imov_raw(:),'omitnan');
Iw = imwarp(Imov_raw, tform, 'OutputView', Rout, 'FillValues', fillv);
end

function [Hc, Tc, iTc] = to_centered(H_abs, cx, cy)
Tc  = [1 0 cx; 0 1 cy; 0 0 1];
iTc = [1 0 -cx; 0 1 -cy; 0 0 1];
Hc = iTc * H_abs * Tc;
end

%% ============================================================
%% ===================== 基础算子/工具 ==========================
%% ============================================================
function [filename, frames] = read_bin_window(bin_dir, target_name, W, H, start_frame, window_len)
if ~exist(bin_dir,'dir'), error('目录不存在：%s', bin_dir); end
if ~isempty(target_name) && exist(fullfile(bin_dir,target_name),'file')
    filename = fullfile(bin_dir,target_name);
else
    bins = dir(fullfile(bin_dir,'*.bin'));
    if isempty(bins), error('目录中未找到 .bin：%s', bin_dir); end
    [~,idx] = max([bins.datenum]);
    filename = fullfile(bin_dir, bins(idx).name);
    fprintf('未指定文件名，自动选择最新的：%s\n', filename);
end

pix_per_frame = W*H;
bytes_per_frame = pix_per_frame*2;
info = dir(filename);
n_frames_total = floor(info.bytes / bytes_per_frame);
if n_frames_total <= 0, error('文件像素数不足一帧：%s', filename); end

end_frame = start_frame + window_len - 1;
if end_frame > n_frames_total
    warning('请求区间 [%d,%d] 超过总帧数 %d，回退到末尾窗口。', start_frame, end_frame, n_frames_total);
    end_frame = n_frames_total;
    start_frame = max(1, end_frame-window_len+1);
end
takeN = end_frame-start_frame+1;
fprintf('文件共 %d 帧；读取 [%d,%d]（%d帧）。\n', n_frames_total, start_frame, end_frame, takeN);

offset_bytes = (start_frame-1)*bytes_per_frame;
read_count = takeN*pix_per_frame;

fid = fopen(filename,'rb');
if fid<0, error('无法打开：%s', filename); end
cobj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fseek(fid, offset_bytes, 'bof');
raw = fread(fid, read_count, 'uint16=>single');
raw = raw(1:min(numel(raw),read_count));
frames = reshape(raw, H, W, takeN);
end

function [dx, dy] = estimate_translation_phasecorr(moving, ref, pyr_levels)
moving = double(moving); ref = double(ref);
dx = 0; dy = 0;
for lvl = pyr_levels:-1:1
    s = 1/(2^(lvl-1));
    mv = imresize(moving, s, 'bilinear');
    rf = imresize(ref,    s, 'bilinear');

    mv_w = imtranslate(mv, [dx*s, dy*s], 'FillValues', median(mv(:)));

    F1 = fft2(rf); F2 = fft2(mv_w);
    R  = F1 .* conj(F2);
    R  = R ./ (abs(R) + eps);
    c  = real(ifft2(R));

    [rows, cols] = size(rf);
    [py, px] = find(c == max(c(:)), 1, 'first');
    midX = floor(cols/2) + 1;
    midY = floor(rows/2) + 1;
    ddx = px - midX;
    ddy = py - midY;

    dx = dx + ddx / s;
    dy = dy + ddy / s;
end
end

function [Ix, Iy] = grad_xy(I)
I = double(I);
kx = [-0.5 0 0.5];
ky = kx';
Ix = conv2(I, kx, 'same');
Iy = conv2(I, ky, 'same');
end

function [Ix, Iy, Ixx, Ixy, Iyy] = grad_hessian(I)
I = double(I);
kx = [-0.5 0 0.5];
ky = kx';
Ix = conv2(I, kx, 'same');
Iy = conv2(I, ky, 'same');

kxx = [1 -2 1];
kyy = kxx';
Ixx = conv2(I, kxx, 'same');
Iyy = conv2(I, kyy, 'same');

Ixy = conv2(Ix, ky, 'same');
end

function W = huber_weights(r, delta)
a = abs(r);
W = ones(size(r));
idx = a > delta;
W(idx) = delta ./ (a(idx) + eps);
end

function Expd = se2_exp(dth, dtx, dty)
J = [0 -1; 1 0];
if abs(dth) < 1e-8
    R = eye(2);
    V = eye(2) + 0.5*dth*J;
else
    c = cos(dth); s = sin(dth);
    R = [c -s; s c];
    A = s/dth;
    B = (1-c)/dth;
    V = A*eye(2) + B*J;
end
rho = [dtx; dty];
t = V * rho;
Expd = [R t; 0 0 1];
end

function I = robust_vis(img)
img = double(img);
vals = img(isfinite(img));
if isempty(vals), vals = img(:); end
lo = prctile(vals, 0.5);
hi = prctile(vals, 99.5);
if ~isfinite(lo), lo = min(vals); end
if ~isfinite(hi), hi = max(vals); end
if hi <= lo, hi = lo + 1; end
I = (img - lo) / (hi - lo);
I = min(max(I,0),1);
end
