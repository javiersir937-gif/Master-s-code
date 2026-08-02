%% chapter4_riblet_validation_bomb.m
% Chapter 4: STL-based mid-body extraction + per-Mach wall-data plotting
% Adds:
%  - Prints run parameters to Command Window (and saves a .txt report)
%  - Creates STL-aligned wall CSV templates for each Mach case (x, R, circumference + NaNs)
%  - Generates separate figures per Mach (tau_w, Cf, nu_w, u_tau, sPlus if present)

clear; clc; close all;

%% ===================== USER SETTINGS =====================
cfg = struct();

cfg.objectName  = "Bomb";

% STL (match your uploaded filename/case)
cfg.stlFile     = "testd.STL";
cfg.scaleFactor = 1.0;       % if STL units are meters, set 1e3 (to convert to mm)

% Radius profile extraction
cfg.nbins         = 250;
cfg.outlierHiPrct = 90;

% Mid-body selection
cfg.xFracRange    = [0.20 0.85];
cfg.radiusTolFrac = 0.12;

% Mach cases + wall CSV names
cfg.machCases = 0.80:0.05:1.40;
cfg.wallFiles = buildWallFileList("wall_bomb_", cfg.machCases);

% Wall coordinate mode:
%  "axial" expects a column x (already axial)
%  "xyz_project" expects X,Y,Z columns and projects onto STL axis
cfg.wallCoordMode = "axial";


% Wall x scaling (templates are in mm; set 1000 if your CFD exports x in meters)
cfg.wallScaleFactor = 1.0;

% Console output control (user requested: ONLY p, q, tau_w, s+, h+ per Mach)
cfg.printRunConditionsToConsole = false;
cfg.printInfoLines = false;

% Riblet sizing in wall units (unitless)
cfg.riblet_sPlus0 = 15.0;
cfg.riblet_hPlus0 = 7.5;
cfg.riblet_subsonicK = 0.05;   % mild compressibility scaling
cfg.riblet_shockK    = 0.25;   % Mach>=1 shock/BL interaction proxy strength
cfg.riblet_maxShockIndex = 5.0;

% Plot control
cfg.plotFullDiagnostics = false; % keep legacy plots available but off by default

% With many Mach increments, overlays get busy
cfg.makeCombinedOverlays = false;

% Outputs
cfg.exportProfileCSV = true;
cfg.profileCSVName   = "bomb_midbody_profile.csv";

cfg.saveRunReport = true;
cfg.reportFile    = "run_conditions_bomb.txt";

cfg.saveFigs = true;
cfg.figDir   = "figs_bomb";
% NEW: STL-aligned wall template generation
cfg.createWallTemplates = true;
cfg.overwriteTemplates  = true;  % WARNING: overwrites existing wall_*.csv files

%% ===================== RUN =====================
run_validation(cfg);

%% ===================== MAIN =====================
function run_validation(cfg)

    if cfg.saveFigs && ~isfolder(cfg.figDir)
        mkdir(cfg.figDir);
    end

    % ---- STEP 1: STL -> robust profile ----
    [V,~] = readSTL_any(cfg.stlFile);
    V = V .* cfg.scaleFactor;

    prof = extractBodyProfileFromSTL(V, cfg.nbins, cfg.outlierHiPrct);

    % Mid-body selection
    mid = selectMidBody(prof, cfg.xFracRange, cfg.radiusTolFrac);
    if isempty(mid.binIdx)
        error("Mid-body selection returned empty. Try widening xFracRange or increasing radiusTolFrac.");
    end

    % Export mid-body geometry CSV (useful for CAD pattern layout)
    if cfg.exportProfileCSV
        Tgeom = table(mid.x(:), mid.R(:), 2*pi*mid.R(:), ...
            'VariableNames', {'x_axial_mm','R_mm','circumference_mm'});
        writetable(Tgeom, cfg.profileCSVName);
    end

    % NEW: Create STL-aligned wall templates for each Mach case
    if cfg.createWallTemplates
        writeWallTemplatesFromSTL(mid, cfg);
    end

    % Plot STL radius profile + mid-body bins
    figure('Name', cfg.objectName + " STL Radius Profile", 'NumberTitle','off');
    grid on; hold on;
    plot(prof.xc, prof.Rb, '-', 'LineWidth', 1.6);
    plot(mid.x, mid.R, 'o', 'LineWidth', 1.0);
    xlabel('x along inferred principal axis (mm)');
    ylabel('robust radius R(x) (mm)');
    legend('Robust radius profile','Mid-body bins','Location','best');
    titleBelowXLabel(gca, cfg.objectName + ": STL-derived robust radius profile");
    if cfg.saveFigs
        exportgraphics(gcf, fullfile(cfg.figDir, cfg.objectName + "_stl_radius_profile.png"));
    end

    % ---- STEP 2: Print run parameters ----
    printRunConditions(cfg, prof, mid);

    % ---- STEP 3: Per Mach case sweep ----
    results = struct([]);

    % User requested: print ONLY p, q, tau_w, s+, h+ per Mach increment
    for i = 1:numel(cfg.machCases)
        M = cfg.machCases(i);
        f = string(cfg.wallFiles(i));

        % Ensure a CSV exists for every Mach (create STL-aligned template if missing)
        if ~isfile(f) && cfg.createWallTemplates
            writeSingleWallTemplateFromSTL(mid, f, cfg);
        end

        if isfile(f)
            wall = readWallCSV(f, cfg.wallScaleFactor);
        else
            wall = struct('x',[],'X',[],'Y',[],'Z',[],'p_w',[],'q_inf',[],'tau_w',[],'cf',[],'rho_w',[],'mu_w',[],'nu_w',[],'sPlus',[],'hPlus',[]);
        end

        % Ensure an axial coordinate exists (x_axial) if we actually have wall points
        if (~isempty(wall.x) || (~isempty(wall.X)&&~isempty(wall.Y)&&~isempty(wall.Z)))
            wall.x_axial = computeWallAxialX(wall, cfg.wallCoordMode, prof);
        else
            wall.x_axial = [];
        end

        % Aggregate wall data into STL mid-body bins
        agg = aggregateWallToBins(wall, prof, mid, "median");

        % Derived fields
        d = computeDerivedFields(agg);

        % Build/patch pressure + dynamic pressure + tau_w if missing
        fs = freestreamSeaLevel(M);
        [p_use, q_use, tau_use] = fillPQTauProfiles(agg, d, fs, M);

        % Compute unitless riblet sizing (wall units) and allow CSV-provided overrides
        rib = validatedRibletsUnitless(M, agg.x, p_use, cfg);
        sPlus = rib.sPlus;
        hPlus = rib.hPlus;
        if any(~isnan(agg.sPlus))
            sPlus = agg.sPlus;
        end
        if isfield(agg,'hPlus') && any(~isnan(agg.hPlus))
            hPlus = agg.hPlus;
        end

        % Print ONLY requested values (use medians to give one value per Mach)
        fprintf('Mach %.2f | p=%.4g Pa | q=%.4g Pa | tau_w=%.4g Pa | s_plus=%.4g | h_plus=%.4g\n', ...
            M, medOrNaN(p_use), medOrNaN(q_use), medOrNaN(tau_use), medOrNaN(sPlus), medOrNaN(hPlus));

        % Plot (one figure per Mach increment)
        plotMachFigures(cfg, M, agg.x, p_use, tau_use, sPlus, hPlus);

        % Store for optional overlays
        results(end+1).Mach = M; %#ok<AGROW>
        results(end).x      = agg.x;
        results(end).tau_w  = tau_use;
        results(end).sPlus  = sPlus;
        results(end).hPlus  = hPlus;
        % medians for end-of-run vs-Mach comparison plots
        results(end).p_med      = medOrNaN(p_use);
        results(end).q_med      = medOrNaN(q_use);
        results(end).tau_med    = medOrNaN(tau_use);
        results(end).sPlus_med  = medOrNaN(sPlus);
        results(end).hPlus_med  = medOrNaN(hPlus);
    end

    if cfg.makeCombinedOverlays && ~isempty(results)
        plotOverlays(cfg, results);
    end

    % End-of-run: compare requested variables vs Mach (medians over mid-body)
    if ~isempty(results)
        plotVsMach(cfg, results);
    end
end

%% ===================== NEW: WALL TEMPLATE WRITER =====================
function writeWallTemplatesFromSTL(mid, cfg)
% Creates/overwrites wall CSV files using STL-derived mid-body x/R/circumference (in mm).
% Flow columns are NaN placeholders meant to be replaced by Fluent/CFD exports.

    for i = 1:numel(cfg.machCases)
        f = string(cfg.wallFiles(i));
        if strlength(f)==0
            continue;
        end

        if isfile(f) && ~cfg.overwriteTemplates
            if ~isfield(cfg,'printInfoLines') || cfg.printInfoLines
                fprintf('[INFO] Wall file exists (not overwritten): %s\n', f);
            end
            continue;
        end

        writeSingleWallTemplateFromSTL(mid, f, cfg);
    end
end

function writeSingleWallTemplateFromSTL(mid, filename, cfg)
% Writes a single STL-aligned wall template CSV for a given Mach file name.

    T = table();
    T.x_mm = mid.x(:);
    T.R_mm = mid.R(:);
    T.circumference_mm = 2*pi*mid.R(:);

    n = height(T);

    % placeholders for CFD / theory values
    T.p_w   = nan(n,1);
    T.q_inf = nan(n,1);
    T.tau_w = nan(n,1);
    T.cf    = nan(n,1);
    T.rho_w = nan(n,1);
    T.mu_w  = nan(n,1);
    T.nu_w  = nan(n,1);

    % Riblet sizing in wall units (unitless). If you fill these from validation/correlation,
    % the main loop will use them instead of the default correlation.
    T.sPlus = nan(n,1);
    T.hPlus = nan(n,1);

    writetable(T, filename);
    if nargin < 3 || ~isfield(cfg,'printInfoLines') || cfg.printInfoLines
        fprintf('[INFO] Wrote STL-aligned wall template: %s\n', filename);
    end
end

function files = buildWallFileList(prefix, machCases)
% Builds filenames like wall_M080.csv, wall_M085.csv, ... wall_M140.csv

    files = strings(size(machCases));
    for k = 1:numel(machCases)
        tag = sprintf('M%03d', round(machCases(k)*100));
        files(k) = prefix + tag + '.csv';
    end
end

%% ===================== PRINTING =====================
function printRunConditions(cfg, prof, mid)
    lines = strings(0);
    lines(end+1) = "================ RUN CONDITIONS ================";
    lines(end+1) = "Object: " + cfg.objectName;
    lines(end+1) = "Timestamp: " + string(datetime("now"));
    lines(end+1) = "";

    lines(end+1) = "---- STL / geometry ----";
    lines(end+1) = "STL file: " + string(cfg.stlFile);
    lines(end+1) = "Scale factor: " + num2str(cfg.scaleFactor);
    lines(end+1) = sprintf("Estimated length (principal axis): %.3f mm", prof.L);
    lines(end+1) = "nbins: " + num2str(cfg.nbins);
    lines(end+1) = "outlierHiPrct: " + num2str(cfg.outlierHiPrct) + " %";
    lines(end+1) = "";

    lines(end+1) = "---- Mid-body selection ----";
    lines(end+1) = sprintf("xFracRange: [%.2f, %.2f]", cfg.xFracRange(1), cfg.xFracRange(2));
    lines(end+1) = sprintf("radiusTolFrac: %.3f", cfg.radiusTolFrac);
    lines(end+1) = sprintf("Mid-body bins: %d", numel(mid.binIdx));
    lines(end+1) = sprintf("Mid-body x range: [%.3f, %.3f] mm", min(mid.x), max(mid.x));
    lines(end+1) = sprintf("Mid-body R median: %.3f mm", median(mid.R,'omitnan'));
    lines(end+1) = "";

    lines(end+1) = "---- Wall cases ----";
    lines(end+1) = "wallCoordMode: " + string(cfg.wallCoordMode);
    lines(end+1) = "createWallTemplates: " + string(cfg.createWallTemplates);
    lines(end+1) = "overwriteTemplates: " + string(cfg.overwriteTemplates);
    for k = 1:numel(cfg.machCases)
        lines(end+1) = sprintf("Mach %.2f -> %s", cfg.machCases(k), string(cfg.wallFiles(k)));
    end
    lines(end+1) = "================================================";
    if isfield(cfg,'printRunConditionsToConsole') && cfg.printRunConditionsToConsole
        fprintf("\n%s\n\n\", strjoin(lines, newline));
    end
    if cfg.saveRunReport
        fid = fopen(cfg.reportFile, 'w');
        if fid < 0
            warning("Could not write report file: %s", cfg.reportFile);
            return;
        end
        fwrite(fid, strjoin(lines, newline));
        fclose(fid);
        if ~isfield(cfg,'printInfoLines') || cfg.printInfoLines
            fprintf("[INFO] Wrote run report: %s\n\", cfg.reportFile);
            if cfg.exportProfileCSV
                fprintf("[INFO] Wrote mid-body geometry CSV: %s\n\", cfg.profileCSVName);
            end
        end
    end
end

%% ===================== PLOTTING =====================
function plotMachFigures(cfg, M, x_mm, p_use, tau_use, sPlus, hPlus)
% One figure per Mach increment (x in mm).

    tag = string(sprintf('%s_M%03d', char(cfg.objectName), round(M*100)));

    % --- Riblet wall-unit sizing (unitless) ---
    figure('Name', tag + " riblets", 'NumberTitle','off');
    grid on; hold on;
    plot(x_mm, sPlus, '-', 'LineWidth', 1.6);
    plot(x_mm, hPlus, '-', 'LineWidth', 1.6);
    xlabel('x (mm)');
    ylabel('Riblet wall units (unitless)');
    titleBelowXLabel(gca, cfg.objectName + sprintf(': Riblet wall units (Mach %.2f)', M));
    legend({'s^+ (unitless)','h^+ (unitless)'}, 'Location','best');
    if cfg.saveFigs
        exportgraphics(gcf, fullfile(cfg.figDir, tag + "_riblets.png"));
    end

    % --- Wall shear stress (separate plot; do NOT combine with riblets) ---
    if any(~isnan(tau_use))
        figure('Name', tag + " tau_w", 'NumberTitle','off');
        grid on;
        plot(x_mm, tau_use, '-', 'LineWidth', 1.6);
        xlabel('x (mm)');
        ylabel('\tau_w (Pa)');
        titleBelowXLabel(gca, cfg.objectName + sprintf(': \tau_w vs x (Mach %.2f)', M));
        if cfg.saveFigs
            exportgraphics(gcf, fullfile(cfg.figDir, tag + "_tau_w.png"));
        end
    end
end
function plotOverlays(cfg, results)
    Ms = [results.Mach];

    if any(arrayfun(@(r) any(~isnan(r.tau_w)), results))
        figure('Name', cfg.objectName + " tau_w overlay", 'NumberTitle','off'); grid on; hold on;
        for k = 1:numel(results)
            if any(~isnan(results(k).tau_w))
                plot(results(k).x, results(k).tau_w, 'LineWidth', 1.2);
            end
        end
        xlabel('x (mm)'); ylabel('\tau_w (Pa)');
        titleBelowXLabel(gca, cfg.objectName + ": \tau_w overlays");
        legend(arrayfun(@(m) sprintf("M=%.2f", m), Ms, 'UniformOutput', false), 'Location','best');
        if cfg.saveFigs, exportgraphics(gcf, fullfile(cfg.figDir, cfg.objectName + "_overlay_tau_w.png")); end
    end
end
function plotVsMach(cfg, results)
% Creates end-of-run comparison plots vs Mach for the variables you requested
% (p, q, tau_w, s^+, h^+). Uses medians over the STL mid-body region.

    Ms = [results.Mach];
    pM   = [results.p_med];
    qM   = [results.q_med];
    tauM = [results.tau_med];
    sM   = [results.sPlus_med];
    hM   = [results.hPlus_med];

    tag = cfg.objectName + "_vsMach";

    figure('Name', cfg.objectName + " summary vs Mach", 'NumberTitle','off');

    % --- Pressures ---
    subplot(3,1,1); grid on; hold on;
    haveP = any(~isnan(pM));
    haveQ = any(~isnan(qM));
    if haveP, plot(Ms, pM, '-o', 'LineWidth', 1.4); end
    if haveQ, plot(Ms, qM, '-o', 'LineWidth', 1.4); end
    xlabel('Mach'); ylabel('Pressure (Pa)');
    titleBelowXLabel(gca, cfg.objectName + " : p and q vs Mach (median)");
    if haveP && haveQ
        legend({'p_w (Pa)','q_{\infty} (Pa)'}, 'Location','best');
    elseif haveP
        legend({'p_w (Pa)'}, 'Location','best');
    elseif haveQ
        legend({'q_{\infty} (Pa)'}, 'Location','best');
    end

    % --- Wall shear stress ---
    subplot(3,1,2); grid on;
    plot(Ms, tauM, '-o', 'LineWidth', 1.4);
    xlabel('Mach'); ylabel('\tau_w (Pa)');
    titleBelowXLabel(gca, cfg.objectName + " : \tau_w vs Mach (median)");
    % --- Riblet wall units ---
    subplot(3,1,3); grid on; hold on;
    plot(Ms, sM, '-o', 'LineWidth', 1.4);
    plot(Ms, hM, '-o', 'LineWidth', 1.4);
    xlabel('Mach'); ylabel('Wall units (unitless)');
    titleBelowXLabel(gca, cfg.objectName + " : s^+ and h^+ vs Mach (median)");
    legend({'s^+','h^+'}, 'Location','best');

    if cfg.saveFigs
        exportgraphics(gcf, fullfile(cfg.figDir, tag + "_summary.png"));
    end
end


%% ===================== WALL CSV HELPERS =====================
function wall = readWallCSV(filename, wallScaleFactor)
    if nargin < 2, wallScaleFactor = 1.0; end

    T = readtable(filename);
    wall = struct();

    % Axial coordinate (expected mm). Apply scale if your export is meters.
    wall.x = getCol(T, ["x_mm","x","Xaxial","x_axial","x_m","x (m)"], true);
    if ~isempty(wall.x)
        wall.x = wall.x .* wallScaleFactor;
    end

    wall.X = getCol(T, ["X","xCoord","x_coord"], true);
    wall.Y = getCol(T, ["Y","yCoord","y_coord"], true);
    wall.Z = getCol(T, ["Z","zCoord","z_coord"], true);

    wall.p_w   = getCol(T, ["p_w","p","pressure","p (Pa)"], true);
    wall.q_inf = getCol(T, ["q_inf","q","dynPress","q (Pa)"], true);
    wall.tau_w = getCol(T, ["tau_w","tauw","wallShear","wall_shear"], true);
    wall.cf    = getCol(T, ["cf","Cf","skinFriction","skin_friction"], true);

    wall.rho_w = getCol(T, ["rho_w","rhow","rho","density"], true);
    wall.mu_w  = getCol(T, ["mu_w","muw","mu","viscosity"], true);
    wall.nu_w  = getCol(T, ["nu_w","nuw","nu","kinVisc"], true);

    wall.sPlus = getCol(T, ["sPlus","s_plus","s+"], true);
    wall.hPlus = getCol(T, ["hPlus","h_plus","h+"], true);
end

function x_axial = computeWallAxialX(wall, mode, prof)
    mode = lower(string(mode));
    switch mode
        case "axial"
            if ~isempty(wall.x)
                x_axial = wall.x(:);
            elseif ~isempty(wall.X) && ~isempty(wall.Y) && ~isempty(wall.Z)
                P = [wall.X(:) wall.Y(:) wall.Z(:)];
                x_axial = projectToAxis(P, prof.centerOriginal, prof.exOriginal);
            else
                error("wallCoordMode='axial' but no axial x or XYZ columns found.");
            end
        case "xyz_project"
            if isempty(wall.X) || isempty(wall.Y) || isempty(wall.Z)
                error("wallCoordMode='xyz_project' requires XYZ columns (X,Y,Z).");
            end
            P = [wall.X(:) wall.Y(:) wall.Z(:)];
            x_axial = projectToAxis(P, prof.centerOriginal, prof.exOriginal);
        otherwise
            error("Unknown wallCoordMode. Use 'axial' or 'xyz_project'.");
    end
end

function xproj = projectToAxis(P, C, ex)
    ex = ex(:) ./ max(norm(ex), eps);
    PC = P - C(:).';
    xproj = PC * ex;
end

function agg = aggregateWallToBins(wall, prof, mid, method)
    method = lower(string(method));
    edges = prof.edges;

    agg = struct();
    agg.x = prof.xc(mid.binIdx);

    fields = ["p_w","q_inf","tau_w","cf","rho_w","mu_w","nu_w","sPlus","hPlus"];
    for f = fields
        val = wall.(f);
        out = nan(size(agg.x));
        if isempty(val)
            agg.(f) = out;
            continue;
        end

        for k = 1:numel(mid.binIdx)
            b = mid.binIdx(k);
            xlo = edges(b);
            xhi = edges(b+1);
            sel = (wall.x_axial >= xlo) & (wall.x_axial < xhi);

            if ~any(sel), continue; end
            vk = double(val(sel));

            switch method
                case "median"
                    out(k) = median(vk, 'omitnan');
                case "mean"
                    out(k) = mean(vk, 'omitnan');
                otherwise
                    error("Unknown aggregation method. Use 'median' or 'mean'.");
            end
        end
        agg.(f) = out(:);
    end
end

function d = computeDerivedFields(agg)
    d = struct();

    d.tau_w = agg.tau_w;
    if all(isnan(d.tau_w)) && any(~isnan(agg.cf)) && any(~isnan(agg.q_inf))
        d.tau_w = 0.5 .* agg.cf .* agg.q_inf;
    end
    d.cf = agg.cf;

    d.nu_w = agg.nu_w;
    if all(isnan(d.nu_w)) && any(~isnan(agg.mu_w)) && any(~isnan(agg.rho_w))
        d.nu_w = agg.mu_w ./ max(agg.rho_w, eps);
    end

    d.u_tau = nan(size(d.tau_w));
    if any(~isnan(d.tau_w)) && any(~isnan(agg.rho_w))
        d.u_tau = sqrt(max(d.tau_w,0) ./ max(agg.rho_w, eps));
    end
end

function col = getCol(T, names, allowMissing)
    if nargin < 3, allowMissing = false; end
    col = [];
    vnames = T.Properties.VariableNames;

    for i = 1:numel(names)
        nm = string(names(i));
        hit = find(strcmpi(vnames, nm), 1);
        if ~isempty(hit)
            col = double(T.(vnames{hit}));
            return;
        end
    end

    if allowMissing
        col = [];
    else
        error("Missing required column. Tried: %s", strjoin(cellstr(names), ", "));
    end
end


%% ===================== THEORY / RIBLET HELPERS =====================
function fs = freestreamSeaLevel(M)
% Simple ISA sea-level freestream (replace with your chapter assumptions if needed).
    gamma = 1.4;
    R = 287.05287;
    T = 288.15;
    p = 101325;
    rho = p/(R*T);
    a = sqrt(gamma*R*T);
    V = M*a;
    mu = sutherlandMu(T);
    nu = mu/rho;
    q = 0.5*rho*V^2;
    fs = struct('gamma',gamma,'R',R,'T',T,'p',p,'rho',rho,'a',a,'V',V,'mu',mu,'nu',nu,'q',q);
end

function mu = sutherlandMu(T)
% Sutherland's law for air dynamic viscosity (Pa*s).
    mu0 = 1.716e-5;
    T0 = 273.15;
    S  = 110.4;
    mu = mu0*(T/T0)^(3/2)*(T0+S)/(T+S);
end

function [p_use, q_use, tau_use] = fillPQTauProfiles(agg, d, fs, M)
% Fill pressure/dynamic pressure/tau_w profiles from CSV when present; otherwise use theory defaults.

    if nargin < 4, M = NaN; end

    x_mm = agg.x(:);

    % pressure
    p_use = agg.p_w;
    if all(isnan(p_use))
        p_use = fs.p .* ones(size(x_mm));
    end

    % dynamic pressure
    q_use = agg.q_inf;
    if all(isnan(q_use))
        q_use = fs.q .* ones(size(x_mm));
    end

    % wall shear
    tau_use = d.tau_w;
    if all(isnan(tau_use))
        % Turbulent local flat-plate estimate with mild compressibility scaling
        x_m = max(x_mm, 1e-3) ./ 1000;
        Re_x = fs.rho .* fs.V .* x_m ./ max(fs.mu, eps);
        cf0 = 0.0592 ./ max(Re_x, 1).^0.2;
        if ~isnan(M)
            comp = (1 + 0.2*M^2)^(-0.467);
        else
            comp = 1.0;
        end
        cf  = cf0 .* comp;
        tau_use = 0.5 .* cf .* q_use;
    end
end

function rib = validatedRibletsUnitless(M, x_mm, p_profile, cfg)
% Returns unitless riblet spacing/height in wall units (s^+, h^+).
% Includes a simple Mach>=1 shock/BL interaction proxy based on |dp/dx|.

    s0 = cfg.riblet_sPlus0;
    h0 = cfg.riblet_hPlus0;

    sPlus = s0 .* ones(size(x_mm));
    hPlus = h0 .* ones(size(x_mm));

    % Mild compressibility scaling (kept small)
    compFactor = 1 + cfg.riblet_subsonicK * M^2;
    sPlus = sPlus ./ compFactor;
    hPlus = hPlus ./ compFactor;

    % Mach>=1: shock/BL interaction proxy (applied in code only; do not print theory text)
    if M >= 1 && any(~isnan(p_profile))
        x_m = max(x_mm, 1e-3) ./ 1000;
        dpdx = gradient(p_profile(:), x_m);
        denom = max(median(abs(dpdx), 'omitnan'), eps);
        shockIndex = abs(dpdx) ./ denom;
        shockIndex = min(shockIndex, cfg.riblet_maxShockIndex);
        shockFactor = 1 + cfg.riblet_shockK * (M - 1) .* shockIndex;
        sPlus = sPlus ./ shockFactor;
        hPlus = hPlus ./ shockFactor;
    end

    rib = struct('sPlus', sPlus, 'hPlus', hPlus);
end

function v = medOrNaN(x)
    if isempty(x)
        v = NaN; return;
    end
    if all(isnan(x))
        v = NaN;
    else
        v = median(x, 'omitnan');
    end
end

%% ===================== STL HELPERS =====================
function [V,F] = readSTL_any(filename)
    if exist('stlread','file') == 2
        try
            TR = stlread(filename);
            V = TR.Points;
            F = TR.ConnectivityList;
            return;
        catch
        end
    end
    try
        [V,F] = readSTL_binary(filename);
    catch
        [V,F] = readSTL_ascii(filename);
    end
end

function [V,F] = readSTL_binary(filename)
    fid = fopen(filename,'rb');
    if fid < 0, error("Cannot open STL file: %s", filename); end
    cleaner = onCleanup(@() fclose(fid));

    fseek(fid, 80, 'bof');
    nTri = fread(fid, 1, 'uint32');
    if isempty(nTri) || nTri==0
        error("STL empty or not binary: %s", filename);
    end

    V = zeros(nTri*3, 3);
    F = reshape(1:(nTri*3), 3, [])';

    for i = 1:nTri
        fread(fid, 3, 'float32'); % normal
        v1 = fread(fid, 3, 'float32')';
        v2 = fread(fid, 3, 'float32')';
        v3 = fread(fid, 3, 'float32')';
        fread(fid, 1, 'uint16');  % attribute
        base = (i-1)*3;
        V(base+1,:) = v1;
        V(base+2,:) = v2;
        V(base+3,:) = v3;
    end
end

function [V,F] = readSTL_ascii(filename)
    txt = fileread(filename);
    tok = regexp(txt, 'vertex\s+([-\d\.eE+]+)\s+([-\d\.eE+]+)\s+([-\d\.eE+]+)', 'tokens');
    if isempty(tok)
        error("Could not parse STL as ASCII: %s", filename);
    end
    V = cellfun(@(c) [str2double(c{1}), str2double(c{2}), str2double(c{3})], tok, 'UniformOutput', false);
    V = vertcat(V{:});
    if mod(size(V,1),3) ~= 0
        error("ASCII STL parse error (vertex count not multiple of 3): %s", filename);
    end
    nTri = size(V,1)/3;
    F = reshape(1:(nTri*3), 3, [])';
end

function prof = extractBodyProfileFromSTL(V, nbins, outlierHiPrct)
    C = mean(V,1);
    X = V - C;

    [~,~,Vsvd] = svd(X, 'econ');
    ex = Vsvd(:,1);
    ey = Vsvd(:,2);
    ez = Vsvd(:,3);

    T = [ex ey ez];
    Xb = (T' * X')';
    x = Xb(:,1); y = Xb(:,2); z = Xb(:,3);

    L = max(x) - min(x);
    if L <= 0, error("Degenerate STL: zero axial length."); end

    xA = min(x) + 0.10*L;
    xB = min(x) + 0.90*L;
    rA = robustRadiusNear(x,y,z,xA,outlierHiPrct);
    rB = robustRadiusNear(x,y,z,xB,outlierHiPrct);

    if rA > rB
        ex = -ex;
        T  = [-T(:,1) T(:,2) T(:,3)];
        Xb = (T' * X')';
        x = Xb(:,1); y = Xb(:,2); z = Xb(:,3);
    end

    edges = linspace(min(x), max(x), nbins+1);
    xc    = 0.5*(edges(1:end-1) + edges(2:end));
    Rb    = nan(size(xc));

    for i = 1:numel(xc)
        inBin = (x >= edges(i)) & (x < edges(i+1));
        if nnz(inBin) < 30, continue; end

        yi = y(inBin); zi = z(inBin);
        cy = median(yi); cz = median(zi);

        ri = sqrt((yi-cy).^2 + (zi-cz).^2);
        cutoff = prctile(ri, outlierHiPrct);
        ri = ri(ri <= cutoff);

        if ~isempty(ri)
            Rb(i) = median(ri);
        end
    end

    good = ~isnan(Rb);
    if nnz(good) < 5
        error("Radius profile extraction failed. Adjust nbins/outlierHiPrct.");
    end
    Rb = interp1(xc(good), Rb(good), xc, 'linear', 'extrap');

    prof = struct();
    prof.xc = xc(:);
    prof.Rb = Rb(:);
    prof.edges = edges(:);
    prof.L = max(x) - min(x);
    prof.centerOriginal = C(:);
    prof.exOriginal = ex(:);
end

function R = robustRadiusNear(x,y,z,x0,outlierHiPrct)
    L = max(x)-min(x);
    dx = 0.02*L;
    sel = abs(x - x0) <= dx;
    if nnz(sel) < 30, R = 0; return; end
    cy = median(y(sel)); cz = median(z(sel));
    r  = sqrt((y(sel)-cy).^2 + (z(sel)-cz).^2);
    r  = r(r <= prctile(r, outlierHiPrct));
    if isempty(r), R = 0; else, R = median(r); end
end

function mid = selectMidBody(prof, xFracRange, radiusTolFrac)
    x = prof.xc; R = prof.Rb;
    xmin = min(x); xmax = max(x); L = xmax - xmin;

    xlo = xmin + xFracRange(1)*L;
    xhi = xmin + xFracRange(2)*L;

    inX = (x >= xlo) & (x <= xhi);
    if ~any(inX)
        mid = struct('binIdx',[],'x',[],'R',[]);
        return;
    end

    Rmed = median(R(inX));
    inR  = abs(R - Rmed) <= radiusTolFrac*Rmed;
    keep = inX & inR;

    mid = struct();
    mid.binIdx = find(keep);
    mid.x = x(keep);
    mid.R = R(keep);
end

function titleBelowXLabel(ax, str)
% bottomTitle Place a plot title at the bottom INSIDE the axes.
%   This is useful when you want a caption-like title below the data rather
%   than the default top title. The string is placed in normalized axes
%   coordinates at (0.5, 0.03).

    if nargin < 1 || isempty(ax); ax = gca; end
    if nargin < 2; str = ""; end

    % Clear standard top title
    title(ax, "");

    % Remove prior bottom-title text (if any)
    old = findall(ax, 'Tag', 'BottomTitleText');
    if ~isempty(old); delete(old); end

    % Add bottom title inside the plotting area
    text(ax, 0.5, 0.03, char(str), ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontWeight','bold', ...
        'Interpreter','tex', ...
        'Tag','BottomTitleText');
end