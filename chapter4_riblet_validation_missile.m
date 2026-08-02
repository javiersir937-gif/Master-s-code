%% chapter4_riblet_validation_missile.m
% Chapter 4: STL-based mid-body extraction + per-Mach wall-data plotting
% Adds:
%  - Prints run parameters to Command Window (and saves a .txt report)
%  - Creates STL-aligned wall CSV templates for each Mach case (x, R, circumference + NaNs)
%  - Generates separate figures per Mach (tau_w, Cf, nu_w, u_tau, sPlus if present)

clear; clc; close all;

%% ===================== USER SETTINGS =====================
cfg = struct();

cfg.objectName  = "Missile";

% STL (match your uploaded filename/case)
cfg.stlFile     = "testm.STL";
cfg.scaleFactor = 1.0;       % if STL units are mm, set 1e-3

% Radius profile extraction
cfg.nbins         = 250;
cfg.outlierHiPrct = 90;

% Mid-body selection
cfg.xFracRange    = [0.20 0.80];
cfg.radiusTolFrac = 0.10;

% Mach cases + wall CSV names
cfg.machCases = [0.8 1.0 1.2 1.4];
cfg.wallFiles = ["wall_M08.csv","wall_M10.csv","wall_M12.csv","wall_M14.csv"];

% Wall coordinate mode:
%  "axial" expects a column x (already axial)
%  "xyz_project" expects X,Y,Z columns and projects onto STL axis
cfg.wallCoordMode = "axial";

% Outputs
cfg.exportProfileCSV = true;
cfg.profileCSVName   = "missile_midbody_profile.csv";

cfg.saveRunReport = true;
cfg.reportFile    = "run_conditions_missile.txt";

cfg.saveFigs = true;
cfg.figDir   = "figs_missile";

cfg.makeCombinedOverlays = true;

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
            'VariableNames', {'x_axial_m','R_m','circumference_m'});
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
    xlabel('x along inferred principal axis (m)');
    ylabel('robust radius R(x) (m)');
    legend('Robust radius profile','Mid-body bins','Location','best');
    title(cfg.objectName + ": STL-derived robust radius profile");
    if cfg.saveFigs
        exportgraphics(gcf, fullfile(cfg.figDir, cfg.objectName + "_stl_radius_profile.png"));
    end

    % ---- STEP 2: Print run parameters ----
    printRunConditions(cfg, prof, mid);

    % ---- STEP 3: Per Mach case plots ----
    results = struct([]);

    for i = 1:numel(cfg.machCases)
        M = cfg.machCases(i);
        f = string(cfg.wallFiles(i));

        if strlength(f)==0
            fprintf("[INFO] No wall file specified for Mach %.2f (skipping)\n", M);
            continue;
        end
        if ~isfile(f)
            fprintf("[WARN] Missing wall file for Mach %.2f: %s (skipping)\n", M, f);
            continue;
        end

        wall = readWallCSV(f);

        % Ensure an axial coordinate exists (x_axial)
        wall.x_axial = computeWallAxialX(wall, cfg.wallCoordMode, prof);

        % Aggregate wall data into STL mid-body bins
        agg = aggregateWallToBins(wall, prof, mid, "median");

        % Derived fields (tau_w, nu_w, u_tau) + optional sPlus passthrough
        d = computeDerivedFields(agg);
        d.sPlus = agg.sPlus;

        % Print per-case summary
        fprintf("\n[%s | Mach %.2f]\n", cfg.objectName, M);
        fprintf("  Wall file: %s\n", f);
        fprintf("  Mid-body x range: [%.4f, %.4f] m\n", min(agg.x), max(agg.x));

        if any(~isnan(d.tau_w))
            fprintf("  tau_w (Pa) min/median/max: %.3g / %.3g / %.3g\n", ...
                min(d.tau_w,[],'omitnan'), median(d.tau_w,'omitnan'), max(d.tau_w,[],'omitnan'));
        end
        if any(~isnan(d.u_tau))
            fprintf("  u_tau (m/s) min/median/max: %.3g / %.3g / %.3g\n", ...
                min(d.u_tau,[],'omitnan'), median(d.u_tau,'omitnan'), max(d.u_tau,[],'omitnan'));
        end
        if any(~isnan(d.nu_w))
            fprintf("  nu_w (m^2/s) min/median/max: %.3g / %.3g / %.3g\n", ...
                min(d.nu_w,[],'omitnan'), median(d.nu_w,'omitnan'), max(d.nu_w,[],'omitnan'));
        end
        if any(~isnan(d.cf))
            fprintf("  Cf min/median/max: %.3g / %.3g / %.3g\n", ...
                min(d.cf,[],'omitnan'), median(d.cf,'omitnan'), max(d.cf,[],'omitnan'));
        end
        if any(~isnan(d.sPlus))
            fprintf("  s+ min/median/max: %.3g / %.3g / %.3g\n", ...
                min(d.sPlus,[],'omitnan'), median(d.sPlus,'omitnan'), max(d.sPlus,[],'omitnan'));
        end

        % Separate figures per Mach
        plotMachFigures(cfg, M, agg.x, d);

        % Store for overlays
        results(end+1).Mach  = M; %#ok<AGROW>
        results(end).x       = agg.x;
        results(end).tau_w   = d.tau_w;
        results(end).cf      = d.cf;
        results(end).nu_w    = d.nu_w;
        results(end).u_tau   = d.u_tau;
        results(end).sPlus   = d.sPlus;
    end

    if cfg.makeCombinedOverlays && ~isempty(results)
        plotOverlays(cfg, results);
    end
end

%% ===================== NEW: WALL TEMPLATE WRITER =====================
function writeWallTemplatesFromSTL(mid, cfg)
% Creates/overwrites wall CSV files using STL-derived mid-body x/R/circumference.
% Flow columns are NaN placeholders meant to be replaced by Fluent/CFD exports.

    for i = 1:numel(cfg.machCases)
        f = string(cfg.wallFiles(i));
        if strlength(f)==0
            continue;
        end

        if isfile(f) && ~cfg.overwriteTemplates
            fprintf("[INFO] Wall file exists (not overwritten): %s\n", f);
            continue;
        end

        T = table();
        T.x = mid.x(:);
        T.R = mid.R(:);
        T.circumference = 2*pi*mid.R(:);

        n = height(T);

        % placeholders for CFD values
        T.tau_w = nan(n,1);
        T.rho_w = nan(n,1);
        T.mu_w  = nan(n,1);
        T.nu_w  = nan(n,1);
        T.cf    = nan(n,1);
        T.q_inf = nan(n,1);
        T.sPlus = nan(n,1);

        writetable(T, f);
        fprintf("[INFO] Wrote STL-aligned wall template: %s\n", f);
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
    lines(end+1) = sprintf("Estimated length (principal axis): %.6f m", prof.L);
    lines(end+1) = "nbins: " + num2str(cfg.nbins);
    lines(end+1) = "outlierHiPrct: " + num2str(cfg.outlierHiPrct) + " %";
    lines(end+1) = "";

    lines(end+1) = "---- Mid-body selection ----";
    lines(end+1) = sprintf("xFracRange: [%.2f, %.2f]", cfg.xFracRange(1), cfg.xFracRange(2));
    lines(end+1) = sprintf("radiusTolFrac: %.3f", cfg.radiusTolFrac);
    lines(end+1) = sprintf("Mid-body bins: %d", numel(mid.binIdx));
    lines(end+1) = sprintf("Mid-body x range: [%.6f, %.6f] m", min(mid.x), max(mid.x));
    lines(end+1) = sprintf("Mid-body R median: %.6f m", median(mid.R,'omitnan'));
    lines(end+1) = "";

    lines(end+1) = "---- Wall cases ----";
    lines(end+1) = "wallCoordMode: " + string(cfg.wallCoordMode);
    lines(end+1) = "createWallTemplates: " + string(cfg.createWallTemplates);
    lines(end+1) = "overwriteTemplates: " + string(cfg.overwriteTemplates);
    for k = 1:numel(cfg.machCases)
        lines(end+1) = sprintf("Mach %.2f -> %s", cfg.machCases(k), string(cfg.wallFiles(k)));
    end
    lines(end+1) = "================================================";

    fprintf("\n%s\n\n", strjoin(lines, newline));

    if cfg.saveRunReport
        fid = fopen(cfg.reportFile, 'w');
        if fid < 0
            warning("Could not write report file: %s", cfg.reportFile);
            return;
        end
        fwrite(fid, strjoin(lines, newline));
        fclose(fid);
        fprintf("[INFO] Wrote run report: %s\n", cfg.reportFile);
        if cfg.exportProfileCSV
            fprintf("[INFO] Wrote mid-body geometry CSV: %s\n", cfg.profileCSVName);
        end
    end
end

%% ===================== PLOTTING =====================
function plotMachFigures(cfg, M, x, d)
    tag = sprintf("%s_M%.2f", cfg.objectName, M);

    if any(~isnan(d.tau_w))
        figure('Name', tag + " tau_w", 'NumberTitle','off'); grid on;
        plot(x, d.tau_w, '-', 'LineWidth', 1.6);
        xlabel('x (m)'); ylabel('\tau_w (Pa)');
        title(cfg.objectName + ": \tau_w vs x (Mach " + num2str(M) + ")");
        if cfg.saveFigs, exportgraphics(gcf, fullfile(cfg.figDir, tag + "_tau_w.png")); end
    end

    if any(~isnan(d.cf))
        figure('Name', tag + " Cf", 'NumberTitle','off'); grid on;
        plot(x, d.cf, '-', 'LineWidth', 1.6);
        xlabel('x (m)'); ylabel('C_f');
        title(cfg.objectName + ": C_f vs x (Mach " + num2str(M) + ")");
        if cfg.saveFigs, exportgraphics(gcf, fullfile(cfg.figDir, tag + "_Cf.png")); end
    end

    if any(~isnan(d.nu_w))
        figure('Name', tag + " nu_w", 'NumberTitle','off'); grid on;
        plot(x, d.nu_w, '-', 'LineWidth', 1.6);
        xlabel('x (m)'); ylabel('\nu_w (m^2/s)');
        title(cfg.objectName + ": \nu_w vs x (Mach " + num2str(M) + ")");
        if cfg.saveFigs, exportgraphics(gcf, fullfile(cfg.figDir, tag + "_nu_w.png")); end
    end

    if any(~isnan(d.u_tau))
        figure('Name', tag + " u_tau", 'NumberTitle','off'); grid on;
        plot(x, d.u_tau, '-', 'LineWidth', 1.6);
        xlabel('x (m)'); ylabel('u_\tau (m/s)');
        title(cfg.objectName + ": u_\tau vs x (Mach " + num2str(M) + ")");
        if cfg.saveFigs, exportgraphics(gcf, fullfile(cfg.figDir, tag + "_u_tau.png")); end
    end

    if any(~isnan(d.sPlus))
        figure('Name', tag + " sPlus", 'NumberTitle','off'); grid on;
        plot(x, d.sPlus, '-', 'LineWidth', 1.6);
        xlabel('x (m)'); ylabel('s^+');
        title(cfg.objectName + ": s^+ vs x (Mach " + num2str(M) + ")");
        if cfg.saveFigs, exportgraphics(gcf, fullfile(cfg.figDir, tag + "_sPlus.png")); end
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
        xlabel('x (m)'); ylabel('\tau_w (Pa)');
        title(cfg.objectName + ": \tau_w overlays");
        legend(arrayfun(@(m) sprintf("M=%.2f", m), Ms, 'UniformOutput', false), 'Location','best');
        if cfg.saveFigs, exportgraphics(gcf, fullfile(cfg.figDir, cfg.objectName + "_overlay_tau_w.png")); end
    end
end

%% ===================== WALL CSV HELPERS =====================
function wall = readWallCSV(filename)
    T = readtable(filename);
    wall = struct();

    wall.x = getCol(T, ["x","Xaxial","x_axial","x_m","x (m)"], true);

    wall.X = getCol(T, ["X","xCoord","x_coord"], true);
    wall.Y = getCol(T, ["Y","yCoord","y_coord"], true);
    wall.Z = getCol(T, ["Z","zCoord","z_coord"], true);

    wall.tau_w = getCol(T, ["tau_w","tauw","wallShear","wall_shear"], true);
    wall.cf    = getCol(T, ["cf","Cf","skinFriction","skin_friction"], true);
    wall.q_inf = getCol(T, ["q_inf","q","dynPress","q (Pa)"], true);

    wall.rho_w = getCol(T, ["rho_w","rhow","rho","density"], true);
    wall.mu_w  = getCol(T, ["mu_w","muw","mu","viscosity"], true);
    wall.nu_w  = getCol(T, ["nu_w","nuw","nu","kinVisc"], true);

    wall.sPlus = getCol(T, ["sPlus","s_plus","s+"], true);
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

    fields = ["tau_w","cf","q_inf","rho_w","mu_w","nu_w","sPlus"];
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
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

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