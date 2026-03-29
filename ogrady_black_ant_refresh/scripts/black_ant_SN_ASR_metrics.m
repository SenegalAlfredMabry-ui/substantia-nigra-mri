% SN_ASR_regenerate_SN_metrics.m
% -------------------------------------------------------------------------
% PURPOSE
%   Recompute SN neuromelanin intensity + CNR metrics for every participant
%   referenced by SN_ASR_DRAFT_analysis.Rmd.
%   Uses the canonical pipeline block from twin_analysis_tse_scale_draft.m
%   (CNR from raw NM, scaled intensities from corrected NM) but operates over
%   the full SN cohort and writes the refreshed table to the SNHEART tree.
%
% OUTPUT
%   ${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data/analysis/SN_nmdata_autosaved.txt
%   (comma-delimited, same column order expected by SN_ASR.)
%
% REQUIREMENTS
%   - Access to the legacy MOTIP data tree via ${MOTIP_ROOT}/... or the
%     /export${MOTIP_ROOT}/... mirror (script auto-falls back if /Volumes is
%     not mounted).
%   - NIfTI_20140122 toolbox bundled with the legacy scripts.
% -------------------------------------------------------------------------

clear; close all;

% ---------------------- Locate source + toolbox paths ----------------------
legacy_base_primary = '${MOTIP_ROOT:-/path/to/MOTIP2018}';
legacy_base_fallback = '${MOTIP_ROOT:-/path/to/MOTIP2018}';

% Ensure we can reach the historical /Volumes paths (fallback to /export when
% needed so the script still runs on systems without the macOS-style mount).
legacy_tse_root = resolve_existing_path( ...
    fullfile(legacy_base_primary, 'MRI_data', 'TSE'), ...
    fullfile(legacy_base_fallback, 'MRI_data', 'TSE'), ...
    'TSE root');
legacy_scripts_root = resolve_existing_path( ...
    fullfile(legacy_base_primary, 'scripts', 'tse'), ...
    fullfile(legacy_base_fallback, 'scripts', 'tse'), ...
    'TSE scripts');

localpath = convertCharsToStrings(legacy_tse_root);
localpath_char = char(localpath);  % convenient char version for sprintf
addpath(fullfile(legacy_scripts_root, 'NIfTI_20140122'));

% ---------------------- Derive participant list ----------------------------
subject_table_source = fullfile(legacy_scripts_root, 'SN_nmdata_autosaved.txt');
if ~isfile(subject_table_source)
    error('Cannot find subject list seed: %s', subject_table_source);
end
opts = detectImportOptions(subject_table_source, 'Delimiter', ',');
opts = setvartype(opts, 'sub', 'double');
sn_table = readtable(subject_table_source, opts);
subjects = sn_table.sub;
subjects = subjects(~isnan(subjects));
subjects = sort(unique(subjects));

fprintf('Found %d SN participants in %s\n', numel(subjects), subject_table_source);

% ---------------------- Constants and allocation ---------------------------
voxelvol = 0.43 * 0.43 * 3;  % mm^3
results = struct('sub',{}, 'CNR',{}, 'total_intensity',{}, 'intensity',{}, ...
                 'background_intensity',{}, 'total_volume',{}, ...
                 'intensity_minimum',{}, 'intensity_maximum',{}, ...
                 'mean_raw_sn',{}, 'mean_raw_control',{}, 'std_raw_control',{});

% ---------------------- Loop over subjects ---------------------------------
for idx = 1:numel(subjects)
    sub = subjects(idx);
    fprintf('Processing subject %d (%d/%d)\n', sub, idx, numel(subjects));

    try
        pth = build_paths(localpath_char, sub);

        % Load corrected neuromelanin volume (with optional session suffixes).
        corrected_candidates = { ...
            sprintf('%s/SN/corrected_tse_DC_%d.nii', localpath_char, sub), ...
            sprintf('%s/SN/corrected_tse_DC_%d_pre.nii', localpath_char, sub), ...
            sprintf('%s/SN/corrected_tse_DC_%d_post.nii', localpath_char, sub)};
        subjnm = load_first_available(corrected_candidates, 'corrected TSE');
        nm = double(subjnm.img);
        nm(nm == -999) = NaN;

        % Raw TSE for CNR computation.
        raw_candidates = { ...
            sprintf('%s/zeropad_tse_%d.nii', localpath_char, sub), ...
            sprintf('%s/zeropad_tse_%d_pre.nii', localpath_char, sub), ...
            sprintf('%s/zeropad_tse_%d_post.nii', localpath_char, sub)};
        rawsubjnm = load_first_available(raw_candidates, 'raw TSE');
        rawnm = double(rawsubjnm.img);
        rawnm(rawnm == -999) = NaN;

        % SN ROI mask (required)
        subjmask = load_untouch_nii(char(pth.roi_mask));
        mask = double(subjmask.img) > 0;
        if ~any(mask(:))
            error('SN mask is empty for %d', sub);
        end

        % Optional CNR control mask – probe zeropad/non-zeropad variants so we still
        % compute CNR even when one flavor is missing.
        controlmask = [];
        if isfile(char(pth.cnr_mask_zp))
            cm = load_untouch_nii(char(pth.cnr_mask_zp));
            controlmask = double(cm.img) > 0;
            if size(controlmask,3) > 100 && isfile(char(pth.cnr_mask_nozp))
                cm = load_untouch_nii(char(pth.cnr_mask_nozp));
                controlmask = double(cm.img) > 0;
            end
        elseif isfile(char(pth.cnr_mask_nozp))
            cm = load_untouch_nii(char(pth.cnr_mask_nozp));
            controlmask = double(cm.img) > 0;
        end

        % Identify slices containing SN voxels—needed for the robust scaling window.
        whereStruct = find(squeeze(any(any(mask,1),2)));
        if isempty(whereStruct)
            error('No slices contain SN voxels for %d', sub);
        end

        % ---------- CNR from raw NM ----------
        SNvalues = rawnm(mask);
        meanSN = mean(SNvalues,'omitnan');
        if ~isempty(controlmask)
            CP = rawnm(controlmask);
            meanCP = mean(CP,'omitnan');
            stdCP  = std(CP,0,'omitnan');
            individualCNR = (meanSN - meanCP) / stdCP;
        else
            meanCP = NaN;
            stdCP  = NaN;
            individualCNR = NaN;
        end

        % ---------- Robust min/max scaling on corrected image ----------
        % Compute scaled NM intensities relative to background slices (mu +/- 3SD)
        % so stats remain comparable across acquisitions.
        nmDummy = nm;
        nmDummy(mask) = NaN;
        nmDummyT = nmDummy(:,:,whereStruct);
        nmDummyV = reshape(nmDummyT,1,[]);
        muD = mean(nmDummyV,'omitnan');
        sdD = std(nmDummyV,0,'omitnan');
        maxN = muD + 3*sdD;
        minN = muD - 3*sdD;
        nm_scaled = (nm - minN) / (maxN - minN);

        % ---------- Stats from scaled image ----------
        vox_idx = find(mask);
        S = nm_scaled(vox_idx);
        total = sum(S,'omitnan');
        number = numel(S);
        avS = mean(S,'omitnan');
        S_sorted = sort(S(~isnan(S)));
        if isempty(S_sorted)
            minS = NaN;
            maxS = NaN;
        else
            k = min(10, numel(S_sorted));
            minS = mean(S_sorted(1:k),'omitnan');
            maxS = mean(S_sorted(end-k+1:end),'omitnan');
        end
        background = nm_scaled(~mask);
        avbackground = mean(background,'omitnan');
        vol_mm3 = number * voxelvol;

        rec = struct();
        rec.sub = sub;
        rec.CNR = individualCNR;
        rec.total_intensity = total;
        rec.intensity = avS;
        rec.background_intensity = avbackground;
        rec.total_volume = vol_mm3;
        rec.intensity_minimum = minS;
        rec.intensity_maximum = maxS;
        rec.mean_raw_sn = meanSN;
        rec.mean_raw_control = meanCP;
        rec.std_raw_control = stdCP;
        results(end+1) = rec; %#ok<SAGROW>

        fprintf('[%d] DONE: mean=%.4f CNR=%.4f vox=%d vol=%.2f mm^3\n', ...
                sub, avS, individualCNR, number, vol_mm3);

    catch ME
        warning('[%d] Failed: %s | %s', sub, ME.identifier, ME.message);
        rec = struct();
        rec.sub = sub;
        rec.CNR = NaN;
        rec.total_intensity = NaN;
        rec.intensity = NaN;
        rec.background_intensity = NaN;
        rec.total_volume = NaN;
        rec.intensity_minimum = NaN;
        rec.intensity_maximum = NaN;
        rec.mean_raw_sn = NaN;
        rec.mean_raw_control = NaN;
        rec.std_raw_control = NaN;
        results(end+1) = rec; %#ok<SAGROW>
    end
end

% ---------------------- Write refreshed table ------------------------------
% Serialize the struct array to the SNHEART analysis directory so the Rmd and
% QC harness pick up the fresh metrics automatically.
outdir = '${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data/analysis';
if ~isfolder(outdir)
    mkdir(outdir);
end
outpath = fullfile(outdir, 'SN_nmdata_autosaved.txt');
outT = struct2table(results);
writetable(outT, outpath);
fprintf('Wrote %s (%d rows)\n', outpath, height(outT));

% ====================== Helper functions ==================================
% Keep utility routines at the bottom so the main loop stays readable. Each
% helper mirrors the legacy MATLAB utilities used in the MOTIP projects.
function base = resolve_existing_path(primary, fallback, label)
% resolve_existing_path  Prefer /Volumes when present (local macOS checkout) but
%                        fall back to the /export mirror on cluster nodes.
    if isfolder(primary)
        base = primary;
        fprintf('Using %s for %s\n', primary, label);
    elseif isfolder(fallback)
        base = fallback;
        fprintf('Using %s for %s (fallback)\n', fallback, label);
    else
        error('Cannot reach either %s or %s for %s', primary, fallback, label);
    end
end

function nii = load_first_available(candidates, label)
% load_first_available  Iterate through possible filenames and load the first
% existing NIfTI. Throws a descriptive error if none exist.
    for ii = 1:numel(candidates)
        candidate = candidates{ii};
        if isfile(candidate)
            nii = load_untouch_nii(candidate);
            return;
        end
    end
    error('None of the %s candidates exist:\n%s', label, strjoin(candidates, '\n'));
end

function pth = build_paths(localpath, sub)
% build_paths  Return a struct with every subject-specific file path the main
% pipeline touches (ROI, CNR masks). Centralizes naming quirks.
    base = char(localpath);
    pth.roi_mask = sprintf('%s/SN/SN_ROI_KW_%d_ZP.nii', base, sub);
    pth.cnr_mask_zp = sprintf('%s/SN/SN_ROI_CNR_KW_%d_ZP.nii', base, sub);
    pth.cnr_mask_nozp = sprintf('%s/SN/SN_ROI_CNR_KW_%d.nii', base, sub);
end
