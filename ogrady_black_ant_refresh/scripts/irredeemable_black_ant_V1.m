% irredeemable_black_ant_V1.m
% -------------------------------------------------------------------------
% PURPOSE
%   Generate a "known-good" Black ANT metrics table that only trusts the
%   freshly rebuilt probabilistic masks created by build_black_ant_masks.sh.
%   - Reads the subject-space SN/SN-control/PT-control probability maps from
%     MRI_data/TSE/probabilistic_masks inside 77_SNHEART.
%   - Thresholds the SN mask using a hard gcount cutoff (default >10) so we
%     mimic the collaborator's "include voxels with at least 10 votes" rule.
%   - Leaves the SN-control and PT-control probability maps unthresholded so
%     both the thresholded/unthresholded SN variants share the same controls.
%   - Computes neuromelanin intensity on the corrected TSE volume using the
%     standard background scaling (mu ± 3*SD) and computes CNR from the raw
%     zeropad TSE relative to the SN-control mask.
%   - Writes a dated metrics file under MRI_data/analysis so SN_ASR can ingest
%     it directly without touching the legacy autosave.
%   - Records both thresholded (> gcount cutoff) and unthresholded (any
%     positive probability) metrics into the output table.
%
% ENVIRONMENT OVERRIDES
%   IRREDEEMABLE_SUBJECTS       Space-delimited list of subject IDs. When empty,
%                              the script scans the probabilistic_masks folder.
%   IRREDEEMABLE_SN_THRESHOLD   Gcount threshold for the SN mask (default 10).
%   IRREDEEMABLE_LOG            Optional path to a log/diary file. When absent
%                              but SLURM_JOB_ID is present, logs are saved to
%                              MRI_data/analysis/logs so batch-node runs keep
%                              their STDOUT/STDERR transcripts.
%
% OUTPUT
%   MRI_data/analysis/SN_nmdata_irredeemable_<yyyymmdd>_thrXX.txt
% -------------------------------------------------------------------------

clear; close all;

%% ---------------------- Resolve paths -------------------------------------
snheart_candidates = {
    '${SNHEART_ROOT:-/path/to/SNHEART}', ...
    '${SNHEART_ROOT:-/path/to/SNHEART}'};
snheart_root = resolve_existing_root(snheart_candidates, 'SNHEART root');
prob_mask_dir = fullfile(snheart_root, 'MRI_data', 'TSE', 'probabilistic_masks');
analysis_dir = fullfile(snheart_root, 'MRI_data', 'analysis');
if ~isfolder(prob_mask_dir)
    error('Cannot find probabilistic mask directory: %s', prob_mask_dir);
end
if ~isfolder(analysis_dir)
    mkdir(analysis_dir);
end

% When running on the batch node capture logs automatically so STDOUT/ERR are
% preserved even if sbatch emails are truncated.
log_path = getenv('IRREDEEMABLE_LOG');
if isempty(log_path)
    slurm_job = getenv('SLURM_JOB_ID');
    if ~isempty(slurm_job)
        log_dir = fullfile(analysis_dir, 'logs');
        if ~isfolder(log_dir)
            mkdir(log_dir);
        end
        log_path = fullfile(log_dir, sprintf('irredeemable_black_ant_%s.log', slurm_job));
    end
end
diary_cleanup = [];
if ~isempty(log_path)
    log_parent = fileparts(log_path);
    if ~isempty(log_parent) && ~isfolder(log_parent)
        mkdir(log_parent);
    end
    try
        diary(log_path);
        diary_cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
        fprintf('Logging to %s\n', log_path);
    catch ME
        warning('Unable to start diary at %s (%s). Continuing without diary.', log_path, ME.message);
    end
end

motip_candidates = {
    '${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data', ...
    '${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data'};
motip_root = resolve_existing_root(motip_candidates, 'MOTIP MRI_data');
motip_tse_root = fullfile(motip_root, 'TSE');

nifti_toolbox = '${MOTIP_ROOT:-/path/to/MOTIP2018}/scripts/tse/NIfTI_20140122';
if ~isfolder(nifti_toolbox)
    error('Missing NIfTI toolbox at %s', nifti_toolbox);
end
addpath(nifti_toolbox);

%% ---------------------- Subject roster ------------------------------------
subjects = discover_subjects(prob_mask_dir, getenv('IRREDEEMABLE_SUBJECTS'));
if isempty(subjects)
    error('No sn_groupmask_in_*.nii files found under %s', prob_mask_dir);
end

sn_threshold = getenv_num('IRREDEEMABLE_SN_THRESHOLD', 10);
sn_prob_subjects = numel(subjects);
if sn_prob_subjects <= 0
    error('No probabilistic masks discovered; cannot derive SN threshold');
end
sn_threshold_prob = sn_threshold / sn_prob_subjects;
if sn_threshold_prob > 1
    warning('Requested SN gcount cutoff %.2f exceeds cohort size %d; clamping to probability 1.0', ...
            sn_threshold, sn_prob_subjects);
    sn_threshold_prob = 1;
elseif sn_threshold_prob < 0
    sn_threshold_prob = 0;
end

sn_prob_fixed = getenv_num('IRREDEEMABLE_SN_PROB_FIXED', 0.1);
if sn_prob_fixed < 0
    sn_prob_fixed = 0;
elseif sn_prob_fixed > 1
    sn_prob_fixed = 1;
end
voxelvol = 0.43 * 0.43 * 3;  % mm^3

fprintf(['Irredeemable Black ANT V1: %d subjects, SN threshold > %.2f gcount ', ...
         '(probability > %.4f); control masks left unthresholded\n'], ...
        sn_prob_subjects, sn_threshold, sn_threshold_prob);
fprintf('Fixed probability SN cutoff: > %.3f\n', sn_prob_fixed);

results = struct('sub',{}, 'CNR',{}, 'CNR_peduncle',{}, ...
                 'CNR_unthresholded',{}, 'CNR_peduncle_unthresholded',{}, ...
                 'total_intensity',{}, 'intensity',{}, 'background_intensity',{}, ...
                 'total_intensity_unthresholded',{}, 'intensity_unthresholded',{}, ...
                 'background_intensity_unthresholded',{}, 'total_volume',{}, ...
                 'total_volume_unthresholded',{}, 'intensity_minimum',{}, ...
                 'intensity_maximum',{}, 'intensity_minimum_unthresholded',{}, ...
                 'intensity_maximum_unthresholded',{}, 'sn_voxels',{}, ...
                 'sn_voxels_unthresholded',{}, 'sn_prob_mean',{}, ...
                 'sn_prob_mean_unthresholded',{}, 'mean_raw_sn',{}, ...
                 'mean_raw_sn_unthresholded',{}, 'mean_raw_control',{}, ...
                 'mean_raw_control_unthresholded',{}, 'std_raw_control',{}, ...
                 'std_raw_control_unthresholded',{}, 'mean_raw_peduncle',{}, ...
                 'mean_raw_peduncle_unthresholded',{}, 'notes',{});
results_prob = results;

%% ---------------------- Subject loop --------------------------------------
for idx = 1:numel(subjects)
    sub = subjects{idx};
    fprintf('[%d/%d] Subject %s\n', idx, numel(subjects), sub);
    note_entries = {};
    try
        raw_candidates = {
            fullfile(motip_tse_root, sprintf('zeropad_tse_%s.nii', sub)), ...
            fullfile(motip_tse_root, sprintf('zeropad_tse_%s.nii.gz', sub)), ...
            fullfile(motip_tse_root, sprintf('%s_TSE_zeropad.nii', sub)), ...
            fullfile(motip_tse_root, sprintf('%s_TSE_zeropad.nii.gz', sub))};
        corr_candidates = {
            fullfile(motip_tse_root, 'SN', sprintf('corrected_tse_DC_%s.nii', sub)), ...
            fullfile(motip_tse_root, 'SN', sprintf('corrected_tse_DC_%s.nii.gz', sub)), ...
            fullfile(motip_tse_root, 'SN', sprintf('corrected_tse_DC_%s_pre.nii', sub)), ...
            fullfile(motip_tse_root, 'SN', sprintf('corrected_tse_DC_%s_post.nii', sub))};

        raw_path = pick_first_existing(raw_candidates, 'raw TSE');
        corr_path = pick_first_existing(corr_candidates, 'corrected TSE');

        sn_mask_prob = load_mask(prob_mask_dir, 'sn', sub);
        snctrl_prob = load_mask(prob_mask_dir, 'snctrl', sub);
        ptctrl_prob = load_mask(prob_mask_dir, 'ptctrl', sub);

        rawnm = load_nifti_image(raw_path);
        corrnm = load_nifti_image(corr_path);

        sn_mask = sn_mask_prob > sn_threshold_prob;
        if ~any(sn_mask(:))
            error('SN mask empty after applying threshold %.2f', sn_threshold);
        end
        sn_mask_full = sn_mask_prob > 0;
        snctrl_mask = snctrl_prob > 0;
        ptctrl_mask = ptctrl_prob > 0;
        if ~any(snctrl_mask(:))
            note_entries{end+1} = 'snctrl_empty'; %#ok<AGROW>
        end
        if ~any(ptctrl_mask(:))
            note_entries{end+1} = 'ptctrl_empty'; %#ok<AGROW>
        end
        [thr_metrics, thr_notes] = compute_mask_metrics(sn_mask, sn_mask_prob, ...
            rawnm, corrnm, snctrl_mask, ptctrl_mask, voxelvol, 'thr');
        note_entries = [note_entries, thr_notes]; %#ok<AGROW>
        [full_metrics, full_notes] = compute_mask_metrics(sn_mask_full, sn_mask_prob, ...
            rawnm, corrnm, snctrl_mask, ptctrl_mask, voxelvol, 'unthr');
        note_entries = [note_entries, full_notes]; %#ok<AGROW>
        prob_mask = sn_mask_prob > sn_prob_fixed;
        [prob_metrics, prob_notes] = compute_mask_metrics(prob_mask, sn_mask_prob, ...
            rawnm, corrnm, snctrl_mask, ptctrl_mask, voxelvol, 'prob');

        rec = struct();
        rec.sub = str2double(sub);
        rec.CNR = thr_metrics.CNR;
        rec.CNR_peduncle = thr_metrics.CNR_peduncle;
        rec.CNR_unthresholded = full_metrics.CNR;
        rec.CNR_peduncle_unthresholded = full_metrics.CNR_peduncle;
        rec.total_intensity = thr_metrics.total_intensity;
        rec.intensity = thr_metrics.intensity;
        rec.background_intensity = thr_metrics.background_intensity;
        rec.total_intensity_unthresholded = full_metrics.total_intensity;
        rec.intensity_unthresholded = full_metrics.intensity;
        rec.background_intensity_unthresholded = full_metrics.background_intensity;
        rec.total_volume = thr_metrics.total_volume;
        rec.total_volume_unthresholded = full_metrics.total_volume;
        rec.intensity_minimum = thr_metrics.intensity_minimum;
        rec.intensity_maximum = thr_metrics.intensity_maximum;
        rec.intensity_minimum_unthresholded = full_metrics.intensity_minimum;
        rec.intensity_maximum_unthresholded = full_metrics.intensity_maximum;
        rec.sn_voxels = thr_metrics.sn_voxels;
        rec.sn_voxels_unthresholded = full_metrics.sn_voxels;
        rec.sn_prob_mean = thr_metrics.sn_prob_mean;
        rec.sn_prob_mean_unthresholded = full_metrics.sn_prob_mean;
        rec.mean_raw_sn = thr_metrics.mean_raw_sn;
        rec.mean_raw_sn_unthresholded = full_metrics.mean_raw_sn;
        rec.mean_raw_control = thr_metrics.mean_raw_control;
        rec.mean_raw_control_unthresholded = full_metrics.mean_raw_control;
        rec.std_raw_control = thr_metrics.std_raw_control;
        rec.std_raw_control_unthresholded = full_metrics.std_raw_control;
        rec.mean_raw_peduncle = thr_metrics.mean_raw_peduncle;
        rec.mean_raw_peduncle_unthresholded = full_metrics.mean_raw_peduncle;
        note_entries = note_entries(~cellfun('isempty', note_entries));
        rec.notes = strjoin(unique(note_entries), ',');
        results(end+1) = rec; %#ok<SAGROW>

        rec_prob = rec;
        rec_prob.CNR = prob_metrics.CNR;
        rec_prob.CNR_peduncle = prob_metrics.CNR_peduncle;
        rec_prob.total_intensity = prob_metrics.total_intensity;
        rec_prob.intensity = prob_metrics.intensity;
        rec_prob.background_intensity = prob_metrics.background_intensity;
        rec_prob.total_volume = prob_metrics.total_volume;
        rec_prob.intensity_minimum = prob_metrics.intensity_minimum;
        rec_prob.intensity_maximum = prob_metrics.intensity_maximum;
        rec_prob.sn_voxels = prob_metrics.sn_voxels;
        rec_prob.sn_prob_mean = prob_metrics.sn_prob_mean;
        rec_prob.mean_raw_sn = prob_metrics.mean_raw_sn;
        rec_prob.mean_raw_control = prob_metrics.mean_raw_control;
        rec_prob.std_raw_control = prob_metrics.std_raw_control;
        rec_prob.mean_raw_peduncle = prob_metrics.mean_raw_peduncle;
        rec_prob.notes = strjoin(unique(prob_notes), ',');
        results_prob(end+1) = rec_prob; %#ok<SAGROW>

        fprintf('    vox=%d CNR=%.3f intensity=%.3f\n', ...
                thr_metrics.sn_voxels, thr_metrics.CNR, thr_metrics.intensity);
    catch ME
        warning('[%s] Failed: %s | %s', sub, ME.identifier, ME.message);
        rec = struct();
        rec.sub = str2double(sub);
        rec.CNR = NaN;
        rec.CNR_peduncle = NaN;
        rec.CNR_unthresholded = NaN;
        rec.CNR_peduncle_unthresholded = NaN;
        rec.total_intensity = NaN;
        rec.intensity = NaN;
        rec.background_intensity = NaN;
        rec.total_intensity_unthresholded = NaN;
        rec.intensity_unthresholded = NaN;
        rec.background_intensity_unthresholded = NaN;
        rec.total_volume = NaN;
        rec.total_volume_unthresholded = NaN;
        rec.intensity_minimum = NaN;
        rec.intensity_maximum = NaN;
        rec.intensity_minimum_unthresholded = NaN;
        rec.intensity_maximum_unthresholded = NaN;
        rec.sn_voxels = NaN;
        rec.sn_voxels_unthresholded = NaN;
        rec.sn_prob_mean = NaN;
        rec.sn_prob_mean_unthresholded = NaN;
        rec.mean_raw_sn = NaN;
        rec.mean_raw_sn_unthresholded = NaN;
        rec.mean_raw_control = NaN;
        rec.std_raw_control = NaN;
        rec.mean_raw_control_unthresholded = NaN;
        rec.std_raw_control_unthresholded = NaN;
        rec.mean_raw_peduncle = NaN;
        rec.mean_raw_peduncle_unthresholded = NaN;
        rec.notes = strjoin([note_entries {ME.identifier}], ',');
        results(end+1) = rec; %#ok<SAGROW>
        results_prob(end+1) = rec; %#ok<SAGROW>
        continue;
    end
end

%% ---------------------- Write output --------------------------------------
out_table = struct2table(results);
if isempty(out_table)
    error('No output rows were generated.');
end
timestamp = datestr(now, 'yyyymmdd');
out_file = fullfile(analysis_dir, sprintf('SN_nmdata_irredeemable_%s_thr%02d.txt', ...
    timestamp, round(sn_threshold)));
writetable(out_table, out_file);
fprintf('Wrote %s (%d rows)\n', out_file, height(out_table));

out_table_prob = struct2table(results_prob);
if isempty(out_table_prob)
    error('No probability-threshold rows were generated.');
end
prob_slug = round(sn_prob_fixed * 100);
out_file_prob = fullfile(analysis_dir, sprintf('SN_nmdata_irredeemable_%s_prob%02d.txt', ...
    timestamp, prob_slug));
writetable(out_table_prob, out_file_prob);
fprintf('Wrote %s (%d rows)\n', out_file_prob, height(out_table_prob));

%% ====================== Helper functions ==================================
function [metrics, note_entries] = compute_mask_metrics(sn_mask, sn_mask_prob, ...
    raw_vol, corr_vol, snctrl_mask, ptctrl_mask, voxelvol, note_tag)

    metrics = struct('CNR', NaN, 'CNR_peduncle', NaN, ...
        'total_intensity', NaN, 'intensity', NaN, ...
        'background_intensity', NaN, 'total_volume', NaN, ...
        'intensity_minimum', NaN, 'intensity_maximum', NaN, ...
        'sn_voxels', NaN, 'sn_prob_mean', NaN, ...
        'mean_raw_sn', NaN, 'mean_raw_control', NaN, ...
        'std_raw_control', NaN, 'mean_raw_peduncle', NaN);
    note_entries = {};

    if ~any(sn_mask(:))
        note_entries{end+1} = sprintf('%s_sn_empty', note_tag); %#ok<AGROW>
        return;
    end

    sn_vals_raw = raw_vol(sn_mask);
    metrics.mean_raw_sn = mean(sn_vals_raw, 'omitnan');
    ctrl_vals_raw = raw_vol(snctrl_mask);
    metrics.mean_raw_control = mean(ctrl_vals_raw, 'omitnan');
    metrics.std_raw_control = std(ctrl_vals_raw, 0, 'omitnan');
    if ~isnan(metrics.mean_raw_control) && metrics.std_raw_control > 0
        metrics.CNR = (metrics.mean_raw_sn - metrics.mean_raw_control) / metrics.std_raw_control;
    else
        note_entries{end+1} = sprintf('%s_cnr_nan', note_tag); %#ok<AGROW>
    end

    ped_vals_raw = raw_vol(ptctrl_mask);
    metrics.mean_raw_peduncle = mean(ped_vals_raw, 'omitnan');
    if ~isnan(metrics.mean_raw_peduncle) && metrics.std_raw_control > 0
        metrics.CNR_peduncle = (metrics.mean_raw_peduncle - metrics.mean_raw_control) / metrics.std_raw_control;
    end

    metrics.sn_voxels = nnz(sn_mask);
    metrics.total_volume = metrics.sn_voxels * voxelvol;
    metrics.sn_prob_mean = mean(sn_mask_prob(sn_mask), 'omitnan');

    slice_mask = squeeze(any(any(sn_mask, 1), 2));
    slice_idx = find(slice_mask);
    if isempty(slice_idx)
        note_entries{end+1} = sprintf('%s_no_slices', note_tag); %#ok<AGROW>
        return;
    end

    nmDummy = corr_vol;
    nmDummy(sn_mask) = NaN;
    nmDummySlices = nmDummy(:, :, slice_idx);
    background_vec = nmDummySlices(:);
    muD = mean(background_vec, 'omitnan');
    sdD = std(background_vec, 0, 'omitnan');
    if isnan(muD) || isnan(sdD) || sdD == 0
        note_entries{end+1} = sprintf('%s_background_invalid', note_tag); %#ok<AGROW>
        return;
    end
    minN = muD - 3 * sdD;
    maxN = muD + 3 * sdD;
    if ~isfinite(minN) || ~isfinite(maxN) || maxN <= minN
        note_entries{end+1} = sprintf('%s_background_invalid', note_tag); %#ok<AGROW>
        return;
    end

    nm_scaled = (corr_vol - minN) / (maxN - minN);
    sn_scaled = nm_scaled(sn_mask);
    metrics.total_intensity = sum(sn_scaled, 'omitnan');
    metrics.intensity = mean(sn_scaled, 'omitnan');
    background_vals = nm_scaled(~sn_mask);
    metrics.background_intensity = mean(background_vals, 'omitnan');

    sorted_vals = sort(sn_scaled(~isnan(sn_scaled)));
    if isempty(sorted_vals)
        metrics.intensity_minimum = NaN;
        metrics.intensity_maximum = NaN;
    else
        k = min(10, numel(sorted_vals));
        metrics.intensity_minimum = mean(sorted_vals(1:k), 'omitnan');
        metrics.intensity_maximum = mean(sorted_vals(end-k+1:end), 'omitnan');
    end
end

function root = resolve_existing_root(candidates, label)
    for ii = 1:numel(candidates)
        if isfolder(candidates{ii})
            root = candidates{ii};
            fprintf('Using %s for %s\n', root, label);
            return;
        end
    end
    error('Cannot resolve path for %s', label);
end

function subs = discover_subjects(mask_dir, override_str)
    if nargin > 1 && ~isempty(strtrim(override_str))
        tokens = regexp(strtrim(override_str), '\s+', 'split');
        subs = tokens(~cellfun('isempty', tokens));
        return;
    end
    listing = dir(fullfile(mask_dir, 'sn_groupmask_in_*.nii*'));
    subs = unique(arrayfun(@(d) extract_id(d.name), listing, 'UniformOutput', false));
    subs = subs(~cellfun('isempty', subs));
    subs = sort(subs);
end

function id = extract_id(filename)
    tokens = regexp(filename, 'sn_groupmask_in_(.+?)\.nii', 'tokens');
    if isempty(tokens)
        tokens = regexp(filename, 'sn_groupmask_in_(.+?)\.nii\.gz', 'tokens');
    end
    if isempty(tokens)
        id = '';
    else
        id = tokens{1}{1};
    end
end

function val = getenv_num(name, default_val)
    raw = getenv(name);
    if isempty(raw)
        val = default_val;
    else
        val = str2double(raw);
        if isnan(val)
            val = default_val;
        end
    end
end

function path = pick_first_existing(candidates, label)
    for ii = 1:numel(candidates)
        if isfile(candidates{ii})
            path = candidates{ii};
            return;
        end
    end
    error('Missing %s (checked %d candidates)', label, numel(candidates));
end

function data = load_nifti_image(path)
    nii = load_untouch_nii(path);
    data = double(nii.img);
    data(data == -999) = NaN;
end

function prob_data = load_mask(mask_dir, label, sub)
    candidates = {
        fullfile(mask_dir, sprintf('%s_groupmask_in_%s.nii', label, sub)), ...
        fullfile(mask_dir, sprintf('%s_groupmask_in_%s.nii.gz', label, sub))};
    prob_data = [];
    for ii = 1:numel(candidates)
        if isfile(candidates{ii})
            nii = load_untouch_nii(candidates{ii});
            prob_data = double(nii.img);
            return;
        end
    end
    error('Missing %s mask for subject %s', label, sub);
end
