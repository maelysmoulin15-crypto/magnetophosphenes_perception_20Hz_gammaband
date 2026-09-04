%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESSING.M
%
% Trial-level EEG spectral processing for 20-Hz magnetic-field exposure.
%
% This script expects 5-s, average-referenced EEG epochs created by the
% preprocessing stage. For each electrode and trial, it computes:
%   1. Raw 30-80 Hz power
%   2. Specparam-adjusted 30-80 Hz power
%   3. Log-log-adjusted 30-80 Hz power
%
% Processing sequence:
%   - Apply a symmetric Hann taper to the complete 5-s epoch
%   - Apply a rectangular 30-80 Hz DFT-domain mask
%   - Retain the 1-4 s interval
%   - Apply zero-phase IIR notch filters at 40, 50, and 60 Hz
%   - Downsample from 10 kHz to 2 kHz
%   - Estimate the PSD using Welch's method
%   - Integrate raw and aperiodic-adjusted power over 42 retained bins
%
% Requirements:
%   - MATLAB Signal Processing Toolbox
%   - sinusoidfilter.m
%   - MATLAB wrapper for specparam 2.x
%   - Python environment containing specparam 2.x
%
% The input signal units are preserved. If EEG_segment is expressed in
% microvolts, PSD values are in microvolts squared per hertz and integrated
% power values are in microvolts squared.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
clc;
close all;

%% Repository paths

repository_dir = fileparts(mfilename('fullpath'));

epoch_dir = fullfile(repository_dir, 'data', 'epochs');
results_dir = fullfile(repository_dir, 'results', 'processing');
specparam_matlab_dir = fullfile(repository_dir, 'external', 'specparam_mat');

if ~isfolder(epoch_dir)
    error('Epoch directory not found: %s', epoch_dir);
end

if ~isfolder(results_dir)
    mkdir(results_dir);
end

if isfolder(specparam_matlab_dir)
    addpath(specparam_matlab_dir);
end

if exist('sinusoidfilter', 'file') ~= 2
    error(['sinusoidfilter.m was not found. Place it in the repository ', ...
           'root or add its directory to the MATLAB path.']);
end

if exist('specparam', 'file') ~= 2
    error(['The MATLAB specparam wrapper was not found. Place it in ', ...
           'external/specparam_mat or add its directory to the path.']);
end

%% Python environment

% Optionally define the SPECPARAM_PYTHON environment variable before
% starting MATLAB. If it is not defined, MATLAB's current Python
% configuration is used.
python_executable = getenv('SPECPARAM_PYTHON');

if ~isempty(python_executable)
    pyenv('Version', python_executable);
end

try
    importlib_metadata = py.importlib.import_module('importlib.metadata');
    specparam_version = string(importlib_metadata.version('specparam'));
catch ME
    error('Unable to import the Python specparam package: %s', ME.message);
end

if ~startsWith(specparam_version, "2.")
    error('specparam 2.x is required; detected version: %s', ...
        specparam_version);
end

fprintf('[INIT] specparam version: %s\n', specparam_version);

%% Analysis configuration

stimulation_frequency_hz = 20;
intensities_to_process_mt = [0, 5, 50];
expected_trials_per_condition = 5;

original_sampling_rate_hz = 10000;
epoch_duration_s = 5;
retained_interval_s = [1, 4];
downsampling_factor = 5;

gamma_band_hz = [30, 80];
notch_frequencies_hz = [40, 50, 60];
notch_half_width_hz = 1.5;

welch_window_duration_s = 1;
welch_overlap_fraction = 0.50;

verbose = true;
progress_interval = 250;
save_frequency_resolved_spectra = true;

%% Specparam configuration

specparam_settings = struct();
specparam_settings.peak_width_limits = [1, 12];
specparam_settings.max_n_peaks = 6;
specparam_settings.min_peak_height = 0;
specparam_settings.peak_threshold = 2;
specparam_settings.aperiodic_mode = 'fixed';
specparam_settings.verbose = 0;

%% Find and select epoch files

all_files = dir(fullfile(epoch_dir, '*.mat'));

if isempty(all_files)
    error('No MAT files were found in: %s', epoch_dir);
end

% Supported example:
% Subject07_OCCIPITAL_20Hz_Oz_50mT_trial3.mat
% The French prefix "Sujet" is also accepted for compatibility with the
% original dataset.
filename_pattern = [ ...
    '^(?:Subject|Sujet)(?<subject>\d+)_' ...
    '(?<configuration>[A-Za-z]+)_' ...
    '(?<frequency>\d+)Hz_' ...
    '(?<electrode>.+)_' ...
    '(?<intensity>\d+)mT_' ...
    'trial(?<trial>\d+)\.mat$'];

file_metadata = cell(numel(all_files), 1);
selected_file = false(numel(all_files), 1);

for i_file = 1:numel(all_files)
    metadata = regexp(all_files(i_file).name, filename_pattern, ...
        'names', 'once');

    if isempty(metadata)
        continue;
    end

    frequency_hz = str2double(metadata.frequency);
    intensity_mt = str2double(metadata.intensity);

    if frequency_hz == stimulation_frequency_hz && ...
            ismember(intensity_mt, intensities_to_process_mt)
        file_metadata{i_file} = metadata;
        selected_file(i_file) = true;
    end
end

epoch_files = all_files(selected_file);
file_metadata = file_metadata(selected_file);
number_of_files = numel(epoch_files);

if number_of_files == 0
    error(['No files matched the expected filename pattern, frequency, ', ...
           'and intensity selection.']);
end

fprintf('[INIT] Selected %d epoch files.\n', number_of_files);

%% Preallocate trial-level outputs

Subject = nan(number_of_files, 1);
Electrode = strings(number_of_files, 1);
Intensity_mT = nan(number_of_files, 1);
Trial = nan(number_of_files, 1);
OriginalSamplingRate_Hz = nan(number_of_files, 1);
EffectiveSamplingRate_Hz = nan(number_of_files, 1);

RawGammaPower = nan(number_of_files, 1);
SpecparamAdjustedGammaPower = nan(number_of_files, 1);
LogLogAdjustedGammaPower = nan(number_of_files, 1);

SpecparamExponent = nan(number_of_files, 1);
SpecparamOffsetLog10 = nan(number_of_files, 1);
SpecparamOffsetLinear = nan(number_of_files, 1);
SpecparamR2 = nan(number_of_files, 1);

LogLogExponent = nan(number_of_files, 1);
LogLogOffsetLog10 = nan(number_of_files, 1);
LogLogOffsetLinear = nan(number_of_files, 1);

EpochFile = strings(number_of_files, 1);

if save_frequency_resolved_spectra
    Frequency_Hz = cell(number_of_files, 1);
    RawSpectrum = cell(number_of_files, 1);
    SpecparamAdjustedSpectrum = cell(number_of_files, 1);
    LogLogAdjustedSpectrum = cell(number_of_files, 1);
end

%% Process each electrode-level epoch

processing_start = tic;

for i_file = 1:number_of_files

    metadata = file_metadata{i_file};
    file_path = fullfile(epoch_files(i_file).folder, ...
        epoch_files(i_file).name);

    subject = str2double(metadata.subject);
    electrode = string(metadata.electrode);
    intensity_mt = str2double(metadata.intensity);
    trial = str2double(metadata.trial);

    loaded_data = load(file_path);

    if ~isfield(loaded_data, 'EEG_segment')
        error('EEG_segment is missing from: %s', file_path);
    end

    eeg_epoch = double(loaded_data.EEG_segment(:));

    if isfield(loaded_data, 'meta') && ...
            isfield(loaded_data.meta, 'Fs')
        sampling_rate_hz = double(loaded_data.meta.Fs);
    elseif isfield(loaded_data, 'time_segment')
        time_vector = double(loaded_data.time_segment(:));
        sampling_rate_hz = 1 / median(diff(time_vector));
    else
        error('Sampling rate information is missing from: %s', file_path);
    end

    expected_samples = round(epoch_duration_s * ...
        original_sampling_rate_hz);

    if abs(sampling_rate_hz - original_sampling_rate_hz) > 1e-6
        error('Expected a 10-kHz sampling rate in: %s', file_path);
    end

    if numel(eeg_epoch) ~= expected_samples
        error(['Expected exactly %d samples but found %d samples in: %s'], ...
            expected_samples, numel(eeg_epoch), file_path);
    end

    if any(~isfinite(eeg_epoch))
        error('The epoch contains non-finite samples: %s', file_path);
    end

    % Apply a symmetric Hann taper to the complete 5-s epoch.
    epoch_taper = hann(expected_samples, 'symmetric');
    tapered_epoch = eeg_epoch .* epoch_taper;

    % Apply a real, symmetric, rectangular DFT-domain mask. With the
    % present 5-s epochs, the DFT spacing is 0.2 Hz. The 30-Hz bin is
    % rejected, the first retained bin is 30.2 Hz, and 80 Hz is retained.
    gamma_filtered_epoch = sinusoidfilter( ...
        tapered_epoch, ...
        sampling_rate_hz, ...
        gamma_band_hz(1), ...
        0, ...
        gamma_band_hz(2), ...
        0);

    dft_resolution_hz = sampling_rate_hz / expected_samples;
    assert(abs(dft_resolution_hz - 0.2) < 1e-12, ...
        'The expected DFT resolution is 0.2 Hz.');

    % Retain [1, 4) s: exactly 3 s and 30,000 samples at 10 kHz.
    first_retained_sample = ...
        round(retained_interval_s(1) * sampling_rate_hz) + 1;
    last_retained_sample = ...
        round(retained_interval_s(2) * sampling_rate_hz);

    retained_epoch = gamma_filtered_epoch( ...
        first_retained_sample:last_retained_sample);

    assert(numel(retained_epoch) == round(3 * sampling_rate_hz), ...
        'The retained interval must contain exactly 3 s of data.');

    % Combine three second-order IIR notch filters and apply them with
    % filtfilt to obtain a zero-phase response. Each notch has a nominal
    % total bandwidth of 3 Hz.
    combined_numerator = 1;
    combined_denominator = 1;

    for notch_frequency_hz = notch_frequencies_hz
        normalized_frequency = ...
            notch_frequency_hz / (sampling_rate_hz / 2);
        normalized_bandwidth = ...
            (2 * notch_half_width_hz) / (sampling_rate_hz / 2);

        [notch_numerator, notch_denominator] = iirnotch( ...
            normalized_frequency, normalized_bandwidth);

        combined_numerator = conv( ...
            combined_numerator, notch_numerator);
        combined_denominator = conv( ...
            combined_denominator, notch_denominator);
    end

    notch_filtered_epoch = filtfilt( ...
        combined_numerator, combined_denominator, retained_epoch);

    % Downsample from 10 kHz to 2 kHz. MATLAB resample applies an
    % anti-aliasing filter internally.
    downsampled_epoch = resample( ...
        notch_filtered_epoch, 1, downsampling_factor);
    effective_sampling_rate_hz = ...
        sampling_rate_hz / downsampling_factor;

    assert(abs(effective_sampling_rate_hz - 2000) < 1e-6, ...
        'The expected effective sampling rate is 2 kHz.');
    assert(numel(downsampled_epoch) == ...
        round(3 * effective_sampling_rate_hz), ...
        'The downsampled epoch must contain exactly 3 s of data.');

    % Welch PSD: 1-s periodic Hann windows with 50% overlap. Setting the
    % FFT length equal to the window length gives a 1-Hz frequency grid.
    welch_window_samples = round( ...
        welch_window_duration_s * effective_sampling_rate_hz);
    welch_overlap_samples = round( ...
        welch_overlap_fraction * welch_window_samples);
    welch_fft_length = welch_window_samples;

    [power_spectrum, frequencies_hz] = pwelch( ...
        downsampled_epoch, ...
        hann(welch_window_samples, 'periodic'), ...
        welch_overlap_samples, ...
        welch_fft_length, ...
        effective_sampling_rate_hz, ...
        'psd');

    number_of_welch_segments = 1 + floor( ...
        (numel(downsampled_epoch) - welch_window_samples) / ...
        (welch_window_samples - welch_overlap_samples));

    assert(number_of_welch_segments == 5, ...
        'Each epoch must yield exactly five Welch segments.');

    frequency_resolution_hz = median(diff(frequencies_hz));

    assert(max(abs(diff(frequencies_hz) - ...
        frequency_resolution_hz)) < 1e-10, ...
        'The Welch frequency vector is not regularly spaced.');
    assert(abs(frequency_resolution_hz - 1) < 1e-10, ...
        'The expected Welch frequency resolution is 1 Hz.');

    % Define frequency masks. The notch-affected bins remain in the
    % complete PSD and are supplied to specparam. They are excluded only
    % from raw integration, log-log fitting, and adjusted integration.
    gamma_mask = ...
        frequencies_hz >= gamma_band_hz(1) & ...
        frequencies_hz <= gamma_band_hz(2);

    notch_mask = false(size(frequencies_hz));
    for notch_frequency_hz = notch_frequencies_hz
        notch_mask = notch_mask | ...
            abs(frequencies_hz - notch_frequency_hz) <= ...
            notch_half_width_hz;
    end

    integration_mask = gamma_mask & ...
        ~notch_mask & isfinite(power_spectrum);

    log_log_fit_mask = integration_mask & power_spectrum > 0;

    assert(nnz(gamma_mask) == 51, ...
        'Expected 51 Welch bins between 30 and 80 Hz.');
    assert(nnz(integration_mask) == 42, ...
        'Expected 42 bins after notch-region exclusion.');

    % 1. Raw gamma power.
    raw_gamma_power = sum(power_spectrum(integration_mask)) * ...
        frequency_resolution_hz;

    % 2. Specparam-adjusted gamma power. Specparam receives the complete
    % 30-80 Hz grid, including the notch-attenuated bins.
    specparam_fit_mask = ...
        gamma_mask & ...
        isfinite(power_spectrum) & ...
        power_spectrum > 0;

    assert(nnz(specparam_fit_mask) == 51, ...
        'Specparam requires the complete 51-bin gamma spectrum.');

    specparam_frequencies = frequencies_hz(specparam_fit_mask).';
    specparam_spectrum = power_spectrum(specparam_fit_mask).';

    try
        spectral_model = specparam( ...
            specparam_frequencies, ...
            specparam_spectrum, ...
            gamma_band_hz, ...
            specparam_settings, ...
            true);
    catch ME
        error('Specparam failed for %s: %s', file_path, ME.message);
    end

    specparam_offset_log10 = spectral_model.aperiodic_params(1);
    specparam_exponent = spectral_model.aperiodic_params(2);
    specparam_offset_linear = 10.^specparam_offset_log10;
    specparam_r_squared = spectral_model.r_squared;

    specparam_background = nan(size(frequencies_hz));
    specparam_background(specparam_fit_mask) = ...
        10.^(spectral_model.ap_fit(:));

    specparam_adjusted_spectrum = ...
        power_spectrum - specparam_background;
    specparam_adjusted_spectrum( ...
        ~isfinite(specparam_adjusted_spectrum) | ...
        specparam_adjusted_spectrum < 0) = 0;

    specparam_integration_mask = ...
        integration_mask & isfinite(specparam_background);

    assert(nnz(specparam_integration_mask) == 42, ...
        'Expected 42 bins for specparam-adjusted integration.');

    specparam_adjusted_gamma_power = sum( ...
        specparam_adjusted_spectrum(specparam_integration_mask)) * ...
        frequency_resolution_hz;

    % 3. Log-log-adjusted gamma power. Notch-region bins are excluded
    % before fitting log10(power) as a linear function of log10(frequency).
    log_log_coefficients = polyfit( ...
        log10(frequencies_hz(log_log_fit_mask)), ...
        log10(power_spectrum(log_log_fit_mask)), ...
        1);

    log_log_exponent = -log_log_coefficients(1);
    log_log_offset_log10 = log_log_coefficients(2);
    log_log_offset_linear = 10.^log_log_offset_log10;

    log_log_background = nan(size(frequencies_hz));
    log_log_background(gamma_mask) = 10.^polyval( ...
        log_log_coefficients, ...
        log10(frequencies_hz(gamma_mask)));

    log_log_adjusted_spectrum = power_spectrum - log_log_background;
    log_log_adjusted_spectrum( ...
        ~isfinite(log_log_adjusted_spectrum) | ...
        log_log_adjusted_spectrum < 0) = 0;

    log_log_adjusted_gamma_power = sum( ...
        log_log_adjusted_spectrum(integration_mask)) * ...
        frequency_resolution_hz;

    % Store trial-level results.
    Subject(i_file) = subject;
    Electrode(i_file) = electrode;
    Intensity_mT(i_file) = intensity_mt;
    Trial(i_file) = trial;
    OriginalSamplingRate_Hz(i_file) = sampling_rate_hz;
    EffectiveSamplingRate_Hz(i_file) = effective_sampling_rate_hz;

    RawGammaPower(i_file) = raw_gamma_power;
    SpecparamAdjustedGammaPower(i_file) = ...
        specparam_adjusted_gamma_power;
    LogLogAdjustedGammaPower(i_file) = ...
        log_log_adjusted_gamma_power;

    SpecparamExponent(i_file) = specparam_exponent;
    SpecparamOffsetLog10(i_file) = specparam_offset_log10;
    SpecparamOffsetLinear(i_file) = specparam_offset_linear;
    SpecparamR2(i_file) = specparam_r_squared;

    LogLogExponent(i_file) = log_log_exponent;
    LogLogOffsetLog10(i_file) = log_log_offset_log10;
    LogLogOffsetLinear(i_file) = log_log_offset_linear;

    EpochFile(i_file) = string(epoch_files(i_file).name);

    if save_frequency_resolved_spectra
        Frequency_Hz{i_file} = frequencies_hz(integration_mask);
        RawSpectrum{i_file} = power_spectrum(integration_mask);
        SpecparamAdjustedSpectrum{i_file} = ...
            specparam_adjusted_spectrum(integration_mask);
        LogLogAdjustedSpectrum{i_file} = ...
            log_log_adjusted_spectrum(integration_mask);
    end

    if verbose && (i_file == 1 || ...
            mod(i_file, progress_interval) == 0 || ...
            i_file == number_of_files)
        elapsed_s = toc(processing_start);
        fprintf('[PROCESSING] %d/%d files | elapsed: %.1f s\n', ...
            i_file, number_of_files, elapsed_s);
    end
end

%% Build the trial-level results table

results_table = table( ...
    Subject, ...
    Electrode, ...
    Intensity_mT, ...
    Trial, ...
    OriginalSamplingRate_Hz, ...
    EffectiveSamplingRate_Hz, ...
    RawGammaPower, ...
    SpecparamAdjustedGammaPower, ...
    LogLogAdjustedGammaPower, ...
    SpecparamExponent, ...
    SpecparamOffsetLog10, ...
    SpecparamOffsetLinear, ...
    SpecparamR2, ...
    LogLogExponent, ...
    LogLogOffsetLog10, ...
    LogLogOffsetLinear, ...
    EpochFile);

results_table = sortrows(results_table, ...
    {'Subject', 'Electrode', 'Intensity_mT', 'Trial'});

%% Verify trial completeness

[group_id, group_subject, group_electrode, group_intensity] = ...
    findgroups( ...
        results_table.Subject, ...
        results_table.Electrode, ...
        results_table.Intensity_mT);

group_trial_count = splitapply(@numel, results_table.Trial, group_id);
group_unique_trial_count = splitapply( ...
    @(values) numel(unique(values)), results_table.Trial, group_id);
group_minimum_trial = splitapply(@min, results_table.Trial, group_id);
group_maximum_trial = splitapply(@max, results_table.Trial, group_id);

invalid_group = ...
    group_trial_count ~= expected_trials_per_condition | ...
    group_unique_trial_count ~= expected_trials_per_condition | ...
    group_minimum_trial ~= 1 | ...
    group_maximum_trial ~= expected_trials_per_condition;

if any(invalid_group)
    invalid_groups = table( ...
        group_subject(invalid_group), ...
        group_electrode(invalid_group), ...
        group_intensity(invalid_group), ...
        group_trial_count(invalid_group), ...
        group_unique_trial_count(invalid_group), ...
        group_minimum_trial(invalid_group), ...
        group_maximum_trial(invalid_group), ...
        'VariableNames', { ...
        'Subject', ...
        'Electrode', ...
        'Intensity_mT', ...
        'NumberOfRows', ...
        'NumberOfUniqueTrials', ...
        'MinimumTrial', ...
        'MaximumTrial'});

    disp(invalid_groups);
    error(['Every subject/electrode/intensity combination must contain ', ...
           'exactly trials 1 through %d.'], ...
        expected_trials_per_condition);
end

subject_electrode_pairs = unique( ...
    results_table(:, {'Subject', 'Electrode'}));
expected_number_of_groups = ...
    height(subject_electrode_pairs) * numel(intensities_to_process_mt);

if numel(group_trial_count) ~= expected_number_of_groups
    error(['At least one complete subject/electrode/intensity group is ', ...
           'missing from the processed dataset.']);
end

%% Save outputs

csv_output = fullfile(results_dir, 'gamma_power_trial_level.csv');
mat_output = fullfile(results_dir, 'gamma_power_trial_level.mat');

writetable(results_table, csv_output);

processing_parameters = struct();
processing_parameters.stimulation_frequency_hz = ...
    stimulation_frequency_hz;
processing_parameters.intensities_to_process_mt = ...
    intensities_to_process_mt;
processing_parameters.original_sampling_rate_hz = ...
    original_sampling_rate_hz;
processing_parameters.epoch_duration_s = epoch_duration_s;
processing_parameters.retained_interval_s = retained_interval_s;
processing_parameters.downsampling_factor = downsampling_factor;
processing_parameters.gamma_band_hz = gamma_band_hz;
processing_parameters.notch_frequencies_hz = notch_frequencies_hz;
processing_parameters.notch_half_width_hz = notch_half_width_hz;
processing_parameters.welch_window_duration_s = ...
    welch_window_duration_s;
processing_parameters.welch_overlap_fraction = welch_overlap_fraction;
processing_parameters.specparam_settings = specparam_settings;
processing_parameters.specparam_version = specparam_version;

save(mat_output, 'results_table', 'processing_parameters');

if save_frequency_resolved_spectra
    spectra_output = fullfile(results_dir, ...
        'gamma_spectra_trial_level.mat');

    spectra_table = table( ...
        Subject, ...
        Electrode, ...
        Intensity_mT, ...
        Trial, ...
        Frequency_Hz, ...
        RawSpectrum, ...
        SpecparamAdjustedSpectrum, ...
        LogLogAdjustedSpectrum, ...
        EpochFile);

    spectra_table = sortrows(spectra_table, ...
        {'Subject', 'Electrode', 'Intensity_mT', 'Trial'});

    save(spectra_output, 'spectra_table', 'processing_parameters', '-v7.3');
end

fprintf('[SAVE] Trial-level results: %s\n', csv_output);
fprintf('[SAVE] MATLAB results: %s\n', mat_output);

if save_frequency_resolved_spectra
    fprintf('[SAVE] Frequency-resolved spectra: %s\n', spectra_output);
end

fprintf('[DONE] Processing completed successfully.\n');
