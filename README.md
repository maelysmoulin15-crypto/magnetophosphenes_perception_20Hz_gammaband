# EEG gamma-band processing during 20 Hz magnetic-field exposure

This repository provides the MATLAB code used to transform individual EEG
epochs into trial-level gamma-band spectral measures for a study of
magnetophosphene perception during 20 Hz magnetic-field exposure.

The repository currently focuses on the central spectral-processing stage of
the analysis. It implements the complete workflow from a pre-extracted,
average-referenced EEG epoch to raw and aperiodic-adjusted gamma-power
estimates.

- `processing.m` runs the trial-level spectral-processing pipeline and saves
  the resulting gamma-power measures.
- `sinusoidfilter.m` applies the real and symmetric DFT-domain band-pass mask.
  Its two helper routines are included as private subfunctions in the same
  file.
- `data/epochs/` is the default location for the input EEG epochs. EEG data are
  not distributed with this repository.
- `external/specparam_mat/` is the default location for the MATLAB wrapper used
  to call specparam.
- `results/processing/` is created automatically when the script is run.

## Scope of the processing pipeline

For every electrode and trial, `processing.m` performs the following steps:

1. Verifies that the input epoch contains exactly 5 s of finite EEG data
   sampled at 10 kHz.
2. Multiplies the complete epoch by a symmetric Hann window.
3. Applies a real and symmetric rectangular DFT-domain mask between 30 and
   80 Hz.
4. Retains the interval from 1 to 4 s of the original epoch.
5. Applies second-order IIR notch filters at 40, 50, and 60 Hz using
   zero-phase `filtfilt` filtering.
6. Downsamples the signal by a factor of five, from 10 kHz to 2 kHz.
7. Estimates the power spectral density using Welch's method.
8. Calculates raw 30--80 Hz gamma power.
9. Calculates gamma power after local aperiodic adjustment with specparam.
10. Calculates gamma power after complementary log--log aperiodic adjustment.

All three integrated power measures exclude Welch bins located within
plus or minus 1.5 Hz of 40, 50, and 60 Hz. With the 1-Hz Welch grid, 42 of the
initial 51 bins between 30 and 80 Hz are retained for numerical integration.

## Main analysis parameters

| Parameter | Value |
|---|---:|
| Stimulation frequency | 20 Hz |
| Conditions processed by default | 0, 5, and 50 mT |
| Original sampling rate | 10 kHz |
| Input epoch duration | 5 s |
| Retained interval | 1--4 s |
| DFT-domain band-pass mask | 30--80 Hz |
| Notch-filter centers | 40, 50, and 60 Hz |
| Nominal notch bandwidth | 3 Hz total per notch |
| Effective sampling rate | 2 kHz |
| Welch window | 1-s periodic Hann window |
| Welch overlap | 50% |
| Welch frequency resolution | 1 Hz |
| Welch segments per trial | 5 |

The DFT-domain mask rejects frequencies at or below 30 Hz, retains frequencies
strictly above 30 Hz through 80 Hz, and rejects frequencies above 80 Hz. For a
5-s epoch sampled at 10 kHz, the DFT spacing is 0.2 Hz; therefore, the first
retained positive-frequency bin is 30.2 Hz and the 80-Hz bin is retained.

## Specparam settings

The script requires specparam 2.x and uses the following local parameterization
settings over the restricted 30--80 Hz spectrum:

| Setting | Value |
|---|---:|
| Peak-width limits | 1--12 Hz |
| Maximum number of peaks | 6 |
| Minimum peak height | 0 |
| Peak threshold | 2 SD |
| Aperiodic mode | Fixed |

The complete 30--80 Hz Welch grid, including the bins attenuated by temporal
notch filtering, is supplied to specparam. The fitted aperiodic background is
converted back to linear power units and subtracted from the observed PSD.
Negative residual values are set to zero before integration.

The fitted exponent, offset, and goodness of fit describe the local 30--80 Hz
spectrum and should not be interpreted as estimates of the global aperiodic
spectral structure.

## Log--log adjustment

The complementary log--log model is fitted as

```text
log10(power) = intercept + slope * log10(frequency)
```

Bins within plus or minus 1.5 Hz of 40, 50, and 60 Hz are excluded before this
regression. The fitted background is reconstructed in linear power units and
subtracted from the observed PSD. Negative residual values are set to zero
before integration over the same 42 frequency bins used for the raw and
specparam-adjusted measures.

## Requirements

- MATLAB with the Signal Processing Toolbox
- Python environment containing specparam 2.x
- MATLAB wrapper providing the `specparam` function used in `processing.m`
- `sinusoidfilter.m`, included in this repository

The analysis was designed for specparam 2.x. The script checks the installed
Python package version before processing the data.

## Python configuration

The script uses MATLAB's currently configured Python environment by default.
Alternatively, define the `SPECPARAM_PYTHON` environment variable with the path
to the required Python executable before running `processing.m`.

For example, from MATLAB:

```matlab
setenv('SPECPARAM_PYTHON', '/path/to/python');
```

The path is intentionally not hard-coded so that the repository remains
portable across systems.

## Input data

Each input MAT file must contain:

- `EEG_segment`: one 5-s, average-referenced, single-electrode EEG epoch;
- either `meta.Fs` or `time_segment`, allowing the sampling rate to be read or
  calculated.

The default filename convention is:

```text
Subject<subject>_<configuration>_20Hz_<electrode>_<intensity>mT_trial<trial>.mat
```

For example:

```text
Subject07_OCCIPITAL_20Hz_Oz_50mT_trial3.mat
```

The prefix `Sujet` is also accepted for compatibility with the original data
files. By default, the script retains the 0, 5, and 50 mT conditions and expects
five trials for every subject, electrode, and condition.

The epochs must already have been extracted and average-referenced. Stimulation
onset detection, epoch reconstruction, channel referencing, and any other
operations performed on the continuous recordings are outside the scope of
`processing.m`.

## Running the analysis

1. Place the epoch files in `data/epochs/`.
2. Place the MATLAB specparam wrapper in `external/specparam_mat/`, or add it
   to the MATLAB path.
3. Ensure that `sinusoidfilter.m` is located beside `processing.m`.
4. Configure the Python environment if necessary.
5. In MATLAB, change the current directory to the repository root and run:

```matlab
processing
```

The script stops with an explicit error if an epoch has an unexpected sampling
rate or length, contains non-finite values, does not produce the expected
frequency grid, or if a subject/electrode/condition combination does not
contain trials 1 through 5.

## Outputs

The following files are written to `results/processing/`:

### `gamma_power_trial_level.csv`

One row per subject, electrode, condition, and trial, including:

- raw integrated gamma power;
- specparam-adjusted integrated gamma power;
- log--log-adjusted integrated gamma power;
- local specparam aperiodic parameters and goodness of fit;
- local log--log slope and offset parameters;
- sampling-rate information and source filename.

### `gamma_power_trial_level.mat`

MATLAB version of the trial-level results, together with a structure containing
the processing parameters and detected specparam version.

### `gamma_spectra_trial_level.mat`

Frequency-resolved raw, specparam-adjusted, and log--log-adjusted spectra over
the 42 retained bins. This output can be disabled by setting:

```matlab
save_frequency_resolved_spectra = false;
```

## Reproducibility notes

- Welch segments are averaged within each 3-s stimulation epoch and are not
  treated as independent observations.
- Numerical integration is implemented as the sum of the retained PSD bins
  multiplied by the frequency resolution. This avoids bridging the gaps
  created by logical exclusion of the notch-region bins.
- No z-score normalization is applied; signal-amplitude information is
  preserved.
- The script does not silently remove incomplete or invalid epochs. It stops
  and reports the affected file.
- Electrode exclusion, subject-level quality control, aggregation, statistical
  modeling, and multiple-comparison correction are not performed by this
  processing script.

## Data availability

The participant-level EEG data are not distributed with this repository.
External dissemination of these data requires the relevant contractual and
institutional authorizations. Aggregate results supporting the associated
publication are reported in the article and its supplementary materials.

## Citation

If you use this processing pipeline, please cite the associated article:

> *Broadband gamma-band EEG changes during magnetophosphene-inducing 20 Hz magnetic field stimulation*
>
> Preprint: [https://doi.org/10.64898/2026.04.15.718626](https://doi.org/10.64898/2026.04.15.718626)

