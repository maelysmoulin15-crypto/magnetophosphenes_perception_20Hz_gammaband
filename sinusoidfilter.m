function res = sinusoidfilter(fname, sr, hi, ramp1, lo, ramp2)
%SINUSOIDFILTER FFT-domain band-pass filter.
%
%   RES = SINUSOIDFILTER(FNAME, SR, HI, RAMP1, LO, RAMP2) removes the
%   mean of FNAME, applies a high-pass cutoff at HI Hz and a low-pass
%   cutoff at LO Hz. RAMP1 and RAMP2 specify the transition widths in Hz.
%
%   With RAMP1 = RAMP2 = 0, the frequency response is rectangular:
%   frequencies <= HI are rejected and frequencies <= LO are retained by
%   the low-pass stage. For example:
%
%       y = sinusoidfilter(x, 10000, 30, 0, 80, 0);
%
%   applies the 30--80 Hz mask used in the EEG preprocessing pipeline.
%
%   The helper functions MAKEFILT_LOCAL and DATAFILT_LOCAL are private to
%   this file, so no additional filtering-function files are required.
%
%   Original filtering routines: Alexandre Legros, 2009.

validateattributes(fname, {'numeric'}, ...
    {'vector','real','finite','nonempty'}, mfilename, 'fname');
validateattributes(sr, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, 'sr');
validateattributes(hi, {'numeric'}, ...
    {'scalar','real','finite','nonnegative'}, mfilename, 'hi');
validateattributes(lo, {'numeric'}, ...
    {'scalar','real','finite','positive'}, mfilename, 'lo');
validateattributes(ramp1, {'numeric'}, ...
    {'scalar','real','finite','nonnegative'}, mfilename, 'ramp1');
validateattributes(ramp2, {'numeric'}, ...
    {'scalar','real','finite','nonnegative'}, mfilename, 'ramp2');

if hi >= lo
    error('sinusoidfilter:InvalidCutoffs', ...
        'The high-pass cutoff HI must be lower than the low-pass cutoff LO.');
end

fname = fname(:);
fname = fname - mean(fname);

highpass_filter = 1 - makefilt_local( ...
    hi / sr, length(fname), ramp1 / sr);
highpassed = datafilt_local(fname, highpass_filter);

lowpass_filter = makefilt_local( ...
    lo / sr, length(highpassed), ramp2 / sr);
res = datafilt_local(highpassed, lowpass_filter);

end


function yf = makefilt_local(cut, len, ramp)
%MAKEFILT_LOCAL Construct a symmetric low-pass mask in the DFT domain.

if cut >= 0.5
    error('sinusoidfilter:CutoffTooHigh', ...
        'Cutoff frequency must be lower than the Nyquist frequency.');
end

if cut + ramp >= 0.5
    error('sinusoidfilter:RampTooWide', ...
        'Cutoff frequency plus transition width must be below Nyquist.');
end

high1 = len / 2;
high2 = len / 2 - 1;

if floor(len / 2) ~= len / 2
    high2 = floor(len / 2);
end

q = ([0:high1, high2:-1:1] / len).';

if ramp ~= 0
    t = (cut + ramp - q) / ramp;
    t = t .* (t < 0.99999) .* (t > 0.0001);
else
    t = zeros(size(q));
end

yf = t + (q <= cut);

end


function y = datafilt_local(data, filt)
%DATAFILT_LOCAL Apply a real symmetric filter in the DFT domain.

data = data(:);
filt = filt(:);

if length(data) ~= length(filt)
    error('sinusoidfilter:LengthMismatch', ...
        'DATA and FILT must have the same length.');
end

d = fft(data);
d = d .* filt;
y = real(ifft(d));

end
