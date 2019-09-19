function [fvecOneS, datvecPSD, phaseFFT] = CI_fft_analysis_e120(Tau, datvec, varargin)
%CI_FFT_ANALYSIS_E100 Does elementary Fourier analysis for the context of
% CI divergence of trajectory overlap simulation/data. Assumes that the
% data is equally spaced over the time span Tau (which needn't start at
% zero). Outputs 1 sided PSD and Phase.
%
% Version 110: 9/5/19: Allows you to specify along which direction to do
% the fft. 1 is down the columns, 2 is across the rows. Default is 1,
% unless input is just a vector
%
% Tau         : Time span of data = (last time) - (first time)
% datvec      : Squared overlap data input to be fourier transformed
%
% fvecOneS    : Frequency vector, one sided
% overSqPSD   : Power Spectral Density of overlapSq

[rows, cols] = size(datvec);

% Check if using input or default. For default, if a vec not mat use that
% dim, else for default mat use dim=1 to go down columns
if size(varargin) == 1
    dim = varargin{1};
else
    if rows == 1
        dim = 2;
    else
        dim = 1;
    end
end

switch dim
    case 1
        Nsamp     = rows;
        Nrept     = cols; 
        rowPsdIdx = 1:floor(Nsamp/2);           % Get PSD by going up to Nyquist. Here, time is down columns so lives as row
        colPsdIdx = 1:Nrept;                    % Grab all columns since time down rows

        rowTm2Idx = 2:(floor(Nsamp/2) - 1 + mod(Nsamp,2)); % Except for DC and Nyquist (present only for even lengths), must muliply by two (hence 'TiMes 2')
        colTm2Idx = 1:Nrept;
        
    case 2
        Nsamp = cols;
        Nrept = rows;         
        rowPsdIdx = 1:Nrept;
        colPsdIdx = 1:floor(Nsamp/2);

        rowTm2Idx = 1:Nrept;
        colTm2Idx = 2:(floor(Nsamp/2) - 1 + mod(Nsamp,2));

end

dt           = Tau/Nsamp;                        % Sample period
fs           = 1/dt;                             % Sample frequency
df           = fs/Nsamp;                         % Frequency step
fvecOneS     = df*(0:floor(Nsamp/2)-1);          % One sided frequency vector
fvecTwoS     = (-fs/2 : df : fs/2-df)' + mod(Nsamp,2)*df/2;  % Two sided freq vec

datvecFFT    = fft(datvec, [], dim);             % Take FFT
datvecFFTpos = datvecFFT(rowPsdIdx, colPsdIdx);  % Extract the positive frequencies up to Nyquist
datvecPSD    = (abs(datvecFFTpos/Nsamp)).^2;     % Get PSD: extract up to Nyquist, abs, square to get power
datvecPSD(rowTm2Idx, colTm2Idx) = 2*datvecPSD(rowTm2Idx, colTm2Idx); % Multiply by 2, except for first term (DC), or last term if even (Nyquist)

phaseFFT     = atan2(imag(datvecFFTpos), real(datvecFFTpos));        % Get the phase by arctan of the imaginary over real parts of the FFT

end

