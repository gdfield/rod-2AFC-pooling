function ContinuousNoise = GenerateContinuousNoise(NumTrials, PCAInfo, Verbose)

% hard coded value!!!!!!!!
SamplingInterval = 0.001;
NoisePS = sqrt(PCAInfo.NoisePS');
NoisePS = repmat(NoisePS, 1, NumTrials);

% generate continuous noise trials: gaussian white noise filters with square root of power
% spectrum of measured noise : each column is a noise trial NOT the row!!
ContinuousNoise = normrnd(0, 1, length(PCAInfo.MeanSingle), NumTrials);

ContinuousNoiseFFT = fft(ContinuousNoise);
ContinuousNoiseFFT = ContinuousNoiseFFT .* NoisePS;

ContinuousNoise = real(ifft(ContinuousNoiseFFT));
ContinuousNoise = ContinuousNoise';

% set std of continuous noise to measured value
ContinuousNoise = ContinuousNoise * PCAInfo.NoiseNormFactor;

% make sure mean of each trace is zero
ContinuousNoise = ContinuousNoise - repmat(mean(ContinuousNoise,2),1,size(ContinuousNoise,2)); 
