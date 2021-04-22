function RFParams = Rod2AFCPooling(save_path, varargin)

% This is a wrapper function for doing 2AFC discrimination task
% over pools (receptive fields composed) of rods.

p =inputParser;
p.addParameter('noise_to_vary', 'thermals', @ischar);
p.addParameter('RF_size', 100, @isnumeric);
p.addParameter('discrim_method', 'threshold', @ischar);
p.addParameter('RawTimeShifts', [400, 200, 120, 80, 40, 20, 10], @isnumeric);
p.addParameter('FlashStrengths', [0.001 0.003 0.01 0.03 0.06 0.1 0.3 0.6 1.0 3], @isnumeric);
p.addParameter('NoiseFactors', [0.001, 0.01, 0.1 0.3 1 3], @isnumeric);
p.addParameter('TrainSize', 500, @isnumeric);
p.addParameter('TestSize', 400, @isnumeric);
p.addParameter('Verbose', false, @islogical);
p.parse(varargin{:});

Method = p.Results.discrim_method;
RFSize = p.Results.RF_size;
noise_to_vary = p.Results.noise_to_vary;
RawTimeShifts = p.Results.RawTimeShifts;
NumShifts = length(RawTimeShifts);
FlashStrengths = p.Results.FlashStrengths;
NoiseFactors = p.Results.NoiseFactors;
NumConditions = length(NoiseFactors);
TrainSize = p.Results.TrainSize;
TestSize = p.Results.TrainSize;
Verbose = p.Results.Verbose;


% Commented by GDF 2008-03-12
% The discrimination takes place at the level of individual rod responses.
% Discriminants are constructed from the time shifted responses and applied
% across flash strengths.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Load Stuff for simulating a pool of rod responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%load PCAInfo

test = load('PCAInfo');
PCAInfo = test.PCAInfo;
%test = load('CompleteConcatData');
%CompleteConcatData = test.CompleteConcatData;
test = load('IsoSingles');
IsoSingles = test.IsoSingles;
%SinglesCondition = CompleteConcatData.SinglesCondition;
test = load('NoiseCondition');
NoiseCondition = test.NoiseCondition;
test = load('SinglesAveResponse');
SinglesCondition.AverageResponse = test.SinglesAveResponse;
%NoiseCondition = CompleteConcatData.NoiseCondition;
clear CompleteConcatData
clear test
SinglesPCs = PCAInfo.SinglesPCs;
NumDim = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Rotate PCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The the conversion factor between the average interpolated time and real time
MinIndex = zeros(1, length(IsoSingles));
for cells = 1:length(IsoSingles)
    [~, TempMinIndex] = min(IsoSingles(cells).SinglesConditionA.AverageResponse);
    MinIndex(cells) = TempMinIndex;
end
AveMinIndex = mean(MinIndex);
[~, InterpolatedMin] = min(SinglesCondition.AverageResponse);

% actual time [in vector] times conversion factor 
TimeShifts = floor(RawTimeShifts .* (InterpolatedMin/AveMinIndex));    
[NoiseData, ~] = GetGoodEpochData([],NoiseCondition);
RotatedPCs = struct([]);
for shift = 1:NumShifts
	
	clear TempPCs MeanShiftedSingle
	TempPCs(TimeShifts(shift) + 1:2000,:) = SinglesPCs(1:2000 - TimeShifts(shift),:);
    MeanShiftedSingle(TimeShifts(shift) + 1:2000) = PCAInfo.MeanSingle(1:2000 - TimeShifts(shift));

    if Verbose
        figure(1)
        for cnt = 1:NumDim
            if cnt == 1
                subplot(NumDim, 1, 1)
                plot(MeanShiftedSingle, 'r')
                title('shifted single photon response')
            else
                subplot(NumDim, 1, cnt)
                plot(TempPCs(:, cnt))
            end
        end
        pause(1)
    end
	
    %Filter time shifted PCs to avoid discontinuities
    TempPCs = ComponentFilter(TempPCs, length(IsoSingles));
    MeanShiftedSingle = ComponentFilter(MeanShiftedSingle', length(IsoSingles));
    
    % check to make sure the basic shape of the components has not been
    % changed by the filtering -- sanity check on filtering.
    if Verbose
        figure(2)
		for cnt = 1:NumDim
            if cnt == 1
                subplot(NumDim, 1, 1)
                plot(MeanShiftedSingle, 'r')
                title('shifted single photon response post filtering')
            else
                subplot(NumDim, 1, cnt)
                plot(TempPCs(:, cnt))
            end
		end
		pause(1)
    end

	RotatedPCs(shift).PCs = TempPCs;
    RotatedPCs(shift).MeanSingle = MeanShiftedSingle';
    ContNoise_ShiftedSpace = NoiseData * TempPCs;
    RotatedPCs(shift).ContNoiseCov = cov(ContNoise_ShiftedSpace);
    
    % Compute and store an explicit difference of means discriminator
    RotatedPCs(shift).DifferenceOfMeans = PCAInfo.MeanSingle - MeanShiftedSingle';
    BiasThreshold = ((PCAInfo.MeanSingle * PCAInfo.MeanSingle') - (MeanShiftedSingle' * MeanShiftedSingle)) .* 0.5;
    RotatedPCs(shift).BiasThreshold = BiasThreshold;
end

PCAInfo.RotatedPCs = RotatedPCs;
PCAInfo.TimeShifts = TimeShifts;
clear IsoSingles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% linear pooling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NoiseParameters.ContNoiseScale = 1;
NoiseParameters.SinglesVarScale = 1;
NoiseParameters.AmpVar = 0;
NoiseParameters.PoissonFlag = 1.0;
NoiseParameters.ThermalRate = 0.0035;

% Set Discrimination Parameters
DiscriminationParameters.DimensionNumber = 2*NumDim;
DiscriminationParameters.FisherTrainSetSize = TrainSize;
DiscriminationParameters.TrainingFlashStrength = 1.0;
DiscriminationParameters.PoissonFlag = 0;
DiscriminationParameters.TestSetSize = TestSize;   % number of test responses (trials) for a class

% Set simulated Response Parameters
StimulationParameters.TimeShifts = TimeShifts;
StimulationParameters.FlashStrengths = FlashStrengths;
StimulationParameters.NumShifts = length(TimeShifts);   
StimulationParameters.NumFlashes = length(FlashStrengths);


Verbose = 0;
cd(save_path)
RFParams = struct([]);
for RF = 1:NumConditions
    fprintf('Number of rods in receptive field simulation %d \n', RFSize);
    
    switch noise_to_vary
        case 'thermals'; NoiseParameters.ThermalRate = 0.0035 * NoiseFactors(RF);
            
        case 'continuous'; NoiseParameters.ContNoiseScale = NoiseFactors(RF);

        case 'singles'; NoiseParameters.SinglesVarScale = NoiseFactors(RF).^2;

        case 'all'
            NoiseParameters.SinglesVarScale = NoiseFactors(RF).^2;
            NoiseParameters.ContNoiseScale = NoiseFactors(RF);      
            NoiseParameters.ThermalRate = 0.0035 * NoiseFactors(RF);
            
        case 'amp'
            NoiseParameters.AmpVar = NoiseFactors(RF);            
    end

    
    StimulationParameters.RodNumber = RFSize;
    TempPCorrect = RodDiscriminationAndPooling(PCAInfo, StimulationParameters, DiscriminationParameters, NoiseParameters, Method, Verbose);
    RFParams(RF).PCorrect = TempPCorrect;
    RFParams(RF).RFSize = StimulationParameters.RodNumber;

    DetectThresh = zeros(1,StimulationParameters.NumShifts);
    for shift = 1:StimulationParameters.NumShifts
        coef = 0.005;
        fitcoef = nlinfit(FlashStrengths', TempPCorrect(shift, :)', 'cumulative_gaussian', coef);
		FitCorrect = cumulative_gaussian(fitcoef, FlashStrengths');
        figure(32)
		semilogx(FlashStrengths, FitCorrect, 'g', FlashStrengths, TempPCorrect(shift, :), 'k')
        DetectThresh(shift) = norminv(0.75, 0, abs(fitcoef));
        figure(33)
        semilogx(NoiseFactors, DetectThresh(shift), 'ko')
        hold on
    
    end
    RFParams(RF).DetectThresh = DetectThresh;
    RFParams(RF).NoiseParameters = NoiseParameters;
    RFParams(RF).DiscriminationParameters = DiscriminationParameters;
    RFParams(RF).StimulationParameters = StimulationParameters;
    
    save RFParams RFParams
end

