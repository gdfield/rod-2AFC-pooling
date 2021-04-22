function [AResponses, BResponses] = MakeRodResponses(PCAInfo, Shift, FlashStrength, NumRodsInPool, NoiseParameters, Verbose)

BaseSampleNumber = 2000;    % matlab memory issues preclude sampling more than 2000 responses at a time -- this is a work around
NumSubSims = NumRodsInPool ./ BaseSampleNumber;
RemainingRods = mod(NumRodsInPool, BaseSampleNumber);

for cnt = 1:NumSubSims

    [ConditionA, ConditionB] = RotatedPCSim(PCAInfo, Shift, BaseSampleNumber, FlashStrength, NoiseParameters, 0);

    % sum responses
    [ATempData, junk] = GetGoodEpochData([], ConditionA);
    [BTempData, junk] = GetGoodEpochData([], ConditionB);

    if cnt == 1
        ATemp = ATempData;
        BTemp = BTempData;
    else
        ATemp = [ATemp; ATempData];
        BTemp = [BTemp; BTempData];
    end

end

if RemainingRods ~= 0;
    [ConditionA, ConditionB] = RotatedPCSim(PCAInfo, Shift, RemainingRods, FlashStrength, NoiseParameters, 0);
   
    % sum responses
    [ATempData, junk] = GetGoodEpochData([], ConditionA);
    [BTempData, junk] = GetGoodEpochData([], ConditionB);
    
    if exist('ATemp')
        ATemp = [ATemp; ATempData];
        BTemp = [BTemp; BTempData]; 
    else
        ATemp = ATempData;
        BTemp = BTempData;
    end

end

AResponses = ATemp;
BResponses = BTemp;
        