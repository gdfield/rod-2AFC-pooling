function ReturnComponents = ComponentFilter(Components, FreqCutoff)

%
%   function ReturnComponents = ComponentFilter(Components, FreqCutoff)
%
%   This function take a set of principal components or any waves form given
%   in "Compoents."  It makes a condition out of these and then feeds them
%   to ApplyFrequencyCutoff for filtering.  Created to filted principal
%   component that had a discontinuity introduced by time shifting
%
%   Created:  GDF:  11/15/04


% Express components as row vectors
Components = Components';

% create a condition structure so that "ApplyFrequencyCutoff can be used
NumComponents = length(Components(:,1));
TempCondition.EpochData.Data = Components;
TempCondition.EpochData.Offset(1:NumComponents) = 0;
TempCondition.SearchCrit{1} = 'StimDur';
TempCondition.SearchCrit{2} = 'SampInterv';
TempCondition.SearchPara(1) = 10;
TempCondition.SearchPara(2) = 1000;
TempCondition.DecimatePts = 1;
TempCondition.EpochNumbers = 1:NumComponents;
TempCondition.ExcludeEpochs(1:NumComponents) = 0;
TempCondition.UserInfo = [];


[junk, TempCondition] = ApplyFrequencyCutoff([], TempCondition, FreqCutoff);


ReturnComponents = TempCondition.EpochData.Data';

