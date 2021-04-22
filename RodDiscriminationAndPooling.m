function PCorrect = RodDiscriminationAndPooling(PCAInfo, StimulationParameters, DiscriminationParameters, NoiseParameters, Method, Verbose)
%
%   function PCorrect = RodDiscriminationAndPooling(PCAInfo, ResponseParameters, DiscriminationParameters,...
%                                    NoiseParameters, Verbose)
%
%   Discriminate rod responses at the level of individual responses and add
%   the projections to form a "pooled" response.
%   Method = 'differene of means'
%          = 'threshold'
%
%   Created GDF 2008-3


NumRodsInPool = StimulationParameters.RodNumber;
NumTestTrials = DiscriminationParameters.TestSetSize;

% Initialize large matrices
% generate and classify test data sets
AResponseClass = zeros(NumTestTrials, 1);
BResponseClass = zeros(NumTestTrials, 1);

PCorrect = zeros(StimulationParameters.NumShifts, StimulationParameters.NumFlashes);

status_factor = NumTestTrials ./ 20;

for shift = 1:StimulationParameters.NumShifts
    
    
    %% get Discriminat and Bias
    
    % difference of means discrimination
    if strcmp(Method, 'difference of means') || strcmp(Method, 'threshold')
        Discriminant = PCAInfo.RotatedPCs(shift).DifferenceOfMeans;
        Discriminant = Discriminant - mean(Discriminant);
        Discriminant = Discriminant ./ norm(Discriminant);
        
        BiasThreshold = PCAInfo.RotatedPCs(shift).BiasThreshold;

        %NLStats = EstimateNLParameters(PCAInfo, NoiseParameters, Discriminant', shift, DiscriminationParameters.FisherTrainSetSize)
        
        % difference of means discrimination
        
        if 1
            
            temp_filt = PCAInfo.MeanSingle;
            norm_temp_filt = temp_filt ./ norm(temp_filt);
            test_filt = filter(norm_temp_filt, 1,  temp_filt);
            tmp_max = max(test_filt);
            norm_temp_filt = norm_temp_filt ./ tmp_max;
            
            %eliminate poisson fluctuations
            test_NoiseParameters = NoiseParameters;
            test_NoiseParameters.PoissonFlag = 0;
            %eliminate thermal noise
            test_NoiseParameters.ThermalRate = 0;

            %simulate failures and estimate statistics
            test_FlashStrength = 0.0;
            [ConditionA, ConditionB] = RotatedPCSim(PCAInfo, shift, DiscriminationParameters.FisherTrainSetSize, test_FlashStrength, test_NoiseParameters, 0);
            NoiseData = [ConditionA.EpochData.Data; ConditionB.EpochData.Data];

            filtered_noise = filter(norm_temp_filt, 1,  NoiseData, [], 2);
            max_noise = max(filtered_noise, [], 2);
            
            
            %simulate singles
            test_FlashStrength = 1.0;
            [ConditionA, ConditionB] = RotatedPCSim(PCAInfo, shift, DiscriminationParameters.FisherTrainSetSize, test_FlashStrength, test_NoiseParameters, 0);
            ASinglesData = ConditionA.EpochData.Data;
            BSinglesData = ConditionB.EpochData.Data;

            temp_filtA = filter(norm_temp_filt, 1,  ASinglesData, [], 2);
            temp_filtB = filter(norm_temp_filt, 1,  BSinglesData, [], 2);
            temp_Amax = max(temp_filtA, [],2);
            temp_Bmax = max(temp_filtB, [],2);
            
            
            if Verbose
                figure(16); clf; hold on
                [hist_noise, hist_bins] = hist(max_noise, 50);
                plot(hist_bins, hist_noise, 'k')
                [hist_A, hist_bins_A] = hist(temp_Amax, 50);
                plot(hist_bins_A, hist_A, 'r')
                [hist_B, hist_bins_B] = hist(temp_Bmax, 50);
                plot(hist_bins_B, hist_B, 'b')
            end
            
            mean_noise = mean(max_noise);
            sd_noise = std(max_noise);
            mean_single = mean([temp_Amax; temp_Bmax]);
            sd_single = std([temp_Amax; temp_Bmax]);
                       
%             
        end
            
        
        
    elseif strcmp(Method, 'difference of means nonlinear')
        Discriminant = PCAInfo.RotatedPCs(shift).DifferenceOfMeans;
        Discriminant = Discriminant - mean(Discriminant);
        Discriminant = Discriminant ./ norm(Discriminant);
        
        %%
        Discriminant = Discriminant - mean(Discriminant);
        Discriminant = Discriminant ./ norm(Discriminant);

    end
       
    %%
    if Verbose
        % plot discriminant
        figure(1)
        if strcmp(Method,'difference of means') || strcmp(Method,'difference of means nonlinear') || strcmp(Method,'threshold')
            plot(Discriminant)
        elseif (strcmp(Method, 'Fisher') || strcmp(Method, 'FisherNonLinear')) 
            plot(FullDDiscriminant)            
        else
            error('Method must be Fisher or difference of means')
        end
        title('Discriminant')
        xlabel('time')
        ylabel('weight')
        drawnow
    end

    
    for flash = 1:StimulationParameters.NumFlashes
        
        FlashStrength = StimulationParameters.FlashStrengths(flash);
        
        fprintf('Simulating test responses \n')


            %h = waitbar(0, 'Simulating test discrimination responses');
        
        singles_prob = poisspdf(1,FlashStrength);
        failures_prob = poisspdf(0,FlashStrength);
        
        threshold_function = @(x)(failures_prob/sqrt(2.*pi.*mean_noise.^2)).*...  normaization for first term
        exp(-((x-mean_noise).^2) ./ (2*sd_noise.^2))-... gaussian
        (singles_prob/sqrt(2.*pi.*mean_single.^2)).*... normalization for 2nd term
        exp(-((x-mean_single).^2) ./ (2*sd_single.^2)); % gaussian        

        response_threshold = fzero(threshold_function, mean_single+sd_single);
        if isnan(response_threshold) 
            response_threshold = fzero(threshold_function, mean_noise + 2*sd_noise);
        end
        if isnan(response_threshold)|| response_threshold < 0
            response_threshold = mean_single - (2*sd_single);
        end
        
        %plot the fits to the distribuion 
        if Verbose
            figure(21); clf;
            temp_x = [(mean_noise-5*sd_noise):0.0002:1];
            temp_y = abs((failures_prob/sqrt(2.*pi.*sd_noise^2)).*...  normaization for first term
                    exp(-((temp_x-mean_noise).^2) ./ (2*sd_noise.^2))-... gaussian
                    (singles_prob/sqrt(2.*pi.*sd_single^2)).*... normalization for 2nd term
                    exp(-((temp_x-mean_single).^2) ./ (2.*sd_single.^2))); % gaussian
            semilogy(temp_x, temp_y, 'b')
            temp_y_a = (failures_prob/sqrt(2.*pi.*sd_noise^2)).*...  normaization for first term
                    exp(-((temp_x-mean_noise).^2) ./ (2*sd_noise.^2));

            temp_y_b =(singles_prob/sqrt(2.*pi.*sd_single^2)).*... normalization for 2nd term
                    exp(-((temp_x-mean_single).^2) ./ (2.*sd_single.^2)); % gaussian        

            hold on    
            semilogy(temp_x, temp_y_a, 'k', temp_x, temp_y_b, 'r')
            drawnow                 
        end
                   
       
        
       
        for test = 1:NumTestTrials
        
            if mod(test, status_factor) == 0
                fprintf('*')
            end
            if mod(test, NumTestTrials) == 0
                fprintf('\n')
            end

            [ATempResps, BTempResps] = MakeRodResponses(PCAInfo, shift, FlashStrength, NumRodsInPool, NoiseParameters, Verbose);
    
           
            % Discrimination at the level of individual responses and then
            % sum
            if strcmp(Method, 'difference of means')
                AProjections = ATempResps * Discriminant';
                BProjections = BTempResps * Discriminant';
                AProjections = AProjections - ((BiasThreshold * poisspdf(1, FlashStrength)) + (2 * BiasThreshold * poisspdf(2, FlashStrength)));
                BProjections = BProjections - ((BiasThreshold * poisspdf(1, FlashStrength)) + (2 * BiasThreshold * poisspdf(2, FlashStrength)));
            elseif strcmp(Method, 'difference of means nonlinear')
                AProjections = ATempResps * Discriminant';
                BProjections = BTempResps * Discriminant';
                

                % apply nonlinear weighting
                AResp_ASinglesProb = singles_prob * normpdf(AProjections, NLStats.ASinglesMean, NLStats.ASinglesCov);
                AResp_BSinglesProb = singles_prob * normpdf(AProjections, NLStats.BSinglesMean, NLStats.BSinglesCov);
                BResp_BSinglesProb = singles_prob * normpdf(BProjections, NLStats.BSinglesMean, NLStats.BSinglesCov);
                BResp_ASinglesProb = singles_prob * normpdf(BProjections, NLStats.ASinglesMean, NLStats.ASinglesCov);
                AResp_NoiseProb = failures_prob * normpdf(AProjections, 0, NLStats.ContNoiseCov);
                BResp_NoiseProb = failures_prob * normpdf(BProjections, 0, NLStats.ContNoiseCov);
                AProjections = AProjections - ((BiasThreshold * poisspdf(1, FlashStrength)) + (2 * BiasThreshold * poisspdf(2, FlashStrength)));
                BProjections = BProjections - ((BiasThreshold * poisspdf(1, FlashStrength)) + (2 * BiasThreshold * poisspdf(2, FlashStrength)));

                %%
                RespAA_likelihood = AResp_ASinglesProb ./ AResp_NoiseProb;
                RespAB_likelihood = AResp_BSinglesProb ./ AResp_NoiseProb;
                RespBB_likelihood = BResp_BSinglesProb ./ BResp_NoiseProb;
                RespBA_likelihood = BResp_ASinglesProb ./ BResp_NoiseProb;
                
                RespA_likelihood = max([RespAA_likelihood'; RespAB_likelihood']);
                RespB_likelihood = max([RespBB_likelihood'; RespBA_likelihood']);
                
                RespA_likelihood(RespA_likelihood<1) = 0;
                RespB_likelihood(RespB_likelihood<1) = 0;
                
                AProjections = AProjections .* RespA_likelihood';
                BProjections = BProjections .* RespB_likelihood';
                
           
            elseif strcmp(Method, 'threshold')
                
                    
                temp_filt = PCAInfo.MeanSingle;
                norm_temp_filt = temp_filt ./ norm(temp_filt);
                test_filt = filter(norm_temp_filt, 1,  temp_filt);
                tmp_max = max(test_filt);
                norm_temp_filt = norm_temp_filt ./ tmp_max;
                
                
                filtA = filter(norm_temp_filt, 1,  ATempResps, [], 2);
                filtB = filter(norm_temp_filt, 1,  BTempResps, [], 2);
                Amax = max(filtA, [],2);
                Bmax = max(filtB, [],2);
                
                ATempResps(Amax < response_threshold,:) = 0;
                BTempResps(Bmax < response_threshold,:) = 0;
                                
                AProjections = ATempResps * Discriminant';
                BProjections = BTempResps * Discriminant';
                

            else
                error('A supported method was not specified')
            end
                
            if Verbose == 1 && test == 0
                figure(3)
                AMaxDist = max(AProjections);
                AMinDist = min(AProjections);
                BMaxDist = max(BProjections);
                BMinDist = min(BProjections);
                MaxDist = max([AMaxDist, BMaxDist]);
                MinDist = min([AMinDist, BMinDist]);
                
               
                Bins = MinDist:((MaxDist - MinDist)/25):MaxDist;
                
                AHist = hist(AProjections, Bins);
                BHist = hist(BProjections, Bins);
                
                stairs(Bins, AHist, 'k')
                hold on
                stairs(Bins, BHist, 'r')
                title('single trial projections')
                hold off
                drawnow
            end

            AResponseClass(test) = sum(AProjections);
            BResponseClass(test) = sum(BProjections);

        end


        % compute the crossing point between these distributions

        all_noise_A = length(find(AResponseClass == 0));
        all_noise_B = length(find(BResponseClass == 0));
        fprintf('there were %g trials with no SPRs that passed threshold \n', (all_noise_A + all_noise_B))

        if length(find(AResponseClass ~=0)) > 10 && length(find(BResponseClass ~= 0)) > 10
            
            Abar = mean(AResponseClass(AResponseClass ~= 0));
            Bbar = mean(BResponseClass(BResponseClass ~= 0));
            sigA = std(AResponseClass(AResponseClass ~= 0));
            sigB = std(BResponseClass(BResponseClass ~= 0));

            myfun = @(x)(1/sqrt(2.*pi.*sigA.^2)).*...  normaization for first term
            exp(-((x-Abar).^2) ./ (2*sigA.^2))-... gaussian
            (1/sqrt(2.*pi.*sigB.^2)).*... normalization for 2nd term
            exp(-((x-Bbar).^2) ./ (2*sigB.^2)); % gaussian        

            DiscrimThreshold = fzero(myfun, mean([Abar, Bbar]));
            
            if DiscrimThreshold > Abar || DiscrimThreshold < Bbar
                disp('Discrimination threshold out of bounds, setting to midpoint between means')
                DiscriminThreshold = mean([Abar, Bbar]);
            end

            %plot the fits to the distribuion 
            if Verbose
                figure(18); clf;
                temp_x = [(Bbar-3*sigB):0.05:(Abar+3*sigA)];
                temp_y = abs((1/sqrt(2.*pi.*sigA^2)).*...  normaization for first term
                        exp(-((temp_x-Abar).^2) ./ (2*sigA.^2))-... gaussian
                        (1/sqrt(2.*pi.*sigB^2)).*... normalization for 2nd term
                        exp(-((temp_x-Bbar).^2) ./ (2.*sigB.^2))); % gaussian
                plot(temp_x, temp_y)
                temp_y_a = (1/sqrt(2.*pi.*sigA^2)).*...  normaization for first term
                        exp(-((temp_x-Abar).^2) ./ (2*sigA.^2));

                temp_y_b =(1/sqrt(2.*pi.*sigB^2)).*... normalization for 2nd term
                        exp(-((temp_x-Bbar).^2) ./ (2.*sigB.^2)); % gaussian        

                hold on    
                plot(temp_x, temp_y_a, 'k', temp_x, temp_y_b, 'r')
                drawnow  
            end
        else
            DiscrimThreshold = 0;
        end

        
        % compute probability correct
        APCorrect = (length(find(AResponseClass(AResponseClass ~= 0) >= DiscrimThreshold)) + floor(0.5* length(find(AResponseClass == 0)))) ./ NumTestTrials;
        BPCorrect = (length(find(BResponseClass(BResponseClass ~= 0) < DiscrimThreshold)) + ceil(0.5* length(find(BResponseClass == 0)))) ./ NumTestTrials;
        TempPCorrect = mean([APCorrect, BPCorrect]);
        
        if Verbose
            figure(2)
            % get extremes of distributions for plotting purposes
            AMaxDist = max(AResponseClass);
            AMinDist = min(AResponseClass);
            BMaxDist = max(BResponseClass);
            BMinDist = min(BResponseClass);
            MaxDist = max([AMaxDist, BMaxDist]);
            MinDist = min([AMinDist, BMinDist]);
            Bins = MinDist:((MaxDist - MinDist)/25):MaxDist;
            AHist = hist(AResponseClass(AResponseClass ~= 0), Bins);
            BHist = hist(BResponseClass(BResponseClass ~= 0), Bins);

            stairs(Bins, AHist, 'k')
            hold on
            stairs(Bins, BHist, 'r')
            plot([DiscrimThreshold DiscrimThreshold], [0 max([max(AHist), max(BHist)])], 'b')
            legend('classA','class B')
            title('Histogram of class A (black) and class B (red) responses')
            hold off
            drawnow
            
        end

        fprintf('the noise values were thermals: %g, cont: %g, singles: %g \n',NoiseParameters.ThermalRate, NoiseParameters.ContNoiseScale, NoiseParameters.SinglesVarScale)
        fprintf('For flash %g and time shift %g \n', FlashStrength, PCAInfo.TimeShifts(shift))
        fprintf('the probability correct for class A is %g \n', APCorrect)
        fprintf('the probability correct for class B is %g \n', BPCorrect)
        fprintf('The mean probability correct is %g \n \n', TempPCorrect)
        PCorrect(shift, flash) = TempPCorrect;
    end
end
            