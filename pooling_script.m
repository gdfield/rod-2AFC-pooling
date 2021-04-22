
% add path to the code -- will need to change path
addpath(genpath('~/Desktop/rod-2afc-pooling-main/'))

clear sim_params
sim_params.FlashStrengths = [0.003 0.006 0.01 0.03 0.06 0.1 0.3 0.6 1]; % units are Rh*/rod
sim_params.RawTimeShifts = [400 20]; % list of time shifts between flashes in 2AFC task
sim_params.NoiseFactors = [0.1 1]; % scale factors on varying noise parameters
sim_params.RF_size = 100; % rod pool size
% set noise to vary
%   options are: 'thermals
%                'continuous'
%                'singles'
%                'all'
%                'amp'
sim_params.noise_to_vary = 'thermals'; 
% set discrimination method
%   options are: 'difference of means'
%                'threshold'
sim_params.discrim_method = 'threshold';
sim_params.Verbose = false; % plot/print more stuff
sim_params.TrainSize = 500; % number of responses to simulate for training
sim_params.TestSize = 400; % number of resposnes to simulate for testing


RFParams = Rod2AFCPooling('~/Desktop/', sim_params);
