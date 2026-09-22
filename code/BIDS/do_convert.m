% the original dataset includes data from 40 subjects
% we are just using a subset here to demonstrate the principle
P = setup_paths();   % FieldTrip location comes from setup_paths_config.m
filenames = {
  'S1001_run01.nirs'
  'S1003_run01.nirs'
  'S1004_run01.nirs'
  'S1005_run01.nirs'
  };


%%

for i=1:length(filenames)
  filename = filenames{i};
  
  cfg = [];
  cfg.dataset = fullfile('original', filename);
  % cfg.dataset = filename;
  cfg.method = 'convert';
  
  % the following settings relate to the directory structure and file names
  cfg.bidsroot = 'bids';
  cfg.sub = filename(1:5);
  cfg.ses = [];
  cfg.run = [];
  cfg.task = 'listenandrepeat';
  cfg.datatype = 'nirs';
  
  % the following settings relate to the dataset_description.json
  cfg.dataset_description.Name                = 'Defenderfer 2019; fNIRS data files for event-related vocoding/background noise study';
  cfg.dataset_description.Authors             = 'Defenderfer, Jessica; Buss, Aaron ';
  cfg.dataset_description.DatasetDOI          = 'http://dx.doi.org/10.17632/4cjgvyg5p2.1';
  cfg.dataset_description.License             = 'CC BY 4.0';
  cfg.dataset_description.ReferencesAndLinks  = {'http://www.fieldtriptoolbox.org/example/nirs_speech/'}; % this can be a list
  cfg.dataset_description.BIDSVersion         = 'BEP030'; % this does not correspnd to an official version, but a BIDS Extension Proposal. See http://bids.neuroimaging.io/bep030
  
  data2bids(cfg)
end

