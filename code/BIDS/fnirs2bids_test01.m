% fnirs2bids_test01
P = setup_paths();  % FieldTrip location comes from setup_paths_config.m
filenames = {
'01_task-resting01_nirs.snirf'
'02_task-resting01_nirs.snirf'
'03_task-resting01_nirs.snirf'
'04_task-resting01_nirs.snirf'
'05_task-resting01_nirs.snirf'
'06_task-resting01_nirs.snirf'
'07_task-resting01_nirs.snirf'
'08_task-resting01_nirs.snirf'
'09_task-resting01_nirs.snirf'
'10_task-resting01_nirs.snirf'
'11_task-resting01_nirs.snirf'
'12_task-resting01_nirs.snirf'
'13_task-resting01_nirs.snirf'
'14_task-resting01_nirs.snirf'
'15_task-resting01_nirs.snirf'
'16_task-resting01_nirs.snirf'
'17_task-resting01_nirs.snirf'
'18_task-resting01_nirs.snirf'
'19_task-resting01_nirs.snirf'
'20_task-resting01_nirs.snirf'
'21_task-resting01_nirs.snirf'
'22_task-resting01_nirs.snirf'
'23_task-resting01_nirs.snirf'
'24_task-resting01_nirs.snirf'
'25_task-resting01_nirs.snirf'
'26_task-resting01_nirs.snirf'
'27_task-resting01_nirs.snirf'
'28_task-resting01_nirs.snirf'
'29_task-resting01_nirs.snirf'
'30_task-resting01_nirs.snirf'
  };

%% .snirf to BIDS conversion

for i=1:length(filenames)
  filename = filenames{i};
  
  cfg = [];
  cfg.dataset = fullfile('..\data\original', filename);
  % cfg.dataset = filename;
  cfg.method = 'convert';
  
  % the following settings relate to the directory structure and file names
  cfg.task                                  = 'resting01';
  cfg.bidsroot                              = '..\data\resting01';
  cfg.sub                                   = filename(1:2);
  cfg.ses                                   = [];
  cfg.run                                   = [];
  cfg.suffix                                = 'nirs';
  
  % the following settings relate to the dataset_description.json
  cfg.dataset_description.Name                = 'Guevara 2025; fNIRS data files for Resident physician mental workload assessment';
  cfg.dataset_description.Authors             = 'Guevara, Edgar; Torres-Cuevas, Gerardo-Enrique; Martínez-Jiménez, Mario Aurelio';
  cfg.dataset_description.DatasetDOI          = 'http://dx.doi.org/10.5281/zenodo.15186570';
  cfg.dataset_description.License             = 'CC BY 4.0';
  cfg.dataset_description.ReferencesAndLinks  = {'http://dx.doi.org/10.5281/zenodo.15186570'}; % this can be a list
  cfg.dataset_description.BIDSVersion         = 'BEP030'; % this does not correspnd to an official version, but a BIDS Extension Proposal. See http://bids.neuroimaging.io/bep030

  % Various
  cfg.InstitutionName                       = 'Universidad Autonoma de San Luis Potosi';
  cfg.InstitutionAddress                    = 'Av Chapultepec 1570. Priv. del Pedregal, San Luis Potosi 78294, Mexico';
  cfg.InstitutionalDepartmentName           = 'Faculty of Science';
  cfg.Manufacturer                          = 'Artinis Medical Systems';
  cfg.ManufacturersModelName                = 'Brite MKII';
  cfg.DeviceSerialNumber                    = '24231';
  
  data2bids(cfg)
end
fprintf('\nConversion of %s done!\n\n',cfg.task)

% EOF
