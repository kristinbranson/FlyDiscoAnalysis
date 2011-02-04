function RegisterTrxWrapper(expdir,protocol,varargin)

[settingsdir,registration_paramsfilestr,dataloc_paramsfilestr] = ...
    myparse(varargin,...
            'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
            'registration_paramsfilestr','registration_params.txt',...
            'dataloc_paramsfilestr','dataloc_params.txt');

% read the parameters
registration_paramsfile = fullfile(settingsdir,protocol,registration_paramsfilestr);
registration_params = ReadParams(registration_paramsfile);
dataloc_paramsfile = fullfile(settingsdir,protocol,dataloc_paramsfilestr);
dataloc_params = ReadParams(dataloc_paramsfile);

% full path to experiment directory
expdir_read = fullfile(dataloc_params.rootreaddir,expdir);
expdir_write = fullfile(dataloc_params.rootwritedir,expdir);

% full path to registration data file
registrationfile = fullfile(expdir_write,dataloc_params.registrationfilestr);

% full path to output trx file
outtrxfile = fullfile(expdir_write,dataloc_params.trxfilestr);

% full path to input trx file
ctraxfile = fullfile(expdir_write,dataloc_params.ctraxfilestr);

% full path to annotation file
annfile = fullfile(expdir_write,dataloc_params.annfilestr);

% full path to movie file
moviefile = fullfile(expdir_read,dataloc_params.moviefilestr);

RegisterTrx(ctraxfile,'annname',annfile,'moviename',moviefile,...
  'outregistrationfile',registrationfile,...
  'detectregistrationparams',struct2paramscell(registration_params),...
  'outtrxname',outtrxfile);

if ~exist(outtrxfile,'file'),
  error('Registered trx file %s not created',outtrxfile);
end

if ~exist(registrationfile,'file'),
  error('Registration data file %s not created',registrationfile);
end

fprintf('Successfully registered trx in experiment %s\n',expdir);