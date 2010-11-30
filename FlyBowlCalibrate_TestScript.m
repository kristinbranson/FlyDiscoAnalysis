% FlyBowlCalibrate_TestScript

expdir = 'E:\Data\FlyBowl\Test1\GMR_12E07_AE_01_TrpA_Rig1Plate01BowlA_20101013T095913';
protocol = '0.1';
protocolFile = 'FlyBowlProtocol.txt';

trx = FlyBowlCalibrate(expdir,protocol,'protocolFile',protocolFile);