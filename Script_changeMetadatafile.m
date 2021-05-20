% running ReSaveMetadata

savefilename = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210506_testingMetadataMods/VNC_JRC_SS52985_RigA_20210503T152107/Metadata.xml';
metadata = ReadMetadataFile('/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210506_testingMetadataMods/VNC_JRC_SS52985_RigA_20210503T152107/Metadata.xml.bak.bak.xml');

[success] = ResaveMetadata(metadata,savefilename);