copied from FlyDiscoAnalysis/settings/20210806_flybubble_LED

modified the ctrax results movie params for the 20220414 LED protocol

20220609
changed VNC_barcode to VNC2__barcode in incoming checks

20220913
added jab classifiers

202201005
added parameters for APT tracking

20230213 from - 20220913_flybubble_LED_VNC2
adding parameters for apt results movie

20230307
updated the lbl file, same name in apt_params.txt 
Same tracker just updated the slots for tracking and training to work with APT update

20230518
added Alice updated walk jaaba classifier -> scores_walk2.mat

20230622
updated the singularity image pointer for cluster upgrade

20230907
added missing parameter to datalocs 
flytrackerbgstr,movie-bg.mat

20240204
creating new settings directory for multiBubble assay - data collected using Jin's data aquisition
analysis-protocol turned off all stages, to be turn back on as set up param files. 
datalocs - changed metedata capitalization, parent calibration file
autochecksincoming - include required_fns and experiment specific fields, added flytracker parent calibration field
ran tests with auto_incoming and tracking, both ran correctly

20240207
turned on registration in analysis-protocol
created and added templates for registration mark and indicator LEDS (these will need to be updated for new plate) 

