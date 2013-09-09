changes on 20130909:

Updated JAABA projects.

Added a few new per-frame stats:

du_ctr,attemptedcopulation,any,5,500
absdv_cor,attemptedcopulation,any,5,500

velmag_ctr,backup,any,5,500

velmag_ctr,bodyturn,any,5,500

velmag_ctr,crabwalkall,any,5,500

velmag_ctr,crabwalkextreme,any,5,500

du_ctr,jump,any,5,500
absdv_cor,jump,any,5,500

velmag_ctr,pivotcenter,any,5,500

velmag_ctr,pivottail,any,5,500

absdv_cor,righting,any,5,500
du_ctr,righting,any,5,500

absdv_cor,stop,any,5,500
du_ctr,stop,any,5,500

absdv_cor,move,any,5,500
du_ctr,move,any,5,500
velmag_ctr,move,any,5,500

absdv_cor,touch,any,5,500
du_ctr,touch,any,5,500

du_ctr,winggrooming,any,5,500
absdv_cor,winggrooming,any,5,500


Updated FlyBowlPlotPerFrameStats2 to not fail when there are no frames for the plotted condition. 

Updated FlyBowlJAABADetect to work with the latest version of JAABA (commit 523e9108f8f27cb42f49ee6118ea166570f7cb43)

changes on 20130722:

automatic_checks_incoming:

made some files desired instead of required (FlyBowlDataCaptureParams_EP*.txt,Log.txt,QuickStats.png,QuickStats.txt,ufmf_diagnostics.txt,ufmf_log.txt)
added more debug info to the output file.
recompiled

register_trx:

recompiled.
added more debug info to the output file.

classify_sex:

recompiled
added more debug info to the output file. 

wing tracking:

created FlyBowlTrackWings wrapper function

compute per-frame features:

Added debug info to each per-frame mat file. 
Added compute functions for wing features. 
Added wing features to per-frame features to compute -- now includes everything except nflies_close_*

jaabadetect:

Added debug info file. 
Updated classifier list from 
/groups/branson/home/robiea/Projects_data/JAABA/ProjectFiles_jab/Pivot_tail2_v8p4.jab
to
/groups/branson/home/robiea/Projects_data/JAABA/ProjectFiles_jab/Pivot_tail2_v9p1.jab

compute per-frame stats:


changes on 20121105:

removed
nearwall,max_dist2wall,15
nearfly,max_dcenter,6.5
from fly conditions. these weren't being used in stats_perframedfeatures.txt, so hopefully it doesn't change anything. 



protocol 20120706:

changes on 20120819:

Added recompiled FlyBowlPlotPerFrameStats2. 

changes on 20120801:

Fixed issue in ClassifySex if both classes are unlikely when learning the HMM parameters. 

changes on 20120730:

Fixed bug in computing magveldiff when the a fly was closest on the last frame of its trajectory. Recompiled ComputePerFrameFeatures and FixPerFrameFeatures. 

changes on 20120725:

Changed to ComputePerFrameStats2 instead of ComputePerFrameStats: new averaging scheme.
Compiled ComputePerFrameStats2. 

changes on 20120706:

added usemediandt to registration parameters. this means that we will just use the median framerate rather than the per-frame timestamps. 
modified registration code to accept usemediandt parameter. 
fixed bug in compute_dv_tail, compute_du_tail: was using theta instead of theta_mm. 
fixed bug in compute_magveldiff: was not normalizing by dt
changed per-frame features computed (added the new social and locomotion features)

recompiled stuff with 2012a:
AutomaticChecks_Incoming
RegisterTrx
ClassifySex2
ComputePerFrameFeatures


*********************************************************************

protocol 20120621:

20120621: small bug fix in FlyBowlComputePerFrameFeatures to not delete per-frame feature mat files that are trx files. 
20120705: Recompiled FlyBowlComputePerFrameStats to use new version of Trx. 
