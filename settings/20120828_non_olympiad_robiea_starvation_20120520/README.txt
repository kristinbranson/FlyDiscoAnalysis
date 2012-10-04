protocol 20120828_non_olympiad_robiea_starvation_20120520:

copied from 20120530_non_olympiad_robiea_starvation_20120520 then merged with 20120706. 

Differences from 20120706:

automatic checks are less stringent.
behavior detection parameters are not divided by sex.
extra perframefns: absdphi, absdv_ctr, absthetadiff_center, absphisideways

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
