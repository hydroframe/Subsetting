
lappend   auto_path $env(PARFLOW_DIR)/bin
package   require parflow
namespace import Parflow::*

pfset     FileVersion    4

#Converting from pfb to sa
#set   perm  [pfload -silo NA_1km.out.perm_x.silo]
#pfsave $perm	-sa NA_1km.out.perm.sa


#Converting sa to pfb and silo
#set   mask  [pfload -sa pf_mask_LakesSinks.txt]
#pfsetgrid {4442 3256 1} {0.0 0.0 0.0} {1000.0 1000.0 2.0} $mask
#pfsave $mask -silo pf_mask_LakesSinks.silo
#pfsave $mask -pfb  pf_mask_LakesSinks.pfb

#set   slopey  [pfload -sa LW.slopey_mod.sa]
#pfsetgrid {41 41 1} {0.0 0.0 0.0} {1000 1000 2} $slopey
#pfsave $slopey -silo LW.slopey_mod.silo
#pfsave $slopey -pfb  LW.slopey_mod.pfb

#Converting sa to pfb and silo
set   slopex  [pfload -sa $out_name.slopex.sa]
pfsetgrid {$ncol $nrow 1} {0.0 0.0 0.0} {1000.0 1000.0 1000.0} $slopex
pfsave $slopex -silo $out_name.slopex.silo
pfsave $slopex -pfb  $out_name.slopex.pfb

set   slopey  [pfload -sa $out_name.slopey.sa]
pfsetgrid {$ncol $nrow 1} {0.0 0.0 0.0} {1000.0 1000.0 1000.0} $slopey
pfsave $slopey -silo $out_name.slopey.silo
pfsave $slopey -pfb  $out_name.slopey.pfb