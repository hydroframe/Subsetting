
set tcl_precision 17

#
# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*

set CONUS_X [pfload $innamex.pfb]
set new_slope_X [pfgetsubbox $CONUS_X $il $jl 0 $iu $ju 1]
pfsave $new_slope_X -pfb $outnamex.pfb
pfsave $new_slope_X -silo $outnamex.silo

set CONUS_Y [pfload $innamey.pfb]
set new_slope_Y [pfgetsubbox $CONUS_Y $il $jl 0 $iu $ju 1]
pfsave $new_slope_Y -pfb $outnamey.pfb
pfsave $new_slope_Y -silo $outnamey.silo
