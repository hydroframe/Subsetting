
set tcl_precision 17

#
# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*

set CONUS_X [pfload Str5Ep0_unsmth.mx0.5.mn5.sec0.5.stan_slopex.pfb]
set new_slope_X [pfgetsubbox $CONUS_X 282 1389 0 607 1673 1]
pfsave $new_slope_X -pfb sub_Str5Ep0_unsmth.mx0.5.mn5.sec0.5.stan_slopex.pfb
pfsave $new_slope_X -silo sub_Str5Ep0_unsmth.mx0.5.mn5.sec0.5.stan_slopex.silo

set CONUS_Y [pfload Str5Ep0_unsmth.mx0.5.mn5.sec0.5.stan_slopey.pfb]
set new_slope_Y [pfgetsubbox $CONUS_Y 282 1389 0 607 1673 1]
pfsave $new_slope_Y -pfb sub_Str5Ep0_unsmth.mx0.5.mn5.sec0.5.stan_slopey.pfb
pfsave $new_slope_Y -silo sub_Str5Ep0_unsmth.mx0.5.mn5.sec0.5.stan_slopey.silo

