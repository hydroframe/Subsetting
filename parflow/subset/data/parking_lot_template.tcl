

#-----------------------------
# IMPORT PARFLOW TCL PACKAGE
#-----------------------------
set tcl_precision               17
lappend auto_path               $env(PARFLOW_DIR)/bin
package require                 parflow
namespace import                Parflow::*
pfset FileVersion               4


#-----------------
# SET PROCESSORS
#-----------------
pfset Process.Topology.P        2
pfset Process.Topology.Q        1
pfset Process.Topology.R        1


#------------------
# SPECIFY RUNNAME
#------------------
set runname spinup_test


#-----------------------
# CREATE OUTPUT DIRECTORY
#-----------------------
file mkdir ./pf_runs
cd "./pf_runs"
file mkdir $runname
cd $runname


#-----------------------
# COPY PARFLOW INPUTS
#-----------------------
file copy -force ../../input_files/sfs_slopex.pfb .
file copy -force ../../input_files/sfs_slopey.pfb .
file copy -force ../../input_files/sfs.pfsol .
file copy -force ../../input_files/sfs_PME.pfb .


#-----------------------
# COPY CLM INPUTS
#-----------------------
# file copy -force "../clm_input/sfs_clmin.dat" .
# file copy -force "../clm_input/sfs_vegp.dat"  .
# file copy -force "../clm_input/sfs_vegm.dat"  .
puts "Files Copied"


#---------------------
# COMPUTATIONAL GRID
#---------------------
pfset ComputationalGrid.Lower.X   0.0
pfset ComputationalGrid.Lower.Y   0.0
pfset ComputationalGrid.Lower.Z   0.0
pfset ComputationalGrid.NX        64
pfset ComputationalGrid.NY        128
pfset ComputationalGrid.NZ        1
pfset ComputationalGrid.DX        1000.0
pfset ComputationalGrid.DY        1000.0
pfset ComputationalGrid.DZ        1000.0

#------------------------------------------------------
# TIMING (time units is set by units of permeability)
#  - units: hr
#------------------------------------------------------
pfset TimingInfo.BaseUnit        1.0
pfset TimingInfo.StartCount      0.0
pfset TimingInfo.StartTime       0.0
pfset TimingInfo.StopTime        2000.0
pfset TimingInfo.DumpInterval    50.0
pfset TimeStep.Type              Constant
pfset TimeStep.Value             1.0


#---------------------------------------------------------------
# TIME CYCLES
#  - units: depends on base unit (set in TIMING)
#---------------------------------------------------------------
#pfset Cycle.Names                       "constant rainrec"
pfset Cycle.Names                       "constant"
pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length     10000000
pfset Cycle.constant.Repeat             -1
#pfset Cycle.rainrec.Names               "rain rec"
#pfset Cycle.rainrec.rain.Length         2.
#pfset Cycle.rainrec.rec.Length          300.
#pfset Cycle.rainrec.Repeat              -1


#-------------------------------------
# INITIAL CONDITIONS: WATER PRESSURE
#-------------------------------------
pfset ICPressure.Type             			HydroStaticPatch
pfset ICPressure.GeomNames                  domain
pfset Geom.domain.ICPressure.Value          0
pfset Geom.domain.ICPressure.RefGeom        domain
pfset Geom.domain.ICPressure.RefPatch       bottom


#-------------------------------------
# SPINUP KEYS
#-------------------------------------
#pfset OverlandFlowSpinUp  		1
#pfset OverlandSpinupDampP1 		1.0
#pfset OverlandSpinupDampP2 		0.00001


#----------------------
# SPECIFY TOPO SLOPES
#----------------------
pfset TopoSlopesX.Type        "PFBFile"
pfset TopoSlopesX.GeomNames   "domain"
pfset TopoSlopesX.FileName    sfs_slopex.pfb
pfset TopoSlopesY.Type        "PFBFile"
pfset TopoSlopesY.GeomNames   "domain"
pfset TopoSlopesY.FileName    sfs_slopey.pfb

#--------------------------------
# BOUNDARY CONDITIONS: PRESSURE
#--------------------------------
pfset BCPressure.PatchNames                     "land top bottom"
pfset Solver.EvapTransFile                      False
#pfset Solver.EvapTrans.FileName                 sfs_PME.pfb

# LAND
pfset Patch.land.BCPressure.Type		        FluxConst
pfset Patch.land.BCPressure.Cycle		        "constant"
pfset Patch.land.BCPressure.alltime.Value	    0.0

# TOP
pfset Patch.top.BCPressure.Type                 OverlandFlow
pfset Patch.top.BCPressure.Cycle                "rainrec"
pfset Patch.top.BCPressure.rain.Value           -0.05
pfset Patch.top.BCPressure.rec.Value            0.0000
pfset Patch.top.BCPressure.Cycle                "constant"
pfset Patch.top.BCPressure.alltime.Value        -0.5

# BOTTOM
pfset Patch.bottom.BCPressure.Type		        FluxConst
pfset Patch.bottom.BCPressure.Cycle		        "constant"
pfset Patch.bottom.BCPressure.alltime.Value	    0.0


#---------------------------
# NAMES OF GEOMETRY INPUTS
#---------------------------
pfset GeomInput.Names                  "domaininput"
pfset GeomInput.domaininput.GeomName   domain
pfset GeomInput.domaininput.InputType  SolidFile
pfset GeomInput.domaininput.GeomNames  domain
pfset GeomInput.domaininput.FileName   sfs.pfsol
pfset Geom.domain.Patches              "land top bottom"


#--------------------
# SUBSURFACE LAYERS
#--------------------
pfset Solver.Nonlinear.VariableDz   True
pfset dzScale.GeomNames             domain
pfset dzScale.Type                  nzList
pfset dzScale.nzListNumber          1
pfset Cell.0.dzScale.Value          0.5
pfset Cell.1.dzScale.Value          0.5
pfset Cell.2.dzScale.Value          0.5
pfset Cell.3.dzScale.Value          0.5
pfset Cell.4.dzScale.Value          0.5


#-----------------------------
# PERMEABILITY (units: m/hr)
#-----------------------------
pfset Geom.Perm.Names                  "domain"
pfset Geom.domain.Perm.Type            Constant
pfset Geom.domain.Perm.Value           0.02849
pfset Perm.TensorType                  TensorByGeom
pfset Geom.Perm.TensorByGeom.Names     "domain"
pfset Geom.domain.Perm.TensorValX      1.0d0
pfset Geom.domain.Perm.TensorValY      1.0d0
pfset Geom.domain.Perm.TensorValZ      1.0d0


#------------------------
# RELATIVE PERMEABILITY
#------------------------
pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames          "domain"
pfset Geom.domain.RelPerm.Alpha        1.
pfset Geom.domain.RelPerm.N            2.


#-------------------
# SPECIFIC STORAGE
#-------------------
pfset SpecificStorage.Type                Constant
pfset SpecificStorage.GeomNames           "domain"
pfset Geom.domain.SpecificStorage.Value   1.0e-4


#-----------
# POROSITY
#-----------
pfset Geom.Porosity.GeomNames           "domain"
pfset Geom.domain.Porosity.Type         Constant
pfset Geom.domain.Porosity.Value        0.39738


#-------------
# SATURATION
#-------------
pfset Phase.Saturation.Type               VanGenuchten
pfset Phase.Saturation.GeomNames          "domain"
pfset Geom.domain.Saturation.Alpha        1.0
pfset Geom.domain.Saturation.N            2.
pfset Geom.domain.Saturation.SRes         0.2
pfset Geom.domain.Saturation.SSat         1.0


#-----------------------
# MANNINGS COEFFICIENT
#-----------------------
pfset Mannings.Type                  "Constant"
pfset Mannings.GeomNames             "domain"
pfset Mannings.Geom.domain.Value     0.0000044


#---------
# PHASES
#---------
pfset Phase.Names                       "water"
pfset Phase.water.Density.Type	        Constant
pfset Phase.water.Density.Value	        1.0
pfset Phase.water.Viscosity.Type	    Constant
pfset Phase.water.Viscosity.Value	    1.0


#----------------
# PHASE SOURCES
#----------------
pfset PhaseSources.water.Type                     Constant
pfset PhaseSources.water.GeomNames                domain
pfset PhaseSources.water.Geom.domain.Value        0.0


#----------
# GRAVITY
#----------
pfset Gravity				1.0


#---------
# DOMAIN
#---------
pfset Domain.GeomName     domain


#-----------
# MOBILITY
#-----------
pfset Phase.water.Mobility.Type        Constant
pfset Phase.water.Mobility.Value       1.0


#---------------
# CONTAMINANTS
#---------------
pfset Contaminants.Names			""


#--------------
# RETARDATION
#--------------
pfset Geom.Retardation.GeomNames           ""


#--------
# WELLS
#--------
pfset Wells.Names                         ""


#------------------------------------------------------
# EXACT SOLUTION SPECIFICATION FOR ERROR CALCULATIONS
#------------------------------------------------------
pfset KnownSolution                                    NoKnownSolution


#------------------------------------------
# DISTRIBUTE INPUTS
#------------------------------------------
pfdist sfs_slopex.pfb
pfdist sfs_slopey.pfb
pfdist sfs_PME.pfb


#------------------------
# SET SOLVER PARAMETERS
#------------------------
pfset Solver                                             Richards
pfset Solver.MaxIter                                     2500
pfset Solver.Nonlinear.MaxIter                           80
pfset Solver.TerrainFollowingGrid                        True
pfset Solver.Nonlinear.ResidualTol                       1e-5
pfset Solver.Nonlinear.EtaChoice                         EtaConstant
pfset Solver.Nonlinear.EtaValue                          1e-3
pfset Solver.Nonlinear.UseJacobian                       True
pfset Solver.Nonlinear.DerivativeEpsilon                 1e-16
pfset Solver.Nonlinear.StepTol                           1e-25
pfset Solver.Linear.KrylovDimension                      500
pfset Solver.Linear.MaxRestarts                          8
pfset Solver.MaxConvergenceFailures                      5
pfset Solver.Linear.Preconditioner                       PFMG
pfset Solver.Linear.Preconditioner.PCMatrixType          FullJacobian


# PmE flux
pfset Solver.EvapTransFile    False
pfset Solver.EvapTrans.FileName  UC_PME.1layer.flux.pfb

#----------
# OUTPUTS
#----------
pfset Solver.WriteSiloPressure                        False
pfset Solver.WriteSiloSpecificStorage                 False
pfset Solver.WriteSiloMannings                        False
pfset Solver.WriteSiloMask                            False
pfset Solver.WriteSiloSlopes                          False
pfset Solver.WriteSiloSubsurfData                     False
pfset Solver.WriteSiloPressure                        False
pfset Solver.WriteSiloSaturation                      False
pfset Solver.WriteSiloEvapTrans                       False
pfset Solver.WriteSiloEvapTransSum                    False
pfset Solver.WriteSiloOverlandSum                     False
pfset Solver.WriteSiloCLM                             False
pfset Solver.WriteCLMBinary                           False
pfset Solver.PrintSubsurfData                         True
pfset Solver.PrintMask                                True
pfset Solver.PrintVelocities                          False
pfset Solver.PrintSaturation                          True
pfset Solver.PrintPressure                            True
pfset Solver.PrintSubsurfData                         True


#-----------------
# RUN SIMULATION
#-----------------
puts    $runname
pfrun   $runname


#-----------------------
# UNDISTRIBUTE OUTPUTS
#-----------------------
pfundist $runname
pfundist sfs_slopex.pfb
pfundist sfs_slopey.pfb
#pfundist sfs_PME.pfb
puts "ParFlow run complete"