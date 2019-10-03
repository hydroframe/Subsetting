#  This runs the tilted-v catchment problem
#  similar to that in Kollet and Maxwell (2006) AWR

set tcl_precision 17

set runname parking_lot_test

#
# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*

pfset FileVersion 4

pfset Process.Topology.P        2
pfset Process.Topology.Q        1
pfset Process.Topology.R        1

#---------------------------------------------------------
# Copy necessary files
#---------------------------------------------------------

#Slope files
file copy -force ../pfbs/Upper_CO.slopex.rivth1500.pfb .
file copy -force ../pfbs/Upper_CO.slopey.rivth1500.pfb .
#Solid file
file copy -force ../pfbs/Upper_CO_v2.pfsol .
#PME file
file copy -force ../pfbs/UC_PME.1layer.flux.pfb .

#---------------------------------------------------------
# Computational Grid
#---------------------------------------------------------
pfset ComputationalGrid.Lower.X           0.0
pfset ComputationalGrid.Lower.Y           0.0
pfset ComputationalGrid.Lower.Z           0.0

pfset ComputationalGrid.NX                608
pfset ComputationalGrid.NY                896
pfset ComputationalGrid.NZ                1

pfset ComputationalGrid.DX	         	1000.0
pfset ComputationalGrid.DY              1000.0
pfset ComputationalGrid.DZ	            1000.0

#---------------------------------------------------------
# The Names of the GeomInputs
#---------------------------------------------------------

pfset GeomInput.Names                 "domaininput"

pfset GeomInput.domaininput.GeomName  domain

pfset GeomInput.domaininput.InputType  SolidFile
pfset GeomInput.domaininput.GeomNames  domain
pfset GeomInput.domaininput.FileName   Upper_CO_v2.pfsol

pfset Geom.domain.Patches             "land top  bottom"

#--------------------------------------------
# variable dz assignments
#------------------------------------------
pfset Solver.Nonlinear.VariableDz   True
pfset dzScale.GeomNames            domain
pfset dzScale.Type            nzList
pfset dzScale.nzListNumber       1

#pfset dzScale.Type            nzList
#pfset dzScale.nzListNumber       3
pfset Cell.0.dzScale.Value 0.5
pfset Cell.1.dzScale.Value 0.5
pfset Cell.2.dzScale.Value 0.5
pfset Cell.3.dzScale.Value 0.5
pfset Cell.4.dzScale.Value 0.5

#-----------------------------------------------------------------------------
# Permeability (values in m/hr)
#-----------------------------------------------------------------------------
pfset Geom.Perm.Names                 "domain"

# Values in m/hour

pfset Geom.domain.Perm.Type             Constant
pfset Geom.domain.Perm.Value            0.02849


pfset Perm.TensorType               TensorByGeom

pfset Geom.Perm.TensorByGeom.Names  "domain"

pfset Geom.domain.Perm.TensorValX  1.0d0
pfset Geom.domain.Perm.TensorValY  1.0d0
pfset Geom.domain.Perm.TensorValZ  1.0d0

#-----------------------------------------------------------------------------
# Specific Storage
#-----------------------------------------------------------------------------

pfset SpecificStorage.Type            Constant
pfset SpecificStorage.GeomNames       "domain"
pfset Geom.domain.SpecificStorage.Value 1.0e-4

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------

pfset Phase.Names "water"

pfset Phase.water.Density.Type	        Constant
pfset Phase.water.Density.Value	        1.0

pfset Phase.water.Viscosity.Type	Constant
pfset Phase.water.Viscosity.Value	1.0

#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------

pfset Contaminants.Names			""

#-----------------------------------------------------------------------------
# Retardation
#-----------------------------------------------------------------------------

pfset Geom.Retardation.GeomNames           ""

#-----------------------------------------------------------------------------
# Gravity
#-----------------------------------------------------------------------------

pfset Gravity				1.0

#-----------------------------------------------------------------------------
# Setup timing info
#-----------------------------------------------------------------------------

#
pfset TimingInfo.BaseUnit        0.1
pfset TimingInfo.StartCount      0.0
pfset TimingInfo.StartTime       0.
pfset TimingInfo.StopTime        300.
#pfset TimingInfo.StopTime	72.0
pfset TimingInfo.DumpInterval    -1
pfset TimeStep.Type              Constant
pfset TimeStep.Value             0.1
#pfset TimeStep.Type Growth
#pfset TimeStep.InitialStep 0.01
#pfset TimeStep.GrowthFactor 1.4
#pfset TimeStep.MaxStep 1.0
#pfset TimeStep.MinStep 0.0001

#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------
pfset Geom.Porosity.GeomNames           "domain"

pfset Geom.domain.Porosity.Type         Constant
pfset Geom.domain.Porosity.Value        0.39738
#pfset Geom.domain.Porosity.Value        0.00000001

#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------

pfset Domain.GeomName domain

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------

pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames      "domain"


pfset Geom.domain.RelPerm.Alpha    1.
pfset Geom.domain.RelPerm.N        2.

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

pfset Phase.Saturation.Type              VanGenuchten
pfset Phase.Saturation.GeomNames         "domain"

pfset Geom.domain.Saturation.Alpha        1.0
pfset Geom.domain.Saturation.N            2.
pfset Geom.domain.Saturation.SRes         0.2
pfset Geom.domain.Saturation.SSat         1.0

#----------------------------------------------------------------------------
# Mobility
#----------------------------------------------------------------------------
pfset Phase.water.Mobility.Type        Constant
pfset Phase.water.Mobility.Value       1.0

#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names                         ""

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names "constant rainrec"
pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length      1
pfset Cycle.constant.Repeat             -1

# rainfall and recession time periods are defined here
# rain for 0.2 hour, recession for 30 hours

pfset Cycle.rainrec.Names                 "rain rec"
pfset Cycle.rainrec.rain.Length           2.
pfset Cycle.rainrec.rec.Length            300.
pfset Cycle.rainrec.Repeat                -1

#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------
pfset BCPressure.PatchNames             "land top  bottom"

# zero head boundaries for ocean, sink and lake boundaries
pfset Patch.ocean.BCPressure.Type       DirEquilRefPatch
pfset Patch.ocean.BCPressure.Type       FluxConst

pfset Patch.ocean.BCPressure.Cycle      "constant"
pfset Patch.ocean.BCPressure.RefGeom        domain
pfset Patch.ocean.BCPressure.RefPatch     bottom
pfset Patch.ocean.BCPressure.alltime.Value  0.

pfset Patch.sink.BCPressure.Type       DirEquilRefPatch
pfset Patch.sink.BCPressure.Type       OverlandFlow

pfset Patch.sink.BCPressure.Cycle      "constant"
pfset Patch.sink.BCPressure.RefGeom        domain
pfset Patch.sink.BCPressure.RefPatch     bottom
pfset Patch.sink.BCPressure.alltime.Value  0.
pfset Patch.sink.BCPressure.Cycle		      "rainrec"
pfset Patch.sink.BCPressure.rain.Value	      -0.05
pfset Patch.sink.BCPressure.rec.Value	      0.0000

pfset Patch.lake.BCPressure.Type       DirEquilRefPatch
pfset Patch.lake.BCPressure.Type       OverlandFlow

pfset Patch.lake.BCPressure.Cycle      "constant"
pfset Patch.lake.BCPressure.RefGeom        domain
pfset Patch.lake.BCPressure.RefPatch     bottom
pfset Patch.lake.BCPressure.alltime.Value  0.
pfset Patch.lake.BCPressure.Cycle		      "rainrec"
pfset Patch.lake.BCPressure.rain.Value	      -0.05
pfset Patch.lake.BCPressure.rec.Value	      0.0000

#no flow boundaries for the land borders and the bottom
pfset Patch.land.BCPressure.Type		      FluxConst
pfset Patch.land.BCPressure.Cycle		      "constant"
pfset Patch.land.BCPressure.alltime.Value	      0.0

pfset Patch.bottom.BCPressure.Type		      FluxConst
pfset Patch.bottom.BCPressure.Cycle		      "constant"
pfset Patch.bottom.BCPressure.alltime.Value	      0.0


## overland flow boundary condition with rainfall then nothing
pfset Patch.top.BCPressure.Type		      OverlandFlow
pfset Patch.top.BCPressure.Cycle		      "rainrec"
pfset Patch.top.BCPressure.rain.Value	      -0.05
pfset Patch.top.BCPressure.rec.Value	      0.0000
pfset Patch.top.BCPressure.Cycle		      "constant"
pfset Patch.top.BCPressure.alltime.Value	      0.0000

# PmE flux
pfset Solver.EvapTransFile    True
pfset Solver.EvapTrans.FileName  UC_PME.1layer.flux.pfb

#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

pfset TopoSlopesX.Type "PFBFile"
pfset TopoSlopesX.GeomNames "domain"

pfset TopoSlopesX.FileName Upper_CO.slopex.rivth1500.pfb


#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------

pfset TopoSlopesY.Type "PFBFile"
pfset TopoSlopesY.GeomNames "domain"

pfset TopoSlopesY.FileName Upper_CO.slopey.rivth1500.pfb

#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------

# set water table to be at the bottom of the domain, the top layer is initially dry
pfset ICPressure.Type                                   HydroStaticPatch
pfset ICPressure.GeomNames                              domain
pfset Geom.domain.ICPressure.Value                      0.0

pfset Geom.domain.ICPressure.RefGeom                    domain
pfset Geom.domain.ICPressure.RefPatch                   bottom

#---------
##  Distribute inputs
#---------

pfset ComputationalGrid.NZ                1 

pfdist Upper_CO.slopex.rivth1500.pfb
pfdist Upper_CO.slopey.rivth1500.pfb
pfdist UC_PME.1layer.flux.pfb

#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

#pfset TopoSlopesX.Type "Constant"
#pfset TopoSlopesX.GeomNames "domain"
#pfset TopoSlopesX.Geom.domain.Value 0.01
#pfset TopoSlopesX.Geom.domain.Value -0.01
#pfset TopoSlopesX.Geom.domain.Value 0.0

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------


#pfset TopoSlopesY.Type "Constant"
#pfset TopoSlopesY.GeomNames "domain"
#pfset TopoSlopesY.Geom.domain.Value 0.001
#pfset TopoSlopesY.Geom.domain.Value -0.01
#pfset TopoSlopesY.Geom.domain.Value 0.0

#---------------------------------------------------------
# Mannings coefficient
#---------------------------------------------------------

pfset Mannings.Type "Constant"
pfset Mannings.GeomNames "domain"
pfset Mannings.Geom.domain.Value 2.e-6
pfset Mannings.Geom.domain.Value 0.0000044
#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------

pfset PhaseSources.water.Type                         Constant
pfset PhaseSources.water.GeomNames                    domain
pfset PhaseSources.water.Geom.domain.Value        0.0

#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------

pfset KnownSolution                                    NoKnownSolution

#-----------------------------------------------------------------------------
# Set solver parameters
#-----------------------------------------------------------------------------
pfset Solver                                             Richards
pfset Solver.MaxIter                                    100000

#pfset Solver.TerrainFollowingGrid                     True
pfset Solver.TerrainFollowingGrid                    True


pfset Solver.Nonlinear.MaxIter                           2000
pfset Solver.Nonlinear.ResidualTol                       1e-5
pfset Solver.Nonlinear.EtaChoice                         EtaConstant
pfset Solver.Nonlinear.EtaValue                          1e-3
pfset Solver.Nonlinear.UseJacobian                       False
pfset Solver.Nonlinear.UseJacobian                       True
pfset Solver.Nonlinear.DerivativeEpsilon                 1e-16
pfset Solver.Nonlinear.StepTol                       1e-25

pfset Solver.Linear.KrylovDimension                      500
pfset Solver.Linear.MaxRestarts                           8
pfset Solver.MaxConvergenceFailures                       5

pfset Solver.Linear.Preconditioner                       PFMG
#pfset Solver.Linear.Preconditioner.MGSemi.MaxIter        1
#pfset Solver.Linear.Preconditioner.MGSemi.MaxLevels      100

pfset Solver.Linear.Preconditioner.PCMatrixType     FullJacobian
pfset Solver.WriteSiloPressure                        False
pfset Solver.PrintSubsurfData True
pfset Solver.PrintMask True
pfset Solver.PrintVelocities   False
pfset Solver.PrintSaturation   False
pfset Solver.PrintPressure   True
#Writing output (no binary except Pressure, all silo):
pfset Solver.PrintSubsurfData                         True 
#pfset Solver.PrintLSMSink                        True 
pfset Solver.PrintSaturation                          True
pfset Solver.WriteCLMBinary                           False

pfset Solver.WriteSiloSpecificStorage                  False
pfset Solver.WriteSiloMannings                         False
pfset Solver.WriteSiloMask                             False
pfset Solver.WriteSiloSlopes                          False
pfset Solver.WriteSiloSubsurfData                     False
pfset Solver.WriteSiloPressure                        False
pfset Solver.WriteSiloSaturation                       False
pfset Solver.WriteSiloEvapTrans                       False
pfset Solver.WriteSiloEvapTransSum                     False
pfset Solver.WriteSiloOverlandSum                     False
pfset Solver.WriteSiloCLM                             False
#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------

pfrun $runname
#pfundist $runname
#pfundist Upper_CO.slopex.rivth1500.pfb
#pfundist Upper_CO.slopey.rivth1500.pfb
#pfundist UC_PME.1layer.flux.pfb

#puts "ParFlow run Complete"

#pfwritedb $runname
