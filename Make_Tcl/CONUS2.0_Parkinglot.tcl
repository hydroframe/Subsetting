#  This runs the tilted-v catchment problem
#  similar to that in Kollet and Maxwell (2006) AWR

set tcl_precision 17

set runname CONUS2_Parkinglot

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
# Computational Grid
#---------------------------------------------------------
pfset ComputationalGrid.Lower.X           0.0
pfset ComputationalGrid.Lower.Y           0.0
pfset ComputationalGrid.Lower.Z           0.0

pfset ComputationalGrid.NX                4442
pfset ComputationalGrid.NY                3256
#pfset ComputationalGrid.NX                400
#pfset ComputationalGrid.NY                300
pfset ComputationalGrid.NZ                1

pfset ComputationalGrid.DX	         1000
pfset ComputationalGrid.DY               1000
pfset ComputationalGrid.DZ	            1000.0

#---------------------------------------------------------
# The Names of the GeomInputs
#---------------------------------------------------------
pfset GeomInput.Names                 "domaininput"

pfset GeomInput.domaininput.GeomName  domain

pfset GeomInput.domaininput.InputType  SolidFile
pfset GeomInput.domaininput.GeomNames  domain
pfset GeomInput.domaininput.FileName   CONUS2.0_test1.pfsol

pfset Geom.domain.Patches             "ocean land top lake sink bottom"

#--------------------------------------------
# variable dz assignments
#------------------------------------------
pfset Solver.Nonlinear.VariableDz   True
pfset dzScale.GeomNames            domain
pfset dzScale.Type            nzList
pfset dzScale.nzListNumber       1

#pfset dzScale.Type            nzList
#pfset dzScale.nzListNumber       3
pfset Cell.0.dzScale.Value 0.0001

#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------

pfset Geom.Perm.Names                 "domain"


pfset Geom.domain.Perm.Type            Constant
pfset Geom.domain.Perm.Value           0.01

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
pfset TimingInfo.BaseUnit        0.10
pfset TimingInfo.StartCount      0
pfset TimingInfo.StartTime       0.0
pfset TimingInfo.StopTime        1.0
pfset TimingInfo.DumpInterval    -1
pfset TimeStep.Type              Constant
pfset TimeStep.Value             0.1

#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------

pfset Geom.Porosity.GeomNames          "domain"

pfset Geom.domain.Porosity.Type          Constant
pfset Geom.domain.Porosity.Value         0.00000001

#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------

pfset Domain.GeomName domain

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------

pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames          "domain"

pfset Geom.domain.RelPerm.Alpha         1.0
pfset Geom.domain.RelPerm.N             2.

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

pfset Phase.Saturation.Type              VanGenuchten
pfset Phase.Saturation.GeomNames         "domain"

pfset Geom.domain.Saturation.Alpha        1.0
pfset Geom.domain.Saturation.N            2.
pfset Geom.domain.Saturation.SRes         0.2
pfset Geom.domain.Saturation.SSat         1.0



#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names                           ""

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names "constant rainrec"
pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length      1
pfset Cycle.constant.Repeat             -1

# rainfall and recession time periods are defined here
# rain for 1 hour, recession for 2 hours

pfset Cycle.rainrec.Names                 "rain rec"
pfset Cycle.rainrec.rain.Length           2
pfset Cycle.rainrec.rec.Length            300
pfset Cycle.rainrec.Repeat                -1

#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------
pfset BCPressure.PatchNames                   "ocean land top lake sink bottom"

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

#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

pfset TopoSlopesX.Type PFBFile
pfset TopoSlopesX.GeomNames "domain"

pfset TopoSlopesX.FileName CONUS.slopex.pfb


#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------

pfset TopoSlopesY.Type PFBFile
pfset TopoSlopesY.GeomNames "domain"

pfset TopoSlopesY.FileName CONUS.slopey.pfb

#---------
##  Distribute inputs
#---------

pfset ComputationalGrid.NZ                1 

pfdist CONUS.slopex.pfb
pfdist CONUS.slopey.pfb

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
pfset Solver.MaxIter                                     2500

pfset Solver.Nonlinear.MaxIter                           50
pfset Solver.Nonlinear.ResidualTol                       1e-6
#pfset Solver.Nonlinear.ResidualTol                       1e-7

pfset Solver.Nonlinear.EtaChoice                         EtaConstant

pfset Solver.Nonlinear.EtaValue                          0.01
pfset Solver.Nonlinear.UseJacobian                       True
pfset Solver.Nonlinear.DerivativeEpsilon                 1e-16
pfset Solver.Nonlinear.StepTol				 1e-30
pfset Solver.Nonlinear.Globalization                     LineSearch
pfset Solver.Linear.KrylovDimension                      20
pfset Solver.Linear.MaxRestart                           2

pfset Solver.Linear.Preconditioner                       PFMGOctree
#pfset Solver.Linear.Preconditioner                       PFMG

pfset Solver.PrintSubsurf				False
pfset  Solver.Drop                                      1E-30
pfset Solver.AbsTol                                     1E-9

pfset Solver.WriteSiloSubsurfData True
pfset Solver.WriteSiloPressure True
pfset Solver.WriteSiloSlopes True

pfset Solver.WriteSiloSaturation True
pfset Solver.WriteSiloConcentration True

#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------

# set water table to be at the bottom of the domain, the top layer is initially dry
pfset ICPressure.Type                                   HydroStaticPatch
pfset ICPressure.GeomNames                              domain
pfset Geom.domain.ICPressure.Value                      0.0

pfset Geom.domain.ICPressure.RefGeom                    domain
pfset Geom.domain.ICPressure.RefPatch                   bottom

#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#-----------------------------------------------------------------------------

pfrun $runname
pfundist $runname
