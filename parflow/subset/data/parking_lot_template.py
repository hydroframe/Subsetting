from parflow.tools import Run


def get_parking_lot_model_boundary(patch_name):
    bcs = dict(land={'Type': 'FluxConst', 'Cycle': 'constant', 'alltime.Value': 0.0},
               top={'Type': 'OverlandKinematic', 'Cycle': 'rainrec', 'rain.Value': -0.05, 'rec.Value': 0.0},
               bottom={'Type': 'FluxConst', 'Cycle': 'constant', 'alltime.Value': 0.0},
               sink={'Type': 'DirEquilRefPatch', 'Cycle': 'constant', 'RefGeom': 'domain', 'RefPatch': 'sink',
                      'alltime.Value': 0.0},
               lakes={'Type': 'DirEquilRefPatch', 'Cycle': 'constant', 'RefGeom': 'domain', 'RefPatch': 'lake',
                      'alltime.Value': 0.0},
               ocean={'Type': 'FluxConst', 'Cycle': 'constant', 'alltime.Value': 0.0})
    return bcs.get(patch_name)


def get_parking_lot_model(out_name, slopex_file, slopey_file, solid_file, NX, NY, NZ=1):
    parking_lot_template = Run(out_name)

    parking_lot_template 
    parking_lot_template.FileVersion = 4

    parking_lot_template.Process.Topology.P = 2
    parking_lot_template.Process.Topology.Q = 1
    parking_lot_template.Process.Topology.R = 1

    # ---------------------------------------------------------
    # Computational Grid
    # ---------------------------------------------------------
    parking_lot_template.ComputationalGrid.Lower.X = 0.0
    parking_lot_template.ComputationalGrid.Lower.Y = 0.0
    parking_lot_template.ComputationalGrid.Lower.Z = 0.0

    parking_lot_template.ComputationalGrid.NX = NX
    parking_lot_template.ComputationalGrid.NY = NY
    parking_lot_template.ComputationalGrid.NZ = NZ

    parking_lot_template.ComputationalGrid.DX = 1000
    parking_lot_template.ComputationalGrid.DY = 1000
    parking_lot_template.ComputationalGrid.DZ = 1000.0

    # ---------------------------------------------------------
    # The Names of the GeomInputs
    # ---------------------------------------------------------
    parking_lot_template.GeomInput.Names = 'domaininput'

    parking_lot_template.GeomInput.domaininput.GeomName = 'domain'

    parking_lot_template.GeomInput.domaininput.InputType = 'SolidFile'
    parking_lot_template.GeomInput.domaininput.GeomNames = 'domain'
    parking_lot_template.GeomInput.domaininput.FileName = solid_file

    parking_lot_template.Geom.domain.Patches = 'land top bottom'

    # --------------------------------------------
    # variable dz assignments
    # ------------------------------------------
    parking_lot_template.Solver.Nonlinear.VariableDz = True
    parking_lot_template.dzScale.GeomNames = 'domain'
    parking_lot_template.dzScale.Type = 'nzList'
    parking_lot_template.dzScale.nzListNumber = 1

    parking_lot_template.Cell['0'].dzScale.Value = 0.0001

    # -----------------------------------------------------------------------------
    # Perm
    # -----------------------------------------------------------------------------

    parking_lot_template.Geom.Perm.Names = 'domain'

    parking_lot_template.Geom.domain.Perm.Type = 'Constant'
    parking_lot_template.Geom.domain.Perm.Value = 0.000001

    parking_lot_template.Perm.TensorType = 'TensorByGeom'

    parking_lot_template.Geom.Perm.TensorByGeom.Names = 'domain'

    parking_lot_template.Geom.domain.Perm.TensorValX = 1.0
    parking_lot_template.Geom.domain.Perm.TensorValY = 1.0
    parking_lot_template.Geom.domain.Perm.TensorValZ = 1.0

    # -----------------------------------------------------------------------------
    # Specific Storage
    # -----------------------------------------------------------------------------

    parking_lot_template.SpecificStorage.Type = 'Constant'
    parking_lot_template.SpecificStorage.GeomNames = 'domain'
    parking_lot_template.Geom.domain.SpecificStorage.Value = 1.0e-4

    # -----------------------------------------------------------------------------
    # Phases
    # -----------------------------------------------------------------------------

    parking_lot_template.Phase.Names = 'water'

    parking_lot_template.Phase.water.Density.Type = 'Constant'
    parking_lot_template.Phase.water.Density.Value = 1.0

    parking_lot_template.Phase.water.Viscosity.Type = 'Constant'
    parking_lot_template.Phase.water.Viscosity.Value = 1.0

    # -----------------------------------------------------------------------------
    # Contaminants
    # -----------------------------------------------------------------------------

    parking_lot_template.Contaminants.Names = ''

    # -----------------------------------------------------------------------------
    # Retardation
    # -----------------------------------------------------------------------------

    parking_lot_template.Geom.Retardation.GeomNames = ''

    # -----------------------------------------------------------------------------
    # Gravity
    # -----------------------------------------------------------------------------

    parking_lot_template.Gravity = 1.0

    # -----------------------------------------------------------------------------
    # Setup timing info
    # -----------------------------------------------------------------------------

    #
    parking_lot_template.TimingInfo.BaseUnit = 0.01
    parking_lot_template.TimingInfo.StartCount = 0
    parking_lot_template.TimingInfo.StartTime = 0.0
    parking_lot_template.TimingInfo.StopTime = 200.0
    parking_lot_template.TimingInfo.DumpInterval = 2
    parking_lot_template.TimeStep.Type = 'Constant'
    parking_lot_template.TimeStep.Value = 0.01

    # -----------------------------------------------------------------------------
    # Porosity
    # -----------------------------------------------------------------------------

    parking_lot_template.Geom.Porosity.GeomNames = 'domain'

    parking_lot_template.Geom.domain.Porosity.Type = 'Constant'
    parking_lot_template.Geom.domain.Porosity.Value = 0.00000001

    # -----------------------------------------------------------------------------
    # Domain
    # -----------------------------------------------------------------------------

    parking_lot_template.Domain.GeomName = 'domain'

    # -----------------------------------------------------------------------------
    # Relative Permeability
    # -----------------------------------------------------------------------------

    parking_lot_template.Phase.RelPerm.Type = 'VanGenuchten'
    parking_lot_template.Phase.RelPerm.GeomNames = 'domain'

    parking_lot_template.Geom.domain.RelPerm.Alpha = 1.0
    parking_lot_template.Geom.domain.RelPerm.N = 2.

    # ---------------------------------------------------------
    # Saturation
    # ---------------------------------------------------------

    parking_lot_template.Phase.Saturation.Type = 'VanGenuchten'
    parking_lot_template.Phase.Saturation.GeomNames = 'domain'

    parking_lot_template.Geom.domain.Saturation.Alpha = 1.0
    parking_lot_template.Geom.domain.Saturation.N = 2.
    parking_lot_template.Geom.domain.Saturation.SRes = 0.2
    parking_lot_template.Geom.domain.Saturation.SSat = 1.0

    # -----------------------------------------------------------------------------
    # Wells
    # -----------------------------------------------------------------------------
    parking_lot_template.Wells.Names = ''

    # -----------------------------------------------------------------------------
    # Time Cycles
    # -----------------------------------------------------------------------------
    parking_lot_template.Cycle.Names = 'constant rainrec'
    parking_lot_template.Cycle.constant.Names = 'alltime'
    parking_lot_template.Cycle.constant.alltime.Length = 1
    parking_lot_template.Cycle.constant.Repeat = -1

    # rainfall and recession time periods are defined here
    # rain for 1 hour, recession for 2 hours

    parking_lot_template.Cycle.rainrec.Names = 'rain rec'
    parking_lot_template.Cycle.rainrec.rain.Length = 100
    parking_lot_template.Cycle.rainrec.rec.Length = 5000000
    parking_lot_template.Cycle.rainrec.Repeat = -1

    # ---------------------------------------------------------
    # Mannings coefficient
    # ---------------------------------------------------------

    parking_lot_template.Mannings.Type = 'Constant'
    parking_lot_template.Mannings.GeomNames = 'domain'
    parking_lot_template.Mannings.Geom.domain.Value = 0.0000044
    # -----------------------------------------------------------------------------
    # Phase sources:
    # -----------------------------------------------------------------------------

    parking_lot_template.PhaseSources.water.Type = 'Constant'
    parking_lot_template.PhaseSources.water.GeomNames = 'domain'
    parking_lot_template.PhaseSources.water.Geom.domain.Value = 0.0

    # -----------------------------------------------------------------------------
    # Exact solution specification for error calculations
    # -----------------------------------------------------------------------------

    parking_lot_template.KnownSolution = 'NoKnownSolution'

    # -----------------------------------------------------------------------------
    # Set solver parameters
    # -----------------------------------------------------------------------------

    parking_lot_template.Solver = 'Richards'
    parking_lot_template.Solver.MaxIter = 10000000

    parking_lot_template.Solver.Nonlinear.MaxIter = 50
    parking_lot_template.Solver.Nonlinear.ResidualTol = 1e-6

    parking_lot_template.Solver.Nonlinear.EtaChoice = 'EtaConstant'

    parking_lot_template.Solver.Nonlinear.EtaValue = 0.01
    parking_lot_template.Solver.Nonlinear.UseJacobian = False
    parking_lot_template.Solver.Nonlinear.DerivativeEpsilon = 1e-16
    parking_lot_template.Solver.Nonlinear.StepTol = 1e-30
    parking_lot_template.Solver.Nonlinear.Globalization = 'LineSearch'
    parking_lot_template.Solver.Linear.KrylovDimension = 20
    parking_lot_template.Solver.Linear.MaxRestart = 2

    parking_lot_template.Solver.Linear.Preconditioner = 'PFMGOctree'
    parking_lot_template.Solver.Linear.Preconditioner.PCMatrixType = 'PFSymmetric'


    parking_lot_template.Solver.PrintSubsurf = False
    parking_lot_template.Solver.Drop = 1E-30
    parking_lot_template.Solver.AbsTol = 1E-9

    # ---------------------------------------------------------
    # Initial conditions: water pressure
    # ---------------------------------------------------------

    # set water table to be at the bottom of the domain, the top layer is initially dry
    parking_lot_template.ICPressure.Type = 'HydroStaticPatch'
    parking_lot_template.ICPressure.GeomNames = 'domain'
    parking_lot_template.Geom.domain.ICPressure.Value = 0.0

    parking_lot_template.Geom.domain.ICPressure.RefGeom = 'domain'
    parking_lot_template.Geom.domain.ICPressure.RefPatch = 'bottom'

    # -----------------------------------------------------------------------------
    # Boundary Conditions: Pressure
    # -----------------------------------------------------------------------------
    parking_lot_template.BCPressure.PatchNames = 'land top bottom'

    # no flow boundaries for the land borders and the bottom
    parking_lot_template.Patch.land.BCPressure.Type = 'FluxConst'
    parking_lot_template.Patch.land.BCPressure.Cycle = 'constant'
    parking_lot_template.Patch.land.BCPressure.alltime.Value = 0.0

    parking_lot_template.Patch.bottom.BCPressure.Type = 'FluxConst'
    parking_lot_template.Patch.bottom.BCPressure.Cycle = 'constant'
    parking_lot_template.Patch.bottom.BCPressure.alltime.Value = 0.0

    ## overland flow boundary condition with rainfall then nothing
    parking_lot_template.Patch.top.BCPressure.Type = 'OverlandKinematic'
    parking_lot_template.Patch.top.BCPressure.Cycle = 'rainrec'
    parking_lot_template.Patch.top.BCPressure.rain.Value = -0.05
    parking_lot_template.Patch.top.BCPressure.rec.Value = 0.0000

    parking_lot_template.TopoSlopesX.Type = 'PFBFile'
    parking_lot_template.TopoSlopesX.GeomNames = 'domain'
    parking_lot_template.TopoSlopesX.FileName = slopex_file

    parking_lot_template.TopoSlopesY.Type = 'PFBFile'
    parking_lot_template.TopoSlopesY.GeomNames = 'domain'
    parking_lot_template.TopoSlopesY.FileName = slopey_file

    return parking_lot_template