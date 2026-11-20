# YASARA MACRO
# TOPIC:       3. Molecular Dynamics
# TITLE:       Running a steered molecular dynamics simulation in water
# REQUIRES:    Dynamics
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro sets up and runs a steered simulation to pull a ligand from its binding site. As soon as they have been pulled apart completely and atoms start crossing periodic boundaries, the macro works no longer correctly. You can set the pulling 'acceleration' below, which is not changed during the simulation. So this macro can currently not tell you how strong you need to pull, it can only tell you whether the chosen acceleration is strong enough or not to pull the ligand off.

# ADATTAMENTO S2 - versione corretta per una replica YASARA il piu' possibile
# fedele alla SMD descritta da Reilly et al. (2025).
#
# PUNTI CHIAVE:
# 1) La SMD di Reilly e' stata eseguita su Chain A (singolo protomero S2) + ligando.
#    Per la replica piu' fedele, preparare quindi una scena contenente il solo protomero A
#    (oltre all'eventuale ligando) prima di eseguire questa macro.
# 2) La forza dissociativa e' implementata come potenziale armonico instabile riferito alla
#    distanza assoluta r tra i CA terminali, cioe' |F| = |k|*r, con |k|=5.0 kJ/mol/nm^2.
#    NON viene sottratta la distanza iniziale r0: questo evita che la forza sia nulla al tempo
#    zero e che cambi segno in caso di contrazione termica.
# 3) La forza viene convertita in accelerazione, come richiesto da YASARA, e applicata a ogni
#    simulation step ai due CA terminali in direzioni opposte.
# 4) Restano differenze rispetto a Reilly: solvente esplicito + PME e termostato YASARA
#    invece di OBC2 + Langevin Middle. Il protocollo va quindi descritto come adattato.

# Include library functions
include md_library

# Parameter section - adjust as needed, but NOTE that some changes only take
# effect if you start an entirely new simulation, not if you continue an existing one. 
# ====================================================================================

# The structure to simulate must be present with a .pdb or .sce extension.
# If a .sce (=YASARA scene) file is present, the cell must have been added.
# You can either set the target structure by clicking on Options > Macro > Set target,
# by providing it as command line argument (see docs at Essentials > The command line),
# or by uncommenting the line below and specifying it directly.
#MacroTarget '/home/myself/test'

# Selection of ligand residue(s), adapt as needed. For PDB file 1BFB, you would need
# ligandsel='Res UAB SGN IDS'. If the ligand is a protein, you can simply use its molecule name,
# e.g. 'Mol B'. If you leave ligandsel empty, YASARA will choose selected atoms as the ligand or,
# if none are selected, try to identify the ligand automatically.
ligandsel='CA Res 1 Mol A'

# Selection of receptor residues. By default, the ligand will be pulled away from the center
# of mass of the selected residues. If the receptor has an unusual shape, select residues
# immediately 'below' the ligand, e.g. the active site 'Res Ser 360 Glu 361'. If the ligand is
# also a protein, add the receptor's molecule name, e.g. 'Mol A'. 
# If you provide an arrow to indicate the pulling direction, then select all receptor atoms,
# YASARA will pick the right subset to update the pulling direction automatically. To verify that
# your selection is right, just look at the red pulling arrow shown by YASARA during the steered MD.
receptorsel='CA Res 339 Mol A'

# Simulation temperature
# Allineata a Reilly et al. (integratore Langevin a 300K nel paper). Se cambi la temperatura,
# adatta anche densità e pressure control più sotto.
temperature='300K'

# Solvent density in g/ml (0.997 is water at 298K; ~0.9965 a 300K)
# If you run at a significantly different temperature, this density needs to
# be adapted. Look at the PressureCtrl documentation.
# For solvents other than water, you have to create your own solvent box
# as described at the FillCellObj documentation and save it as .._solvent.sce
density=0.9965

# Name of solvent molecules used to measure solvent density
solvent='HOH'

# Imposta l'incremento a 0 per avere una forza costante (non usato più: vedi bondforce_const)
accdelta=0

# Valore assoluto della custom bond force di Reilly et al. [kJ/mol/nm^2].
# La costante pubblicata e' -5.0; il segno repulsivo e' implementato applicando
# accelerazioni opposte ai due terminali, con modulo proporzionale alla distanza assoluta r.
bondforce_const=5.0

# Limite opzionale sull'accelerazione [pm/ps^2].
# 0 = nessun cap (raccomandato per la replica del protocollo pubblicato).
# Impostare un valore >0 solo come protezione numerica, dichiarandolo nel Methods.
accmax=0

# Durata della fase di steered MD in nanosecondi (100 ns in Reilly et al.)
pulltime_ns=100

# Time in picoseconds when the pulling starts
# (everything before can be considered equilibration)
pullstart=3

# pH at which the simulation should be run, by default physiological pH 
ph=4.5

# The ion concentration as a mass fraction, here we use 0.9% NaCl (physiological solution)
# NOTA: Reilly et al. hanno usato solvente implicito (OBC2), quindi senza ioni espliciti; qui
# restiamo nell'approssimazione a solvente esplicito discussa in precedenza.
ions='Na,Cl,0.9'

# Multiple timestep: 1.25 femtoseconds for intramolecular and 2*1.25 fs for intermolecular forces
tslist=2,1.25

# Save simulation snapshots every 250000 femtoseconds = 250 picoseconds
saveinterval=250000

# The format used to save the trajectories: YASARA 'sim', GROMACS 'xtc' or AMBER 'mdcrd'.
format='sim'

# Rete di sicurezza secondaria [A]: il criterio di arresto principale ora e' il tempo fisso
# (pulltime_ns, come nel paper), non piu' la distanza. Questo valore alto serve solo a fermare
# la simulazione in caso di comportamento anomalo (es. attraversamento dei confini periodici).
enddis=100

# Use longrange coulomb forces (particle-mesh Ewald)
longrange='Coulomb'

# Forcefield to use (these are all YASARA commands, so no '=' used)
ForceField AMBER14

# Cutoff
Cutoff 8

# Number of simulation steps per screen update and pairlist update.
# Since the macro can apply pulling forces only once per screen update, screen and pairlist are updated each
# simulation step to ensure that forces are applied continuously. This slows down the simulation considerably
# to about 40% of the normal speed. If speed is an issue, you can change the command below to 'SimSteps 2,2'
 # (50% of normal speed), 'SimSteps 3,3' (60% of normal speed) or 'SimSteps 5,5' (80% of normal speed). Going
# higher is not recommended, since acceleration*N will be applied every Nth step, and no acceleration in between.
SimSteps Screen=1,Pairlist=1

# Normally no change required below this point
# ============================================

# Requires 'with distance from local point'
RequireVersion 19.1.27

# Keep the solute from diffusing around and crossing periodic boundaries
CorrectDrift On

# Treat all simulation warnings as errors that stop the macro
WarnIsError On

# Normally no change required below this point
# Temperature
Temp (temperature)

# Calculate total timestep, we want a float, so tslist2 is on the left side
ts=tslist2*tslist1
# Snapshots are saved every 'savesteps'
savesteps=saveinterval/ts
# Tempo target di pulling in picosecondi (coerente con l'unita' di 't' nel ciclo piu' sotto)
pulltime_ps=pulltime_ns*1000

Console Off

# Do we have a target?
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"

Clear
# Do we already have a scene with water or other solvent?
waterscene = FileSize (MacroTarget)_water.sce
solventscene = FileSize (MacroTarget)_solvent.sce
if waterscene
  LoadSce (MacroTarget)_water
elif solventscene 
  LoadSce (MacroTarget)_solvent
else
  # No scene with water present yet
  if solvent!='HOH'
    RaiseError "When using a solvent other than water, the cell cannot be filled automatically. Look at the documentation of the FillCellObj command"
  # Do we have a scene at all?
  scene = FileSize (MacroTarget).sce
  if scene
    LoadSce (MacroTarget)
  else
    # No scene present, assume it's a PDB or YOB file
    for type in 'yob','pdb'
      size = FileSize (MacroTarget).(type)
      if size
        break
    if !size
      RaiseError 'Initial structure not found, expected (MacroTarget).pdb or .yob. Make sure to create a project directory and place the structure there'
    # Load structure
    Load(type) (MacroTarget)
    # Align model with major axes to minimize cell size
    NiceOriAll
    CleanAll
    # Create the simulation cell, make sure that there is enough space to reach 'enddis'
    Cell Auto,Extension=(enddis*2)
    SaveSce (MacroTarget)
  # Fill with water, predict pKas, place counter ions
  Boundary periodic
  Experiment Neutralization
    WaterDensity (density)
    pH (ph)
    Ions (ions)
    pKaFile (MacroTarget).pka
    Speed Fast
  Experiment On
  Wait ExpEnd
  # Save scene with water
  SaveSce (MacroTarget)_water

# Identify the ligand
if ligandsel==''
  # Check if atoms are selected
  selatomlist()=ListAtom selected
  if count selatomlist
    ligandsel='(join selatomlist)'
  else
    ligandsel=FindLigand 0

for type in 'ligand','receptor'
  # Make sure that ligand and receptor are really present
  if (type)sel==''
    (type)atoms=0
  else
    (type)atoms = CountAtom ((type)sel)
  if !(type)atoms
    RaiseError 'This macro requires that the (type) is selected explicitly. The current (type) ((type)sel) was not found'
  # Get central ligand/receptor atom for arrow visualization
  cen()=GroupCenter ((type)sel),Type=Mass
  for atomsel in 'CA','!H'
    (type)atom=ListAtom (atomsel) and ((type)sel) with minimum distance from local point (join cen)
    if (type)atom
      break

# Get ligand mass
ligmasslist()=MassAtom (ligandsel)
ligmass=sum ligmasslist
recmasslist()=MassAtom (receptorsel)
recmass=sum recmasslist

# For S2 terminal pulling the direction is always defined by the instantaneous
# vector joining the C- and N-terminal CA atoms of Chain A.  Ignore any arrow object.
arrowobj=0
if arrowobj
  cellsize1,cellsize2,cellsize3,angle1,angle2,angle3=Cell
  if angle1!=90 or angle2!=90 or angle3!=90
    RaiseError "Steered MD with a user-provided arrow requires a cuboid cell (with all three cell angles at 90 degrees)"
  diag=norm cellsize
  # Create a chain of dummy atoms
  chainobj=BuildGrid 0,(2*diag)
  # Copy global position/orientation of arrow to chain of dummy atoms
  TransferObj (chainobj),(arrowobj),Local=Keep
  # Determine number of atoms that will be used to place the arrow
  neighboratoms = CountAtom (receptorsel)
  if neighboratoms>20
    neighboratoms=20
  # Select all receptor atoms near the chain
  cutoff=0.5
  do
    arrowatomlist()=ListAtom (receptorsel) with distance<(cutoff) from Obj (chainobj)
    cutoff=cutoff+0.5
  while count arrowatomlist<neighboratoms
  # Determine direction vector of arrow in the coordinate system of the simulation cell from chainobj.
  # This requires wall boundaries to prevent the long chainobj from being wrapped around.
  boundary=Boundary
  Boundary Wall
  Longrange None
  PressureCtrl Off
  Sim On
  _,_,_,lastdir()=GroupLine Obj (chainobj)
  DelObj (arrowobj) (chainobj)
  Boundary (boundary)
  Longrange (longrange)

# IMPORTANT: Reilly et al. extracted Chain A before SMD.  This macro selects Mol A terminal CAs.
# For the closest reproduction, the input scene should therefore contain only protein Mol A
# (plus the bound ligand, if present).  If other protein protomers are retained, the force law
# remains correct but the simulation becomes a trimer-context adaptation rather than an exact replica.

# Start the Simulation
TimeStep (tslist)
Sim On
# Make sure all atoms are free to move
FreeAll

# Do a short steepest descent energy minimization
TempCtrl SteepDes
ShowMessage "Starting simulation with a steepest descent minimization..."
# Wait for convergence
for i=1 to 500
  Wait 1
  if SpeedMax<4000
    break
# Continue with 500 simulated annealing steps
ShowMessage "Continuing with 500 simulated annealing steps..."
Temp 0
TempCtrl Anneal
Wait 500

# And now start with the real simulation
HideMessage
Temp (temperature)
# Reset time to 0
Time 0
# Set temperature and pressure control
TempCtrl Rescale
# NOTA: Reilly et al. hanno usato un integratore Langevin; qui si resta con il controllo di
# temperatura nativo di questa macro (Rescale). Se vuoi provare un'opzione piu' vicina a un
# termostato stocastico stile Langevin, controlla la documentazione YASARA del comando TempCtrl.
PressureCtrl SolventProbe,(solvent),(density)

# Save snapshots
trajectfilename='(MacroTarget)00000.sim'
if format!='sim'
  trajectfilename='(MacroTarget).(format)'
Save(format) (trajectfilename),(savesteps)

# Run the steered simulation
MakeTab Default,2,4
disstart=-1
dismax=-1
work=0
progtime=0
step=0
while 1
  t = Time
  t=t/1000
  if arrowobj
    # Keep the arrow direction updated
    _,_,_,dir()=GroupLine (join arrowatomlist)
    if sign (sum (dir*lastdir))<0
      dir()=-dir
    lastdir()=dir
  if t>pullstart
    # Steer the simulation
    ligpos()=GroupCenter (ligandsel),Type=Mass
    if arrowobj
      # Calculate pulling 'dir'ection from arrowatomlist, and 'dis'tance of the ligand from its starting position
      recpos()=GroupCenter (join arrowatomlist),Type=Mass
    else
      # Calculate pulling direction from centers of mass, and 'dis'tance between receptor and ligand
      recpos()=GroupCenter (receptorsel),Type=Mass
    # Absolute N/C-terminal distance r [A].  The published custom bond coefficient
    # is applied to r itself, not to an extension (r-r0).
    disvec()=ligpos-recpos
    rabs=norm disvec
    if disstart==-1
      disstart=rabs
    if !arrowobj
      dir()=disvec/rabs
    # Extension relative to the start is retained only for reporting and the safety stop.
    dis=rabs-disstart
    HideArrowAtom (ligandatom)
    ShowArrow AtAtom,(ligandatom),StepFromAtom,(ligandatom),(dir*10)

    # Criterio di arresto principale: tempo fisso di simulazione, come in Reilly et al. (100 ns)
    if t>=pulltime_ps
      ShowMessage 'Raggiunto il tempo target di (pulltime_ns) ns in (0.00+t) ps, simulazione terminata'
      break
    # Rete di sicurezza secondaria (non presente nel paper): ferma comunque se la distanza esplode
    if dis>enddis
      ShowMessage 'La distanza ha superato (enddis) A in (0.00+t) ps per sicurezza, simulazione terminata'
      break

    # ---- Reilly-like dissociative custom bond force ----
    # With |k| in kJ/mol/nm^2 and rabs in Angstrom:
    #   |F|[pN] = 0.1660539 * |k| * rabs[A]
    # YASARA applies acceleration.  For a selected atom/group of mass m[g/mol]:
    #   a[pm/ps^2] = 100 * |k| * rabs[A] / m[g/mol]
    # Separate accelerations are calculated for the two terminal selections so that
    # the same force magnitude acts on both sides even if their selected masses differ.
    accelerationN=bondforce_const*rabs*100/ligmass
    accelerationC=bondforce_const*rabs*100/recmass
    if accmax>0
      if accelerationN>accmax
        accelerationN=accmax
      if accelerationC>accmax
        accelerationC=accmax

    # Force magnitude reported from the N-terminal selection (same target force on C-term).
    force=1e21/AvoConst*accelerationN*ligmass
    acceleration=accelerationN
    if dis>dismax
      progtime=t
      if (dismax>0)
        # work[EnergyUnit] = dis[A]*1e-10[m/A]*force[pN]*1e-12[N/pN]*JToUnit[EnergyUnit/J]
        work=JToUnit*1e-22*((dis-dismax)*force)+work
      Tabulate (dis),(t),(acceleration),(force)
      dismax=dis
    # Get current steps per screen update
    steps=SimSteps
    ShowMessage 'Steered MD: t=(0.00+t)/(pulltime_ps) ps, r=(0.00+rabs) A, estensione=(0.00+dis) A, forza=(0.00+force) pN'
    # Apply the acceleration on both termini in opposite directions. The 'CorrectDrift' command
    # further above ensures that the system's center of mass does not move. NOTE that this recipe
    # only works as long the complex doesn't cross a periodic boundary.

    # 'ligandsel' = N-terminal CA, 'receptorsel' = C-terminal CA.
    # dir points from C-term to N-term in this inherited macro convention.
    # Apply equal-magnitude opposing forces by using the mass-adjusted accelerations.
    AccelAtom (ligandsel),(dir*accelerationN*steps)
    AccelAtom (receptorsel),(dir*-accelerationC*steps)
  else
    ShowMessage 'Equilibrazione, la steered MD comincia tra (0.00+pullstart-t) picosecondi...'  
  # Proceed with one simulation step
  step=step+1
  Wait 1
Sim off

# Save table and plots
SaveTab Default,(MacroTarget),Format=Text,Columns=4,NumFormat=12.3f,"_Extension[A]____Time[ps]_AccN[pm/ps^2]___Force[pN]"
SavePlot Filename=(MacroTarget)_acc,Default,Width=1600,Height=1200,Title='Steered MD: total work = (0.00+work) (EnergyUnit)',
         XColumn=1,YColumn=3,XLabel='Extension from start [A]',YLabel='Acceleration [pm/ps^2]'
SavePlot Filename=(MacroTarget)_frc,Default,Width=1600,Height=1200,Title='Steered MD: total work = (0.00+work) (EnergyUnit)',
         XColumn=1,YColumn=4,XLabel='Extension from start [A]',YLabel='Force [pN]'

# Exit YASARA if this macro was provided as command line argument in console mode
if runWithMacro and ConsoleMode
  Exit
