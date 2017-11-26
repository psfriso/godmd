!
! CommLine & I/O
!
   integer, parameter :: NFILES = 7
 
   integer unit_i, unit_o,unit_ener, unit_traj,unit_targ
 
   type(commLineOption) :: files(NFILES) = (/&
   commLineOption("-i",    "param",          "formatted",   "old",     "Settings"),&
   commLineOption("-pdbin", "raw PDB file",  "formatted",   "old",     "Initial PDB"),&
   commLineOption("-ener", "energy",         "formatted",   "unknown", "Energies"),&
   commLineOption("-trj", "trajectory.pdb", "formatted",   "unknown", "Trajectory (PDB)"),&
   commLineOption("-pdbtarg",  "target.pdb", "formatted",   "old",     "Target PDB"),&
   commLineOption("-o",    "log",            "formatted",   "unknown", "Calculation Log"),&
   commLineOption("-lig",   "ligand.pdb",    "formatted",   "old",     "Ligand PDB")/)

  INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)! Double precision real number definition processor independent
  INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6) ! Single Precision real number definition processor independent

 ! Coordinates and velocities 
   type(pointDP), allocatable :: r(:), v(:)
   type(point), allocatable :: rsp(:),rspTarg(:) ! Single precision version of r for binary I/O

   type(point) :: rcm ! C of M  
   type(pointDP), allocatable :: rtarg(:)
         
 
! Step Potentials
   type (stepPotInt), allocatable :: stepPts(:,:)
   type (intDataTop), allocatable :: stepPtsDef(:,:)

! Structure
   character(len=4), allocatable :: res(:), atp(:)
   character(len=4), allocatable :: atoms(:)
   character(len=1), allocatable :: chain(:)
   integer, allocatable :: rnum(:)
   integer natom, nres
   integer nmol
   integer, allocatable :: molNum(:), molRes(:), mols(:,:)
   
   ! Pedro
   REAL*8 :: error
   real*8, dimension(3,3) :: U               ! Rotation Matrix
   real*8, dimension(3) :: center1, center2  ! Center of initial and target structures
   real, allocatable :: w(:)
   type(pointDP), allocatable :: rprev(:)
   integer :: iacc, auxiliar
   real *8 :: score, scoreprev
   integer :: accepted, totalTrials
   real :: aux_acc
   logical :: stopping
   real, save :: DRMSd=1.E-4
   integer :: contarTodos
   !
   
 ! Potentials & energies
   real :: evdw,rvdw,xm, rhc, xsum
   real xmassa, calcEkin, ekin, epotfis, epotgo, ekin0
   real epot(MAXTIPINT)
   
! Interaction pairs
   type(intpList), allocatable :: blist(:), nblist(:)

! Collisions
   integer :: ierr
   integer ibloc, iev
   logical, allocatable :: toUpdate(:)
! Time
   real*8 tacact, taccorr, tacum, temps, tacrect
   real*8 tinit, tsetup, tfin
   real*8, allocatable :: tpart(:)
   integer, allocatable :: ipart(:)
!
   integer i,j, ioerr
!
! Old setup
   type(struc) ::  str, strTarg,strLig,newStr
   type(stepPotIntDefList) :: stpList
   type(residueLibrary) :: resLib
   type(ffprm) :: ff
   integer :: numEnB=0
   real::fractionDone=0.00
   real:: initiation
   integer ::axran
!
! PEDRO
   real :: distHC = 0.0
   real*8, allocatable :: dist_INI(:,:)
   real*8, allocatable :: dist_TARG(:,:)
   character(len=50) :: fmtsteppot='(2I4,I2,x,I1,4(2(f7.1,X),L1,X),L1,X,f9.4,f8.3)'
!SELF STOP
  REAL(SGL), ALLOCATABLE, DIMENSION(:,:) :: error_evo      ! rmsd vs time vector
  REAL(SGL) :: slope=0._SGL                                ! rmsd vs time slope
  INTEGER :: errcont=0                                     ! index to write new rmsd value
!SELF STOP2
  REAL(SGL), ALLOCATABLE, DIMENSION(:,:) :: error_long     ! rmsd vs time vector
  REAL(SGL) :: slopeLong=0._SGL                                ! rmsd vs time slope
  INTEGER :: errcontLong=0                                     ! index to write new rmsd value

! DOUBLE WELL DMD
    real, allocatable, dimension (:,:) :: discardedEner
    integer :: saco
    integer :: sacoLimit
! ENCHUFA MAS META-DMD
  REAL(SGL), ALLOCATABLE, DIMENSION(:,:) :: ener_evo  ! pot energy vs time vector
  REAL(SGL) :: slopeEner=0._SGL                       ! slope of pot energy vs time
  INTEGER :: EnerCont=0                               ! index to write new energy value
  LOGICAL :: enchufaMeta
! Normal Modes
  TYPE (pointDP), ALLOCATABLE, DIMENSION(:,:) :: evec
  REAL, ALLOCATABLE, DIMENSION (:) :: evals
  integer :: contarBondedInt=0
  integer :: contarsaltos
   real*8 :: rgyr=0.0
   real*8 :: error0
   real :: errorMin
   LOGICAL :: calcAuxAcc
   LOGICAL :: unForceIt
   INTEGER :: natomBare
   INTEGER :: contarIntLig=0
   INTEGER :: nwritten=0
  