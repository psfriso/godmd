 MODULE paramSet
     use stepPotentials
! Input param
!
   real, save :: &

      temp = 300., &
      tmin = 1.e-22*1.e15
  integer, save :: &
      isolv = 1, &
      nbloc = 1000, &
      idab = 1, &
      igoab = 1, &
      iwr = 1, &
      seed = 2381, &
      TCALC = 1, & ! E Const 0,  T Const 1, TBany 2
      outformat = 0, &
      idims = 0, &
      rst=0, &
      rstv=0
  real*8, save :: &
      tsnap = 1000., &
      tcorr = 100., &
      tact = 10., &
      trect = 100., &
      tini = 0.
      
  real,save :: acceptance=0.70 ! DEFAULT 0.65 NO MODIFICAR MUCHO. (0.5-0.8) 
  real :: xbeta=2.0000 !2
  real :: goener=0.5000
  real :: sigmago=0.1000
  real *8 :: dtol=0.30000
  real :: wellWidth=0.1300            ! DEFAULT 0.13
  real :: scalFactor=40.00   !40
  real :: rgyrPercent=0.500           ! DEFAULT=0.52
  real *8 :: rcut2GO4=0.0             !  depends on radius gyration
  real *8 :: rcut1= 70.00
  real *8 :: mxRcut4go=16.0
  real *8 :: ssectol=dble(0.035)
  real :: minEnerWell=0.05         ! DEFAULT 0.05
 LOGICAL :: StartMD=.FALSE.        ! IGUAL QUE ENER_evo_size
 INTEGER :: error_evo_size=20        ! number of rmsd point to use in linear regression
 INTEGER :: error_long_size=250
 INTEGER :: Ener_evo_size=20         ! DEFAULT 15
 INTEGER , PARAMETER :: NEVECS=5
 real*8 :: hstep=1.0
 REAL , PARAMETER :: DEGO_BASE=0.045 ! 5 porciento del pozo a sacar en cada metaDMD
 REAL, SAVE :: rmsdStop=0.5
LOGICAL :: LIGAND=.FALSE.
REAL :: errorAcceptable=1.00
INTEGER :: nskip=5

LOGICAL :: writingStandard = .TRUE.

CONTAINS
!===============================================================================
 subroutine readInputParamSet (unit_i)
   integer, intent(IN) :: unit_i
!   
   namelist /input/ tsnap,temp,seed,nbloc,rmsdStop,scalFactor,Ener_evo_size,&
   LIGAND,Ener_evo_size,nskip,errorAcceptable,writingStandard
!
   read (unit_i, INPUT)
   ! checking 
   if (TRECT.gt.TSNAP) TRECT = TSNAP
   if (TCORR.gt.TRECT) TCORR = TRECT
   if (TACT.gt.TCORR)  TACT = TCORR

 end subroutine readInputParamSet
!===============================================================================
 subroutine writeInputParamSet (unit_o)
   integer, intent(IN) :: unit_o
   ! Copiem l'arxiu de parametres. Pendent format
   write (unit_o, *)
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | CALCULATION PARAMETERS                                   |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Simulation settings                                      |")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Simulation Time (ps) (Nbloc x TSnap)       |",f12.3," |")') NBLOC * TSNAP / 1.e3
   write (unit_o, '(" | Output structure (fs)             | TSnap  |",f12.3," |")') TSNAP 
   if (IDIMS.eq.1) &
   write (unit_o, '(" | Re-scoring target (fs)            | Trect  |",f12.3," |")') TRECT
   write (unit_o, '(" | Update velocities (fs)            | Tcorr  |",f12.3," |")') TCORR
   write (unit_o, '(" | Update Lists, collision times (fs)| Tact   |",f12.3," |")') TACT
   write (unit_o, '(" | Min. accepted colision time (fs)  | TMin   |",f12.8," |")') TMIN   
   write (unit_o, '(" | Temperature (K)                   | Temp   |",6X,f6.2, " |")') TEMP
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Other                                                    |")')  
   write (unit_o, '(" ------------------------------------------------------------")')
   write (unit_o, '(" | Random generator seed                      |",7X,i5  " |")') seed
   write (unit_o, '(" | IDAB, IGOAB, IWR, ISOLV                    |",4X,4i2," |")') IDAB, IGOAB, IWR, ISOLV
   write (unit_o, '(" ------------------------------------------------------------")')
 end subroutine writeInputParamSet
!===============================================================================
 END MODULE paramSet
