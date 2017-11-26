!===========================================================================
!  DISCRETE 0.2.6a
!
!  Discrete Molecular Dynamics,
!
!
!==========================================================================
program discrete
    use commLine
    use geometry
    use geometryDP
    use stepPotentials
    use intList
    use paramSet
    use dims_utils
    use ls_rmsd
    use least_sq_fit
    use Structure
    use ran_mod
    use doubleDMD
    use ANM
    !
#include "vars.h"
 
  call readCommLine(files, unit_i, unit_o, NFILES)
  call readInputParamSet(unit_i)
 !
 ! SetUp
 ! 
  str = loadStructurePDB(openFn(files, '-pdbin'))
  strTarg= loadStructurePDB(openFn(files, '-pdbtarg'))
  IF(LIGAND) then
    !  write(unit_o,*) "ENTRO"
      strLig= loadStructureLIG(openFn(files, '-lig'))
      newStr = mergeStructures(str, strLig )
  END IF

 !
 call setBonds(str)
 !
 call  setStepPotentials(str)
 !
 !
 nmol=str%nmol
 nres=str%nres
 natomBare=str%natoms
 !
 IF(LIGAND) THEN
     natom=newStr%natoms
 ELSE
     natom=natomBare
 END IF
 !
 
 write (unit_o, '(" Molecules: "i4," Residues:   ",i4," Atoms:  ",i6)') nmol, nres, natomBare
 IF(LIGAND)  write (unit_o, '(" Ligand: "i4," Atoms:  ",i6)') strLig%nmol, strLig%natoms
 IF(LIGAND)  write (unit_o, '(" Total:  "i4," Atoms:  ",i6)') strLig%nmol+nmol, natom

 allocate (mols(nmol, 2), stat = ioerr)
 call errorAllocmem(ioerr, 'NMol')
 allocate (molres(nres), stat = ioerr)
 call errorAllocmem(ioerr, 'NRes')
 allocate (atoms(natom), rnum(natom), res(natom), atp(natom), chain(natom),molnum(natom), &
           stepPtsDef(natomBare,natomBare),stat = ioerr)
 call errorAllocmem(ioerr, 'NAtom') 
 !
 !  Lectura de la estructura final
 !
 allocate (rtarg(natomBare),rspTarg(natomBare), stat = ioerr)
 call errorAllocmem(ioerr, 'Coordinates')
 
 !
 IF (LIGAND) THEN
 atoms=newStr%ats(:) %atomId
 rnum= newStr %ats(:) %resnum
 res= newStr %ats(:) %resIdByAtom
 chain= newStr %ats(:) %chainId
 atp= newStr %ats(:) %atType
 molnum= newStr %ats(:) %molnum
 stepPtsDef=str %idata 
 ELSE
 atoms=str%ats(:) %atomId
 rnum= str %ats(:) %resnum
 res= str %ats(:) %resIdByAtom
 chain= str %ats(:) %chainId
 atp= str %ats(:) %atType
 molnum= str %ats(:) %molnum
 stepPtsDef=str %idata 
 
 END IF
 
 !
 ! Coordinates
 allocate (rsp(natom), r(natom), stat = ioerr)
 call errorAllocmem(ioerr, 'Coordinates')
 IF (LIGAND) THEN
 rsp=newStr % ats(:) % r
 ELSE
 rsp=str % ats(:) % r
 END IF
 
 rspTarg= strTarg % ats(:) %r
!
    do i = 1, natom
        r(i) = SPtoDP(rsp(i))
    enddo
    !
    do i = 1, natomBare
        rtarg(i) = SPtoDP(rspTarg(i))
    enddo

    ! 
 deallocate(rsp,rspTarg)   
 !
 !END SETUP
 !
    call cpu_time(tinit)
    call readCommLine(files, unit_i, unit_o, NFILES)
    call programHeader(unit_o)
    call printFileSummary(files, unit_o)
    !
    call writeInputParamSet(unit_o)
    !
    ! Todos los atomos son CA es una asignacion a una lista entera
    evdw = 0.15
    rvdw = 2.0 
    rhc  = 1.6
    xm   = 0.012
    xsum = 1.0/xm !0.5 * (1./xm + 1./xm)
    distHC=(rhc + rhc)**2
    xmassa = xm*natom !sum(xm(1:natom))
    !
    !
    unit_ener = openFn(files, '-ener')
    !
    ! Necesito tambien las distancias ( en A)
    allocate(dist_INI(natom,natom),dist_TARG(natomBare,natomBare))
    dist_INI=0.0
    dist_TARG=0.0
    ! CALL ERROR!!!
    DO j=1,natom
        DO i=1,natom
            dist_INI(i,j)=calcDistDP(r(i),r(j))
            !
        END DO
    END DO
    !
      DO j=1,natomBare
        DO i=1,natomBare
            !
            dist_TARG(i,j)=calcDistDP(rtarg(i),rtarg(j))
        END DO
    END DO
    !
    
   call decideRCut4GO(natomBare,r,rtarg,xm,rgyrPercent,mxRcut4go,rgyr,rcut2GO4)
   Write(unit_o,*) "PNS rgyr  ", rgyr, " CutOff ", sqrt(rcut2GO4)
    !
    ! SELF STOP
    !
    allocate( error_evo(error_evo_size,2 ) ,error_long(error_long_size,2), stat=ioerr)
    error_long=0.0
    ! CALL ERROR!!
    error_evo=0.0 
    !    !
    allocate( ener_evo(ener_evo_size,2 ) , stat=ioerr)
    ! CALL ERROR!!
    ener_evo=0.0
    !
    !
    allocate (v(natom), stepPts(natom, natom),stat = ioerr)
    call errorAllocmem(ioerr, 'Setup/Natom')
    !
    !
   ! Estructura de la proteina
     do j = 2, natomBare
        do i = 1, j-1
                ! Las particulas interaccionan
            stepPts(i, j) = getStepPotFromDEF(GOENER,&
                sigmago,&
                dtol, &
                wellWidth,&
                dist_INI(i,j),&
                dist_TARG(i, j),&
                stepPtsDef(i, j) % tipInt)
                stepPts(i, j) % dmin = distHC
                stepPts(i, j) % dego = 0.0
                stepPts(i,j) % saltos =10
             !     
              stepPts(j,i)=stepPts(i,j)
        enddo
    enddo
    deallocate (stepPtsDef)
 ! Estructura del ligando
    IF (LIGAND) THEN
         DO j=natomBare+1,natom
             DO i=natomBare,j-1
                 stepPts(j, i) % tipInt = COVB
                 stepPts(i, j) % tipInt = COVB
                 stepPts(j, i) % step(1) %r = dist_INI(i,j) *DBLE(0.95)
                 stepPts(i, j) % step(1) %r = dist_INI(i,j) *DBLE(0.95)
                 stepPts(j, i) % step(2) %r = dist_INI(i,j) *DBLE(1.05)
                 stepPts(i, j) % step(2) %r = dist_INI(i,j) *DBLE(1.05)
                 stepPts(i,j)  % nstep = 0
                 stepPts(j,i)  % nstep = 0
                 stepPts(i,j)  % saltos = 0
                 stepPts(j,i)  % saltos = 0
             END DO
         END DO
    END IF
    !
    ! Interaccion Proteina Ligando
    IF (LIGAND) THEN
         DO j=1,natomBare
             DO i=natomBare+1,natom
          !       write(unit_o,*) calcDistDP(r(i),r(j))
                 IF (calcDistDP(r(i),r(j)) < 8.0 ) THEN
                 stepPts(j, i) % tipInt = 10
                 stepPts(i, j) % tipInt = 10
                 stepPts(j, i) % step(1) =stepPot(dist_INI(i,j) -0.7*wellWidth, -goener * FACTE,.TRUE.)
                 stepPts(j, i) % step(2) =stepPot(dist_INI(i,j) +0.7*wellWidth,  goener * FACTE,.TRUE.)
                 stepPts(i,j)  % nstep = 2
                 stepPts(j,i)  % nstep = 2
                 stepPts(i,j)  % saltos = 0
                 stepPts(j,i)  % saltos = 0
                 contarIntLig=contarIntLig+1
              END IF
             END DO
         END DO
    END IF
    !
    ! a por los cambios en ca dist i i+1
    DO j=1,natomBare-1
          i=j+1
             IF ( stepPts(j, i) % tipInt == COVB) THEN
                IF( (ABS(dist_INI(i,j)-dist_TARG(i,j)))/MIN(dist_INI(i,j),dist_TARG(i,j)) > 0.02_DBL) THEN
                  WRITE(unit_o,*) "UNO QUE SE MUEVE" , i, j,dist_INI(i,j),dist_TARG(i,j)
                  stepPts(j, i) % step(1) %r = REAL(MIN(dist_INI(i,j) ,dist_TARG(i,j)))*0.98
                  stepPts(i, j) % step(1) %r = REAL(MIN(dist_INI(i,j),dist_TARG(i,j)))*0.98
                  stepPts(j, i) % step(2) %r = REAL(MAX(dist_INI(i,j),dist_TARG(i,j)))*1.02
                  stepPts(i, j) % step(2) %r = REAL(MAX(dist_INI(i,j),dist_TARG(i,j)))*1.02
                END IF
             END IF
     END DO
    
    IF(LIGAND) write(unit_o,*) "Interaction Ligand", contarIntLig,  " Energy", contarIntLig*goener, natom
    !    
    call fixDistances(natom, natomBare,dist_INI,dist_TARG,ssectol,stepPts)
    !
    allocate(discardedEner(natomBare,natomBare),stat = ioerr)
    ! CALL ERROR
    discardedEner=0.0    
    !
    ! NM
    !
    ALLOCATE(evec(nevecs,natomBare),evals(nevecs),stat = ioerr)
    DO i=1,nevecs
       evec(i,:)=pointDP(0._DBL,0._DBL,0._DBL)
    ENDDO
    
    call ANM1(natomBare,r,nevecs,evec,evals)
    !
    DO i=1,nevecs
        WRITE (unit_o, '( "EIGENVALUE ",i3, " Valor ", f10.5  )')i,evals(i)
    END DO
    !
    call modulateDEgo(unit_o,natom,natomBare,evec,evals,r,dist_INI,stepPts)
    !
    DEALLOCATE(evec,evals)
    
 
    ! precalcul parells
    allocate(blist(natom), nblist(natom), stat = ioerr)
    call errorAllocMem(ioerr, 'Pair lists')
    
    do i = 1, natom
        blist(i) = allocateintPList(MAX(natom-natomBare,3), ioerr)
        call errorAllocMem(ioerr, 'Bonded List')
    enddo
    
    do j = 2, natomBare
        do i = 1, j - 1
            if (stepPts(i, j) % tipInt .eq. COVB ) then
                stepPts(i, j) % active = .false.
                stepPts(j, i) % active = .false.
! j > i only 
                IF(dist_INI(i,j) >20.0 ) THEN
                    stepPts(i, j) % tipInt=0
                    stepPts(j, i) % tipInt=0
                ELSE
                blist(i) % nats = blist(i) % nats + 1
                blist(i) % iData(blist(i) % nats) = intData(j, blist(j) % nats, stepPts(i, j), xsum, 1.e15, 0.)
              !  write(unit_o,*) i,j
                contarBondedInt=contarBondedInt+1
                ENDIF
             endif
          enddo
    enddo
    
    IF (LIGAND) THEN
       do j = natomBare+1, natom
        do i = natomBare, j - 1
            if (stepPts(i, j) % tipInt .eq. COVB ) then
                stepPts(i, j) % active = .false.
                stepPts(j, i) % active = .false.
! j > i only 
                IF(dist_INI(i,j) >20.0 ) THEN
                    stepPts(i, j) % tipInt=0
                    stepPts(j, i) % tipInt=0
                ELSE
                blist(i) % nats = blist(i) % nats + 1
                blist(i) % iData(blist(i) % nats) = intData(j, blist(j) % nats, stepPts(i, j), xsum, 1.e15, 0.)
              !  write(unit_o,*) i,j
!                contarBondedInt=contarBondedInt+1
                ENDIF
             endif
          enddo
        enddo
    END IF
    
    
    do i = 1, natom
        nblist(i) = allocateintPList(INT(0.6*natom), ioerr)
        call errorAllocMem(ioerr, 'Non Bonded List')
    enddo
        
    call activateStepPot(stepPts, r,LIGAND,dist_INI,dist_TARG ,rcut1,rcut2GO4, natom,natomBare, nblist, xsum)
    !
    DEALLOCATE(dist_TARG,dist_INI)
   
    call writePairList(unit_o,natom,contarBondedInt,stepPts)
    
    deallocate (atp)

    !
    call cpu_time(tsetup)
    ! suma l'energia potencial de la conformacio inicial
    call calcEpot(natom,natomBare, r, stepPts, epotgo, epotfis, epot)
    
    rcm = calcCM(natom, r, xm)
    !
    unit_traj = openFn(files, '-trj')

   
    do i = 1, natom
        r(i) = r(i) - SPtoDP(rcm)
    enddo
           
    iev = 0
    call thermalize(seed, iev, natom, temp, xmassa, v, xm, ekin)
    !DIMS
    !
    allocate( w(natomBare) ,rprev(natomBare), stat=ioerr)
    call errorAllocMem(ioerr, 'DIMS')
   
    call rmsd(natomBare, r(1:natomBare), rtarg(:), 1, U, center2, center1, error0, .FALSE.)
    call pdbsuperpos(natom, U , center1, center2,r)
    write(unit_o,*)"RMSD INICIAL ", error0
  
    call MCWeigth(natomBare,r, rtarg,w)
    score = distanceCA(r,rtarg,w,natomBare)
    auxiliar=seed+iev
    scoreprev=1.E10
    rprev=r(1:natomBare)
    iacc=0
    
    IF(.not.writingStandard) THEN
    IF(.not.LIGAND) call writeSnapshot(unit_traj, ibloc, rtarg, atoms, res, chain, rnum, natomBare)
    !
    call writeSnapshot(unit_traj, ibloc, r, atoms, res, chain, rnum, natom)
    END IF
    !
    ekin0 = ekin
    !
    temps = 0.
    allocate (tpart(natom), ipart(natom), toUpdate(natom), stat = ioerr)
    call errorAllocMem(ioerr, ' Colision times')

    Do i=1,Ener_evo_size
        Ener_evo(i,1)=epotgo
        Ener_evo(i,2)=i
    end do
    !
    saco=0
    slope=1.E10
    error_evo=1.E-10
    error_long=1.E-10
    accepted=0
    totalTrials=0
    stopping=.FALSE.
    DRMSd=5.5E-5
    aux_acc=0.
    ibloc=0
    contarTodos=count(stepPts % nstep == 4 .and.stepPts % active)
    EnerCont=0
    enchufaMeta=.FALSE.
    slopeEner=1.
    contarSaltos=0
    aux_acc=1.0
    errorMin=rmsdStop
    calcAuxAcc=.FALSE.
    nwritten=1
    !
  !  Write (unit_o,'("Numero de potenciales con doble pozo ",I5)')contarTodos
  !  open(unit=55,file='saltos.dat') 
  !  contarSaltos=0 
  !    DO j=2,natomBare
  !       DO i= 1, j-1
  ! 
  !       if (stepPts(i,j) % nstep == 4 .and.stepPts(i,j) % active ) write(55,*) stepPts(i,j) % saltos
  !       if (stepPts(i,j) % nstep == 4 .and.stepPts(i,j) % active ) contarSaltos=contarSaltos+1
  !  
  !           
  !       END DO 
  !   END DO
  !   
  !  close(55)
    contarSaltos=0
    unForceIt=.FALSE.
    !
    !
    !----------------------------------------------------------------------------
    do while (.not.stopping)
        
    !
        axran=-ibloc-seed
        initiation=ran1(axran)
        ibloc=ibloc+1
        tacum = 0.
    !    
        EnerCont=EnerCont+1
        IF (EnerCont > Ener_evo_size ) EnerCont=1
        Ener_evo(EnerCont,2) = epotgo
        Ener_evo(EnerCont,1)= real(ibloc)

        slopeEner=lreg(Ener_evo_size,Ener_evo)
        
        IF (slopeEner < 1.0E-1 .and. abs(slopeEner) < 2.5E-2) enchufaMeta=.TRUE.
        
        if( enchufaMeta ) THEN
            startMD=.TRUE.
            contarSaltos=0
            call metaDMD(natom,natomBare, r, goener,stepPts,discardedEner)
            call actualizeStepPot(stepPts,natom,natomBare, nblist, xsum)
            enchufaMeta=.FALSE.
            saco = saco +1
            Ener_evo(:,2)=0.0
              do j = 2, natomBare - 1
                 do i = 1, j-1
                    IF (stepPts(i,j) %active) then
                        IF ( stepPts(i,j) %nstep==4 ) then
                            IF ( stepPts(i,j) % saltos==0 )then
                                contarSaltos=contarSaltos+1
                            end if
                        end if
                    end if
                end do
            end do

        END IF


        !----------------------------------------------------------------------------
        do while (tacum .lt. TSNAP)
            tacrect = 0
            do while (tacrect .lt. TRECT)
                taccorr = 0.
                iev = 0
                ierr = 0
                
                !----------------------------------------------------------------------------
                do while (taccorr .lt. TCORR)
                    ekin0 = calcEkin(v, xm, natom)
                    call colisioBond(blist, r, v, xm, natom)
                    call colisioNonBond(stepPts, temps, nblist, r, v, xm, ierr, &
                    TMIN, natom, ISOLV)
                    call DMDIntLoop(nblist, r, v, xm, tacact, tpart, ipart, TACT, TMIN, temps, iev, natom, toUpdate)
                    taccorr = taccorr + tacact
                    !call thermalize(seed, iev, natom, temp, xmassa, v, xm, ekin)
                enddo
                ! end do while(taccorr.lt.tcorr)-----------------------------------------------
                call thermalize(seed, iev, natom, temp, xmassa, v, xm, ekin)   
                tacrect = tacrect + taccorr
            enddo
        !
         if ( StartMD ) then ! Restraints do not apply up to 3 (arbitrary) time step to allow relaxation in case
         ! of a terrible structure
            score = distanceCA(r,rtarg,w,natomBare)
            auxiliar=seed+iev
            totalTrials=totalTrials+1
            if (MCCheck(error,auxiliar, xbeta , scoreprev, score)) then
               r(1:natomBare)=rprev
            else
               rprev=r(1:natomBare)
               iacc=iacc+1
               scoreprev=score
               tacum=tacum+tacrect
               accepted=accepted+1
            endif
         else
               iacc=iacc+1
               tacum=tacum+tacrect
         endif
            

        enddo
        ! end DO while(tacum.lt.tsnap) -------------------------------------------------

        IF(calcAuxAcc) aux_acc= accepted*1./totalTrials
      ! SELF_STOP
   
         errcont=errcont+1
         errcontLong=errcontLong+1
         IF (errcont > error_evo_size ) errcont=1
         IF (errcontLong > error_long_size ) errcontLong=1
         error_evo(errcont,2)=error
         error_evo(errcont,1)=real(ibloc)
         error_long(errcontLong,2)=error
         error_long(errcontLong,1)=real(ibloc)
         ! DECIDE TO STOP
         slope=abs(lreg(error_evo_size,error_evo))
         slopeLong=abs(lreg(error_long_size,error_long))
        
         if ( (MOD(real(ibloc),5.0) == 0 .and.StartMD) .or.aux_acc<0.2) then
             xbeta= updateXbeta(xbeta,aux_acc,acceptance)
             accepted=1
             totalTrials=1
             calcAuxAcc=.TRUE.
         END IF
 
         call AreInB(natom,natomBare,stepPts,r,numEnB)
        
         fractionDone=Real(numEnB)/Real(contarTodos)
   
         call rmsd(natomBare, r(1:natomBare),rtarg(:), 1, U, center2, center1, error,.FALSE.)
         
         call pdbsuperpos(natomBare, U , center1, center2,r)
         
         ! ways to stop
         !
         IF ( ibloc > 20000 ) stopping=.TRUE.
      !   IF ( fractionDone > 0.95 .and. error < 0.65 ) stopping=.TRUE.
      !   IF ( error < errormin) stopping=.TRUE.
         IF ( slope < DRMSd .and. error < errorAcceptable .and. unForceIt ) stopping=.TRUE.
         ! OJO
         IF ( slopeLong < DRMSd*5.E-2 .and. error < 0.75*error0 .and. ibloc>error_long_size) stopping=.TRUE.
         
         !
         IF (slope < DRMSD .and. error < 0.45*error0 .and. .not.stopping) THEN
           IF(ibloc> error_evo_size) THEN
             ! WRITE(unit_o,*)'FORCE IT'
             IF(.not.unForceIt) error_evo(errcont,2)=1.E10
             call ForceIT(natomBare,r, rtarg,w)
             xbeta=5.00
             scoreprev=1.E10
             unForceIt=.TRUE.
           ENDIF
         END IF
            
         call calcEpot(natom,natomBare, r, stepPts, epotgo, epotfis, epot)
         !
         call writeEnergies(unit_ener, unit_o, temps, epotgo, ekin0, natom,fractionDone,error,&
              aux_acc,xbeta,saco,sum(discardedEner(:,:)),count(stepPts(:,:) % saltos <0), &
              contarSaltos,.TRUE., .TRUE.)          
        
      IF (mod(ibloc,nskip)==0) THEN
          call writeSnapShot(unit_traj, nwritten, r, atoms, res, chain, rnum, natom)
          nwritten=nwritten+1
      END IF

    enddo
    ! end DO iblock-----------------------------------------------------------------
    close(5)
    call cpu_time(tfin)
    write (unit_o, *)
    write (unit_o, '(" T I M I N G S ")')
    write (unit_o, '(" ============= ")')
    write (unit_o, '(" Setup: " f10.2," s")') tsetup - tinit
    write (unit_o, '(" Traj:  ",f10.2," s (",f10.2," ns/h, ",f7.2," ns/day)")') tfin - tsetup, &
    temps / (tfin - tsetup) * 3.6 / 1000., temps / (tfin - tsetup) * 3.6 / 1000. * 24
    write (unit_o, *)
    
    
    WRITE(unit_o,'("INFO")') 
    WRITE(unit_o,'("NATOM ",I5)') natom
    WRITE(unit_o,'("cuts ", 2f6.2 "   rgyr ",f6.2 )') rcut2GO4,rcut1,rgyr
    WRITE(unit_o,'("n pozos 4", I5)')  contarTodos
    WRITE(unit_o,'("TIME ",f10.2, " min")')  (tfin - tinit)/60.
    WRITE(unit_o,'("RMSD ",f8.2, "   % Done ", f4.3)') error, fractionDone
    WRITE(unit_o,*) "slope ", slope , " slopeLong", slopeLong
    WRITE(unit_o,'(" Discarded Energy ", f20.3 )') sum(discardedEner(:,:))
    
    
end
!============================================================================= 