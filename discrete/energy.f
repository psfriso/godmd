!===============================================================================
   subroutine actualizeStepPot(stepPts,natom,natomBare, nblist, xsum)
   use stepPotentials
   use intList
   integer, intent(IN) :: natom,natomBare
   type(stepPotInt), intent(IN) :: stepPts(natom,natom)
   type(intpList), intent(INOUT) :: nblist(natom)
   real, intent(IN) :: xsum
   integer i, k,j
  
   do i = 1,natomBare-1
   
       do k=1,nblist(i)%nats 
       
           j=nblist(i)%idata(k)%patnum

         if ( stepPts(i,j)%active .and. stepPts(i,j)%nstep ==4 ) THEN

             nblist(i)%iData(k) %stepPt =  stepPts(i,j)
           
       END IF
  
      enddo
  enddo
  
  
   end subroutine actualizeStepPot
!===============================================================================
!===============================================================================
   subroutine activateStepPot(stepPts, r,LIGAND ,dist_INI,dist_TARG,rcut2GO2,rcut2GO4,natom,natomBare, nblist, xsum)
   use stepPotentials
   use geometryDP
   use intList
   integer, intent(IN) :: natom,natomBare
   real*8, intent(IN) :: rcut2GO2,rcut2GO4
   REAL*8, INTENT (IN) :: dist_INI(natom,natom), dist_TARG(natomBare,natomBare)
   type(pointDP), intent(IN) :: r(natom)
   type(stepPotInt), intent(INOUT) :: stepPts(natom,natom)
   type(intpList), intent(INOUT) :: nblist(natom)
   LOGICAL, INTENT(IN) :: LIGAND
   real, intent(IN) :: xsum
   integer i, j
   real*8 :: rij_INI,rij_TARG
   type(pointDP) rj
   logical :: IsActive
!
   where (stepPts%active)
      stepPts%active = .false.
   end where
   
   nblist%nats = 0
   IsActive=.FALSE.
  ! write(*,*) "here 0"

 IF(LIGAND) THEN
 
          
   do j = 2,natom
    !  rj = r(j) ! intentem millorar cache hits a r(i)
      do i = 1,j-1
    !      write(*,*)stepPts(i,j)%tipInt
       secstructlig: if (stepPts(i,j)%tipInt == SSEC ) THEN
   !    write(*,*) "here"
! inline 
         rij_INI = dist_INI(i,j)**2!(r(i)%x-rj%x)**2+(r(i)%y-rj%y)**2+(r(i)%z-rj%z)**2
         rij_TARG =dist_TARG(i,j)**2
         
       nstepsl:  IF( stepPts(i,j)%nstep == 2) THEN
             
           isactive1l:  IF(rij_INI< rcut2GO2) THEN
                stepPts(i,j)%active = .true.
                nblist(i)%nats = nblist(i)%nats + 1
                nblist(j)%nats = nblist(j)%nats + 1
                nblist(i)%iData(nblist(i)%nats) = intData(j, nblist(j)%nats, stepPts(i,j), xsum, 1.e15, 0.)
                nblist(j)%iData(nblist(j)%nats) = intData(i, nblist(i)%nats, stepPts(i,j), xsum, 1.e15, 0.)
             END IF isactive1l
   
           ELSEIF  ( stepPts(i,j)%nstep == 4) THEN
     
    isactive21l:   IF (MOD(real(i+j),2.0)==0 ) THEN
                 
    isactive22l:     IF ((rij_INI<rcut2GO4 ).OR.(rij_TARG<rcut2GO4)) THEN
                                
                     stepPts(i,j)%active = .true.
                     nblist(i)%nats = nblist(i)%nats + 1
                     nblist(j)%nats = nblist(j)%nats + 1
                     nblist(i)%iData(nblist(i)%nats) = intData(j, nblist(j)%nats, stepPts(i,j), xsum, 1.e15, 0.)
                     nblist(j)%iData(nblist(j)%nats) = intData(i, nblist(i)%nats, stepPts(i,j), xsum, 1.e15, 0.)
                END IF isactive22l
        
                 END IF isactive21l
                 
                
         ENDIF nstepsl

         ELSEIF (stepPts(i,j)%tipInt==10)THEN
      !       write(*,*) " ACTIVO ", i,j
       stepPts(i,j)%active = .true.
       nblist(i)%nats = nblist(i)%nats + 1
       nblist(j)%nats = nblist(j)%nats + 1
       nblist(i)%iData(nblist(i)%nats) = intData(j, nblist(j)%nats, stepPts(i,j), xsum, 1.e15, 0.)
       nblist(j)%iData(nblist(j)%nats) = intData(i, nblist(i)%nats, stepPts(i,j), xsum, 1.e15, 0.)
             
      END IF secstructlig
       
  
      enddo
  enddo
 
 ELSE
     
   do j = 2,natomBare
    !  rj = r(j) ! intentem millorar cache hits a r(i)
      do i = 1,j-1
    !      write(*,*)stepPts(i,j)%tipInt
       secstruct: if (stepPts(i,j)%tipInt == SSEC ) THEN
   !    write(*,*) "here"
! inline 
         rij_INI = dist_INI(i,j)**2!(r(i)%x-rj%x)**2+(r(i)%y-rj%y)**2+(r(i)%z-rj%z)**2
         rij_TARG =dist_TARG(i,j)**2
         
       nsteps:  IF( stepPts(i,j)%nstep == 2) THEN
             
           isactive1:  IF(rij_INI< rcut2GO2) THEN
                stepPts(i,j)%active = .true.
                nblist(i)%nats = nblist(i)%nats + 1
                nblist(j)%nats = nblist(j)%nats + 1
                nblist(i)%iData(nblist(i)%nats) = intData(j, nblist(j)%nats, stepPts(i,j), xsum, 1.e15, 0.)
                nblist(j)%iData(nblist(j)%nats) = intData(i, nblist(i)%nats, stepPts(i,j), xsum, 1.e15, 0.)
             END IF isactive1
   
           ELSEIF  ( stepPts(i,j)%nstep == 4) THEN
     
    isactive21:   IF (MOD(real(i+j),2.0)==0 ) THEN
                 
    isactive22:     IF ((rij_INI<rcut2GO4 ).OR.(rij_TARG<rcut2GO4)) THEN
                                
                     stepPts(i,j)%active = .true.
                     nblist(i)%nats = nblist(i)%nats + 1
                     nblist(j)%nats = nblist(j)%nats + 1
                     nblist(i)%iData(nblist(i)%nats) = intData(j, nblist(j)%nats, stepPts(i,j), xsum, 1.e15, 0.)
                     nblist(j)%iData(nblist(j)%nats) = intData(i, nblist(i)%nats, stepPts(i,j), xsum, 1.e15, 0.)
                END IF isactive22
        
                 END IF isactive21
                 
                
         ENDIF nsteps

      END IF secstruct
  
      enddo
  enddo
  
  ENDIF 
  
 end subroutine activateStepPot
!===============================================================================
 
 subroutine thermalize(seed, iev, natom, TEMP, xmassa, v, xm, ekin)
 use geometryDP
 use ran_mod
   integer, intent(IN) :: iev, seed, natom
   real, intent(IN) ::  xmassa, temp, xm
   type(pointDP), intent(OUT) :: v(natom)
   real, intent(OUT) :: ekin
   integer i, kk,NegSeed
   real calcEkin
   type(pointDP) vcm
   real*8 fi, sto
!
   negSeed=seed
   fi=ran1(negSeed)
  ! fi=ran1(seed)
   !   write(*,*)"iev desde thermalize", iev
   vcm = pointDP(0., 0., 0.)
   do i = 1,natom
      kk = (seed + 3 * i + 1 + iev)
!   write(*,*)"kk desde thermalize", kk
      fi=ran1(kk)
      v(i)%x = fi - 0.5
!   write(*,*)"kk desde thermalize", kk
       kk = (seed + 3 * i + 1 + iev)
      fi=ran1(kk)
      v(i)%y = fi - 0.5
!   write(*,*)"kk desde thermalize", kk
      kk = (seed + 3 * i + 1 + iev)
      fi=ran1(kk)
      v(i)%z = fi - 0.5
!
      vcm = vcm + xm*1.d0 * v(i)
   enddo
   vcm = (1.d0/xmassa) * vcm
 !  write(*,*)" vel CM", vcm
   do i = 1,natom
      v(i) = v(i) - vcm
   enddo
   sto = sqrt(1.5 * natom * TEMP / calcEkin(v, xm, natom))
 !  write(*,*) sto
   do i = 1,natom
      v(i) =  sto * v(i)
 !     write(*,*) v(i)

   enddo

   ekin = 1.5 * natom * TEMP ! No cal recalcularla 
   
 end subroutine thermalize
!========================================================================
!========================================================================
   subroutine calcEpot(natom,natomBare, r, stepPts, epotgo, epotfis, epot)
        use geometryDP
        use stepPotentials
        use constants

        integer, intent(IN) :: natom,natomBare
        type(pointDP), intent(IN) :: r(natom)
        type(stepPotInt), intent(IN) :: stepPts(natom, natom)
        real, intent(OUT) :: epotgo, epotfis
        real epot(MAXTIPINT)
        real dist
        integer i, j, k
        !PENDENT TREBALLAR AMB NBLIST
        epotgo = 0.
        epotfis = 0.
        epot = 0.
        do j = 2, natomBare
            do i = 1, j - 1
                if (stepPts(i, j) % active) then
                    dist = sqrt(real(calcDist2DP(r(i), r(j))))
                    k = stepPts(i, j) % nstep
                    do while ((k .gt. 1).and.dist .lt. stepPts(i, j) % step(k) % r)
                        epot(stepPts(i, j) % tipInt) = epot(stepPts(i, j) % tipInt) - stepPts(i, j) % step(k) % e / FACTE
                        k = k - 1
                    enddo
                    if (dist .lt. stepPts(i, j) % step(k) % r) &
                    epot(stepPts(i, j) % tipInt) = epot(stepPts(i, j) % tipInt) - stepPts(i, j) % step(k) % e / FACTE
                endif
            enddo
       enddo
       epotgo = epot(SSEC)
       epotfis = sum(epot) - epotgo
   end subroutine calcEpot
    !========================================================================
 pure function  calcEkin (v, xm, natom) result (ekin)
 use geometryDP
   integer, intent(IN) :: natom
   real ekin
   type(pointDP), intent(IN) :: v(natom)
   real, intent(IN) :: xm
   real*8, parameter :: a2 = 1.d-20*1.d30
   integer i
   ekin = 0.
   do i = 1,natom
      ekin = ekin + 0.5 * xm* real(A2 * dotDP(v(i), v(i)))
   enddo
 end function calcEkin
!===============================================================================

