!!
!!  Program MONTE: Monte-Carlo simulation of polymer chain
!!  Author: Stefan Giovan
!!
!!  Syntax:  monte parmfile simfile [-erg ergfile] [-traj trajfile] [-alex] [-homf] 
!!
!!  Last Update:
!!
!!      03/23/2015: (SMG) Program is released to group.
!!
!!

        program monte
        use forfun !module with low level functions
        use fwlc_parms !module with parameters
        use fwlc !module with functions
        implicit none

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! VARIABLE DEFINITIONS

        ! default flags (these are changed internally when/if required)
        logical :: doalex = .false.     ! calculate alexander polynomial
        logical :: dohomf = .false.     ! calculate homfly polynomial 
        logical :: dotraj = .false.     ! write tajectory to file
        logical :: doerg = .false.      ! write energy file
        logical :: doev = .false.
        logical :: dotwist = .false.
        logical :: verbose = .false.
        logical :: step = .false.

        integer :: nargs, arglen ! number and length of arguments
        character(len=500) :: simfile ! file of runtime values
        character(len=500) :: trajfile ! file containing saved congifurations
        character(len=500) :: ergfile ! file containing saved energies
        character(len=10)  :: arg ! program arguments

        integer :: nseg        ! number of coordinates (segments)
        integer :: nconfig     ! number of configurations to save in total
        integer :: nstep       ! number of moves to make between saved configurations
        integer :: nupdate=500 ! number of attempted trial moves attempted before
                               ! updating step size limits

        integer :: itrial, ntrial   ! total number of trial moves made, to make
        real(dbl) :: pdisp, prept, pcrank, pmove      ! trial move probability
        real(dbl) :: disp_lim, crank_lim              ! step size limits
        integer :: ncrank_attempt=0, ncrank_accept=0  ! # attempted/accepted trial moves
        integer :: ndisp_attempt=0, ndisp_accept=0
        integer :: nrept_attempt=0, nrept_accept=0
        integer :: nrept_genuine=0                    ! AH: # genuine rept moves
        integer :: nrept_err1=0, nrept_err2=0         ! AH: # 4-seg section too long /
                                                      ! 3-seg and 4-seg sections don't adjust
        real(dbl) :: paccept ! trial move acceptance rate
        integer :: movetype, pta, ptb, ndist

        
        real(dbl), allocatable :: rcurrent(:,:), rtrial(:,:) ! current/trial configs
        real(dbl), allocatable :: rground(:,:), coef0(:,:) ! ground state config
        real(dbl) :: lin(3), vec(3)    ! AH
        real(dbl) :: cbend, cstretch, ctwist, lambda    ! energy parameters
        real(dbl) :: utot, utrial      ! total (effective) energy
        real(dbl) :: ubend, ubtrial    ! bending energy
        real(dbl) :: uint, uinttrial   ! intrinsic bend (blow-away) energy
        real(dbl) :: utwist, utwtrial  ! twist energy
        real(dbl) :: ustretch, ustrial ! stretch energy
        real(dbl) :: dia, dlk, cm(3)   ! diameter, DeltaLK, center of mass

        integer :: nx, nx_trial ! number of knot crossings
        real(dbl), allocatable :: xdat(:,:), xdat_trial(:,:)
        real(dbl) :: wr, wr_trial       ! writhe of current/trial configs
        real(dbl) :: delta_wr           ! AH
        real(dbl) :: gs                 ! AH: square radius of gyration 
        real(dbl) :: aknot, aknot_trial ! alexander polynomials
        integer :: glen, glen_trial     ! gauss code length
        integer :: hlen, hlen_trial     ! homfly length
        character(len=len_gauss_max) :: gcode='', gcode_trial='' ! gauss codes
        character(len=len_homfly_max) :: homf='', homf_trial=''  ! homfly polynomials 

        real(dbl) :: rnd    ! random number [0,1)
        integer :: ierr     ! error identifer
        integer :: juxt     ! AH: flag juxtaposed
        integer :: nrngseed ! length of rngseed
        integer, allocatable :: rngseed(:) 
        integer :: i, j, k  ! generic indices
 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SETUP SIMULATION

        ! Number of command line arguments
        nargs = command_argument_count()

        ! read parameter and simulation file
        if (nargs.lt.1) then
          write(*,*) 'ERROR: You must supply the input file'
          stop
        else
          call get_command_argument(1,simfile)
          write(*,*) 'Input file: ',trim(simfile)
          call read_sim_file()
        endif

        ! get other arguments
        do i=2,nargs
          call get_command_argument(i,arg)
          if (trim(arg).eq.'-v') verbose=.true.
          if (trim(arg).eq.'-s') step=.true.
          if (trim(arg).eq.'-alex') doalex=.true.
          if (trim(arg).eq.'-homf') dohomf=.true.
          if (trim(arg).eq.'-traj') then
            call get_command_argument(i+1,trajfile)
            dotraj=.true.
            write(*,*) 'Saving configurations to: ',trim(trajfile)
          elseif (trim(arg).eq.'-erg') then
            call get_command_argument(i+1,ergfile)
            doerg=.true.
            write(*,*) 'Writing energies to: ',trim(ergfile)
          endif
        enddo
        
        if (dia.gt.0._dbl) then
          write(*,*) 'Excluding chain overlap: diameter/b0=', dia
          doev=.true.
        endif

        if (ctwist.gt.0._dbl) then
          write(*,*) 'Tortional energy enabled: DeltaLk = ', dlk 
          dotwist=.true.
        endif

        ntrial = nstep*nconfig
        if (ntrial.lt.10*nupdate) then
          write(*,*) 'Warning: not enough trial moves will be made'
          write(*,*) ' to update step size limits.'
        endif

        ! calculate energy of initial configuration
        call energy(rcurrent,utot,uint,ubend,ustretch,utwist,wr,nx,xdat)  
  
        ! check for segment overlap
        if (doev) then
          ierr=0
          call check_ev(rcurrent,0,ierr)
          if (ierr.gt.0) then
            write(*,*) 'ERROR: Initial configuration violates EV'
            stop
          endif
        endif
        
        ! calculate knot-type of initial configuration
        if (doalex) then
          aknot = alex(nx,xdat)
          write(*,*) 'Initial Alexander Product: ', aknot
        endif

        if (dohomf) then
          call gauss_code(nx,xdat(:,1:nx),glen,gcode)
          call homfly(glen,gcode(1:glen),hlen,homf,.true.) 
          write(*,*) 'Initial homfly: ', homf(1:hlen)
        endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MONTE CARLO ITERATION

        do itrial = 1,ntrial
          if (step) read(*,*)
          ierr=0

          ! UPDATE TRIAL MOVE STEP SIZES
          ! AH: reptation moves don't have step size, nothing to update
          
          ! update crankshaft limit
          if (ncrank_attempt.ge.nupdate) then
            paccept = dble(ncrank_accept)/dble(ncrank_attempt)
            if (paccept>0.7_dbl) then
              crank_lim = crank_lim*1.15_dbl
              !if (verbose) then
              !  write(*,*) 'Crankshaft Limit increased to: ', crank_lim
              !endif
            elseif (paccept<0.5_dbl) then
              crank_lim = crank_lim*0.95_dbl
              !if (verbose) then
              !  write(*,*) 'Crankshaft Limit decreased to: ', crank_lim
              !endif
            endif
            ncrank_accept = 0
            ncrank_attempt = 0
          endif

          ! update displacement limit
          if (ndisp_attempt.ge.nupdate) then
            paccept = dble(ndisp_accept)/dble(ndisp_attempt)
            if (paccept>0.7_dbl) then
              disp_lim = disp_lim*1.15_dbl
              !if (verbose) then
              !  write(*,*) 'Displacement Limit increased to: ', disp_lim
              !endif
            elseif (paccept<0.5_dbl) then
              disp_lim = disp_lim*0.95_dbl
              !if (verbose) then
              !  write(*,*) 'Displacement Limit decreased to: ', disp_lim
              !endif
            endif
            ndisp_accept = 0
            ndisp_attempt = 0
          endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
! UPDATE OUTPUT FILES

          if (mod(itrial,nstep).eq.0) then

            ! recenter
            cm = sum(rcurrent,2)/dble(nseg)
            do i = 1,nseg
              rcurrent(:,i) = rcurrent(:,i)-cm
            enddo

            ! update simfile
            call write_sim_file()

            ! AH: calculate square radius of gyration
            !gs = 0.0_dbl
            !do i = 1,nseg
            !   gs = gs + sum(rcurrent(:,i)**2)
            !enddo
            !gs = gs / dble(nseg)
            
            ! update simfile
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            ! AH: Check if segments are juxtaposed and return knot type of passed configuration
            ! hairpin and straight G-segments are implemented in subroutine wlc_energy in fwlc.f90

            ! juxt=0
            ! call check_juxt_straight(rcurrent,juxt,vec)  
            ! if (juxt.gt.0) then           ! if juxtaposed, define
            !  do j=2,nseg                  ! trial conformation := passed conformation
            !    rtrial(:,j)=rcurrent(:,j)
            !  enddo                     
            !  
            !  ! for straight G-segment:
            !  rtrial(:,1) = rcurrent(:,1) + vec  ! define point 1 of passed conformation
            ! 
            !  ! for hairpin G-segment: 
            !  lin(:) = (rcurrent(:,3) - rcurrent(:,nseg-1)) / 4.0_dbl
            !  rtrial(:,nseg) = rtrial(:,nseg-1) + lin(:)  ! define N, 1, 2 of passed conformation
            !  rtrial(:,1) = rtrial(:,nseg) + lin(:)
            !  rtrial(:,2) = rtrial(:,1) + lin(:)
              
              !nx_trial=0
              !xdat_trial=0
              !wr_trial=0
              !call find_crossings(nseg,rtrial,nx_trial,xdat_trial)
              !aknot_trial = alex(nx_trial,xdat_trial)
              !wr_trial = writhe(nseg,rtrial,nx_trial,xdat_trial)                 
              !delta_wr = wr_trial - wr
              
              !if (verbose) then
              !   write(*,*) 'Iteration Juxtaposed: ', itrial/nstep, &
              !        'Alexander Product of Current, Passed: ', &
              !        aknot, aknot_trial, &
              !        'Delta Wr: ', delta_wr
              !endif

              ! save juxtaposed  configuration to trajectory file

              nx=0
              xdat=0
              wr=0
              call find_crossings(nseg,rcurrent,nx,xdat)
              aknot = alex(nx,xdat)
              wr = writhe(nseg,rcurrent,nx,xdat)

              call gauss_code(nx,xdat(:,1:nx),glen,gcode)
              call homfly(glen,gcode(1:glen),hlen,homf,.true.) 

              if (verbose) then
                 !write(*,*) 'itrial/nstep, aknot, nx, rept moves: attempted, no-fit, no-adjust, genuine, accepted', &
                 !  itrial/nstep, aknot, nx, nrept_attempt, nrept_err1, nrept_err2, nrept_genuine, nrept_accept 
                 write(*,*) 'itrial/nstep, crank_lim, disp_lim, rept:genuine, rept:accepted, alex, homf', &
                      itrial/nstep, crank_lim, disp_lim, nrept_genuine, nrept_accept, aknot, homf(1:hlen)  
              endif

              ! write trajectory file
              if (dotraj) then
                open(14,file=trim(trajfile),position='append')
                write(14,'(3ES25.15)') (rcurrent(:,i),i=1,nseg)
                write(14,*)
                !!!!!!!
                !write(14,'(3ES25.15)') (rtrial(:,i),i=1,nseg)
                !write(14,*)
                !write(14,*)
                !write(14,*)
                !!!!!!!
                close(14)
              endif 
                 
              ! write data file
              if (doerg) then
                open(15,file=trim(ergfile),position='append')
                !if (dotwist) then
                !  write(15,'(7ES15.5)') utot,ubend,uint,ustretch,utwist,real(juxt),aknot_trial
                !else
                write(15,'(I10,F20.12,A,A)') itrial/nstep, aknot, '     ', homf(1:hlen)
                ! space b/w aknot and homf(1:hlen) should be exactly 5 blanks 
                ! to make output compatible with find_chiral_knots.f90
                !endif
                close(15)
              endif

            !endif     ! endif (juxt.gt.0) 
            
          endif        ! endif (mod(itrial,nstep).eq.0)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MAKE A TRIAL MOVE

 800      rtrial = rcurrent
          call random_number(pmove) ! pick a trial move type
          if (pmove.lt.pcrank) then
            ! Do crankshaft
            movetype = 1
            ncrank_attempt = ncrank_attempt+1
            pta = randi(1,nseg)
            ndist = randi(2,nseg/2)
            ptb = pta+ndist
            if (ptb.gt.nseg) ptb=ptb-nseg
            call crankshaft(nseg,rtrial,pta,ptb,randr(-crank_lim,crank_lim))
          elseif (pmove.lt.pcrank+pdisp) then 
            ! Do displacement
            movetype = 2
            ndisp_attempt = ndisp_attempt+1
            pta = randi(1,nseg)
            ndist = randi(2,nseg/2)
            ptb = pta+ndist
            if (ptb.gt.nseg) ptb=ptb-nseg
            call displace(nseg,rtrial,pta,ptb,randr(0._dbl,disp_lim)) 
          elseif (pmove.lt.pcrank+pdisp+prept) then
            ! AH: Do reptation
            movetype = 3
            nrept_attempt = nrept_attempt+1
            pta = randi(1,nseg)                ! pta = is, start of original 3-seg section
            ndist = randi(6,nseg-6)
            ptb = pta+ndist                    ! ptb = il, start of original 4-seg section
            if (ptb.gt.nseg) ptb=ptb-nseg
            ierr=0
            call reptation(nseg,rtrial,pta,ptb,ierr)  ! new subroutine in fwlc.f90
            if (ierr.ne.0) then                ! error => attempted rept move not genuine      
               if (ierr.eq.1) then
                  nrept_err1 = nrept_err1 + 1  ! 4-seg section too long for 3-seg section
               else
                  nrept_err2 = nrept_err2 + 1  ! 3-seg and 4-seg section lengths don't adjust
               endif   
               goto 800                        ! make another trial move 
            endif
            nrept_genuine = nrept_genuine + 1  ! no error => attempted rept move is genuine 
          else   ! should never happen
            write(*,*) 'ERROR: No trial move chosen!'
            stop
          endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ACCEPT/REJECT TRIAL CONFIGURATION

          ! calculate energy of trial configuration
          call energy(rtrial,utrial,uinttrial,ubtrial,&
          &ustrial,utwtrial,wr_trial,nx_trial,xdat_trial)
          ! metropolis criterion        
          call random_number(rnd)
          if (exp(utot-utrial).lt.rnd) then
            !if (verbose) then
            !  write(*,*) 'Iter:', itrial, ' Movetype:', movetype, ' Failed (Energy)'
            !endif
            cycle  ! next Monte Carlo iteration with previous 'current'
          endif
            
          ! Check for segment overlap (ev)
          if (doev) then
            ierr=0
            call check_ev(rtrial,movetype,ierr)
            if (ierr.gt.0) then
              !if (verbose) then
              !  write(*,*) 'Iter:', itrial, ' Movetype:', movetype, ' Failed (EV)'
              !endif
              cycle ! next Monte Carlo iteration with previous 'current'
            endif
          endif

          ! check for changes in knot-type using Alexander polynomial  
          if (doalex) then
            aknot_trial = alex(nx_trial,xdat_trial(:,1:nx_trial))
            if (abs(aknot-aknot_trial).gt.0.001_dbl) then
              ! if (verbose) then
              !   write(*,*) 'Iter:', itrial, ' Movetype:', movetype, ' Failed (ALEX)'
              ! endif
              cycle ! next Monte Carlo iteration with previous 'current'
            endif
          endif

          ! check for changes in knot-type using HOMFLY polynomial
          if (dohomf) then
            ierr=0 
            call check_homf(ierr)
            if (ierr.gt.0) then
              !if (verbose) then
              !  write(*,*) 'Iter:', itrial, ' Movetype:', movetype, ' Failed (HOMF)'
              !endif
              cycle ! next Monte Carlo iteration with previous 'current'
            endif
          endif

          ! TRIAL MOVE HAS BEEN ACCEPTED - SAVE AS NEW CURRENT
          ! if (verbose) then
          !   write(*,*) 'Iter:', itrial, ' Movetype:', movetype, ' PASSED'
          ! endif
          if (movetype.eq.1) then
            ncrank_accept = ncrank_accept+1
          elseif (movetype.eq.2) then
            ndisp_accept = ndisp_accept+1
          elseif (movetype.eq.3) then
            nrept_accept = nrept_accept+1
          endif
          rcurrent = rtrial
          utot = utrial
          ubend = ubtrial
          uint = uinttrial
          ustretch = ustrial
          utwist = utwtrial
          wr = wr_trial
          nx = nx_trial
          xdat = xdat_trial
          glen = glen_trial
          gcode = gcode_trial
            
        end do  ! next Monte Carlo iteration with accepted trial move as new 'current'
        
! END MONTE CARLO ITERATION        
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! SUBROUTINES

        contains

!------------------------------------------------------------------------------

        subroutine read_sim_file()
        implicit none
        integer :: nk, krand(50), sec
        real(dbl) :: psum
        open(13,file=trim(simfile))
        read(13,*) nseg                         ! line 1

        ! allocate
        allocate(rground(3,nseg), rcurrent(3,nseg), rtrial(3,nseg))
        allocate(coef0(6,nseg))
        allocate(xdat(4,nseg*nseg),xdat_trial(4,nseg*nseg))

        read(13,*) nconfig, nstep               ! 2
        read(13,*) dlk, dia                     ! 3
        read(13,*) cbend, cstretch, ctwist      ! 4
        read(13,*) lambda                       ! 5
        read(13,*) pcrank, pdisp, prept         ! 6
        read(13,*) ! blank                      ! 7
        read(13,*) ! blank                      ! 8
        read(13,*) (rground(:,i), i=1,nseg)     ! 9->nseg+8
        read(13,*) ! blank                      ! nseg+9
        read(13,*) ! blank                      ! nseg+10
        read(13,*) nk, krand(1:nk)              ! nseg+11
        read(13,*) ! blank
        read(13,*) crank_lim, disp_lim
        read(13,*) ! blank
        read(13,*) (rcurrent(:,i), i=1,nseg)
        close(13)

        ! set RNG seed
        call random_seed(size=nrngseed)
        allocate(rngseed(nrngseed))
        if (nk.eq.1.and.krand(1).eq.-1) then        
          call system_clock(sec)
          rngseed(:) = abs( mod((sec*181)*((getpid()-83)*359), 104729) )
        elseif (nk.eq.1.and.krand(1).ne.-1) then
          rngseed(:) = krand(1)
        elseif (nk.eq.nrngseed) then 
          rngseed = krand(1:nk)
        elseif (nk.ne.nrngseed) then
          write(*,*) 'RNG seed length has changed for this system'
          write(*,*) 'Change the seed to a scalar and restart'
          stop
        endif
        call random_seed(put=rngseed)

        ! normalize trial move probabalities
        psum = pcrank+pdisp+prept ! sum of trial move probabilities
        pcrank = pcrank/psum ! normalized
        pdisp = pdisp/psum ! normalized
        prept = prept/psum ! normalized

        ! get ground state coef
        call get_coef(nseg,rground,coef0)

        end subroutine read_sim_file

!------------------------------------------------------------------------------

        subroutine write_sim_file()
        implicit none
        integer :: i
        call random_seed(get=rngseed)
        open(13,file=trim(simfile))
        do i = 1,nseg+10
          read(13,*)
        enddo
        write(13,*) nrngseed, rngseed
        write(13,*) '! Trial Move Displacement Limits [crank_lim, disp_lim]'
        write(13,'(2ES25.15)') crank_lim, disp_lim
        write(13,*) '! Current Configuration'
        write(13,'(3ES25.15)') (rcurrent(:,i), i=1,nseg)
        close(13)
        end subroutine write_sim_file

!------------------------------------------------------------------------------

        subroutine energy(r,u,ui,ub,us,utw,wr,nx,xdat)
        implicit none
        real(dbl), intent(in) :: r(3,nseg)
        real(dbl), intent(out) :: u, ui, ub, us, utw
        integer, intent(out) :: nx 
        real(dbl), intent(out) :: wr, xdat(4,nseg*nseg)
        nx=0
        xdat=0
        wr=0
        if (doalex.or.dohomf.or.dotwist) then
          call find_crossings(nseg,r,nx,xdat)
        endif
        if (dotwist) then
          wr = writhe(nseg,r,nx,xdat)
        endif
        call wlc_energy(nseg,r,coef0,cbend,cstretch,ctwist,lambda,dlk,wr,&
       &u,ui,ub,us,utw)
        end subroutine energy

!------------------------------------------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine check_juxt(r,juxt)   ! AH: return juxt=1 if chain passes 
                                        ! through hairpin-triangle N-1,1,3
        implicit none
        real(dbl), intent(in) :: r(3,nseg)
        integer, intent(out) :: juxt
        integer :: i
        real(dbl) :: p1(3), p2(3), p3(3), p4(3), p5(3)
 
        juxt=0
        
        p1=r(:,nseg-1)
        p2=r(:,1)
        p3=r(:,3)
        do i=4,nseg-3  ! check segments i for passing through triangle
           p4=r(:,i)
           p5=r(:,i+1)
           if (segment_intersect_triangle(p1,p2,p3,p4,p5)==1) then
              juxt=1   ! segment i juxtaposes
              return
           endif
        enddo   
        end subroutine check_juxt
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        subroutine check_juxt_straight(r,juxt,vec)  ! AH: return juxt=1 if chain passes
                                                    ! through triangle N, 2, 1+vec 
        implicit none
        real(dbl), intent(in) :: r(3,nseg)
        integer, intent(out) :: juxt
        real(dbl), intent(out) :: vec(3)            ! random vector perp to vector 1 -> 2
                                                    ! of length sqrt(3) 
        integer :: i
        real(dbl) :: p1(3), p2(3), p3(3), p4(3), p5(3)
        real(dbl) :: sqrt3
        real(dbl) :: b1(3), rv(3)
        
        juxt=0

        sqrt3 = 1.73205_dbl
        b1 = r(:,2) - r(:,1)      ! vector 1 -> 2
        rv = randvec()            ! randomly oriented unit vector
        vec=cross(b1,rv)          ! vec = random vector perp b1
        if (sqrt(sum(vec(:)**2)).lt.0.001_dbl) then
           return          
        endif
        vec(:) = sqrt3 * vec(:) / sqrt(sum(vec(:)**2)) ! length of vec = sqrt(3) so that
                                                       ! N, 2, 1+vec = equal sided triangle  
        p1=r(:,nseg)
        p2=r(:,2)
        p3=r(:,1) + vec
        do i=3,nseg-2  ! check segments i for passing through triangle
           p4=r(:,i)
           p5=r(:,i+1)
           if (segment_intersect_triangle(p1,p2,p3,p4,p5)==1) then
              juxt=1   ! segment i juxtaposes
           endif   
        enddo   
        end subroutine check_juxt_straight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine check_ev(r,movetype,ierr)
        implicit none
        integer, intent(in) :: movetype
        real(dbl), intent(in) :: r(3,nseg)
        integer, intent(out) :: ierr
        integer :: i, j
        real(dbl) :: p1(3), p2(3), p3(3), p4(3)
        real(dbl) :: m1(3), m2(3), l1, l2, dist, cutoff
        ierr=0
        do i=1,nseg-2  ! check segments i, j for ev
          do j=i+2,nseg
            if (i.eq.1.and.j.eq.nseg) cycle
            ! TODO: Review possible conditions to avoid unncessary EV checks
!            if (movetype.eq.1) then ! crankshaft move
!              if (i.lt.pta.and.i.ge.ptb) cycle ! i not in [pta,ptb]
!              if (j.ge.pta.and.j.lt.ptb) cycle ! j in [pta,ptb]
!            endif
!            if (movetype.eq.2) then ! stretching move
!              if (i.lt.pta.and.i.ge.ptb) cycle ! i not in [pta,ptb]
!              if (j.ge.pta+1.and.j.lt.ptb-1) cycle ! j in [pta+1,ptb-1]
!            endif
            p1=r(:,i)
            p2=r(:,i+1)
            p3=r(:,j)
            if (j.eq.nseg) then
              p4=r(:,1)
            else
              p4=r(:,j+1)
            endif
            m1=p1+0.5_dbl*(p2-p1) !midpoint of segment 1
            m2=p3+0.5_dbl*(p4-p3) !midpoint of segment 2
            l1 = sqrt(sum((p2-p1)**2)) ! length of segment 1
            l2 = sqrt(sum((p4-p3)**2)) ! length of segment 2
            dist = sqrt(sum((m2-m1)**2)) ! distance between midpoints
            cutoff = 0.5_dbl*(l1+l2) + dia
            if (dist.gt.cutoff) cycle    ! segments cannot overlap 
            if (segdist(p1,p2,p3,p4).lt.dia) then
              ierr=1 ! Segments overlap
              return
            endif
          enddo
        enddo
        end subroutine check_ev

!------------------------------------------------------------------------------
        
        subroutine check_homf(ierr)
        implicit none
        integer, intent(out) :: ierr
        real(dbl) :: xdat_red(4,nseg*nseg)
        character(len=len_gauss_max) :: gcode_red
        integer :: nx_red, glen_red
        ierr=0
        ! calculate gauss code
        call gauss_code(nx_trial,xdat_trial(:,1:nx_trial),glen_trial,gcode_trial)
        ! if gauss code is unchanged then so is knot-type, no need to check
        ! but... 
        if (gcode_trial.ne.gcode) then
          ! check homf
          if (nx_trial.lt.15) then
            ! calculate homfly of the trial configuration
            call homfly(glen_trial,gcode_trial(1:glen_trial),hlen_trial,homf_trial)
          else ! attempt to reduce number of crossings
            call reduce_crossings(nseg,rtrial,nx_red,xdat_red)
            call gauss_code(nx_red,xdat_red(:,1:nx_red),glen_red,gcode_red)
            call homfly(glen_red,gcode_red(1:glen_red),hlen_trial,homf_trial)
          endif
          if (homf_trial.ne.homf) ierr=1 ! back to top
        endif
        end subroutine check_homf

!------------------------------------------------------------------------------
        end program monte        
