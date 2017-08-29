  ! compile:
  ! > gfortran -c forfun.f90
  ! > gfortran -O3 forfun.o find_chiral_knots.f90 -o find_chiral_knots.out

    program find_chiral_knots

    use forfun      ! module with low level functions
    use fun_parms   ! module with parameters

    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!! edit parameters

    character(len=50), parameter :: filein  = 'Unknot_200_rept_1.erg'
    integer, parameter :: ntot = 10100        ! # saved configs = # lines in erg file 
    
    !!!!!!!!!!!!!!!!!!!!!!!!! end edit
    
    character(len=10) :: knot(100)     ! knot type         
    character(len=5)  :: space         ! space between alex and homf in erg file, needs to be 5 characters

    real(dbl) :: alex(ntot), aread     ! alexander polynomial 
    integer :: alex_int(ntot)          ! nearest integer to (alex * 10000)
    
    character(len=100) :: homf(ntot), hread       ! homfly polynomial
    
    integer :: n, nread                ! configurations
    integer :: nreduced                ! ntot-100, actual number of configs considered

    integer :: h, hmax                   
    character(len=100) :: homf_list(100)  ! list of different homf, indexed by h=1,...,hmax
    integer :: alex_list(100)             ! list of corresponding alex, indexed by h=1,...,amax
    integer :: nh(100)                    ! # homf with given h
    integer :: sh                         ! total # homf, should be = nreduced
    
    integer :: a, amax
    integer :: alex_found(100)            ! list of *different* alex, indexed by a=1,...,amax
    
    logical :: found

    nreduced = ntot - 100

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! read alex and homf values
    open(10,file=trim(filein))         
    do n=1,100                         ! disregard first 100 lines in erg file
       read(10,'(I10,F15.5,A5,A)') nread, aread, space, hread      ! this line determines format of erg file
    enddo
    do n=1,nreduced                    ! read remaining
       read(10,'(I10,F15.5,A5,A)') nread, alex(n), space, homf(n)
       alex_int(n) = nint(alex(n)*10000.0_dbl)  ! nearest integer to (alex * 10000)
    enddo
    close(10)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! make lists homf_list(h), alex_list(h), h=1,...,hmax
    hmax=1                             ! total # of different homf
    homf_list(1) = homf(1)             ! 1st found homf = value for n=1
    alex_list(1) = alex_int(1)         ! 1st found alex = value for n=1
    do n=2, nreduced
       found=.false.
       do h=1,hmax                     ! homf(n) already found?
          if(homf(n).eq.homf_list(h)) then
             found=.true.                                       
          endif
       enddo       
       if(.NOT.found) then                ! homf(n) not found
          hmax = hmax + 1
          homf_list(hmax) = homf(n)       ! add homf(n) to list  
          alex_list(hmax) = alex_int(n)   ! add alex_int(n) to list
       endif
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! find nh(h) = # homf with given h
    nh = 0       ! nh(h) = # homf with given h = 1,...,hmax
    sh = 0       ! sh = total # homf, should be = nreduced
    do h=1,hmax
       do n=1, nreduced
           if(homf(n).eq.homf_list(h)) then
              nh(h) = nh(h) + 1
           endif
       enddo
       sh = sh + nh(h)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! make list alex_found(a), a=1,...,amax
    amax=1                           ! total # of different alex
    alex_found(1) = alex_list(1)     ! 1st found alex = alex_list(1)
    do h=2,hmax
       found=.false.
       do a=1,amax                   ! alex_list(h) already found?
           if(alex_list(h).eq.alex_found(a)) then
             found=.true.                                       
          endif
       enddo       
       if(.NOT.found) then                 ! alex_list(h) not found
          amax = amax + 1
          alex_found(amax) = alex_list(h)  ! add alex_list(h) to list
       endif
    enddo        

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! assign knot types knot(h) to alex_list(h), h=1,...,hmax
    do h=1,hmax 
       if(alex_list(h).eq.10000) then
          knot(h)='0.1'
       else if(alex_list(h).eq.90546) then
          knot(h)='3.1'
       else if(alex_list(h).eq.250910) then
          knot(h)='4.1'
       else if(alex_list(h).eq.254574) then
          knot(h)='5.1'
       else if(alex_list(h).eq.492549) then
          knot(h)='5.2'           
       else if(alex_list(h).eq.508063) then
          knot(h)='7.1'
       else if(alex_list(h).eq.97267) then
          knot(h)='8.19'
       else if(alex_list(h).eq.14892) then
          knot(h)='10.124'
       else if(alex_list(h).eq.77771) then
          knot(h)='10.139'
       else if(alex_list(h).eq.4449) then
          knot(h)='12.242'
       else if(alex_list(h).eq.218658) then
          knot(h)='12.725'            
       else                 ! add more knot types as needed
          knot(h)='unknown'
       endif  
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! output
    write(*,*)
    write(*,*) 'total # different alex:', amax
    write(*,*) 'total # different homf:', hmax 
    write(*,*) 'nreduced, total # homf:', nreduced, sh
    write(*,*)
    write(*,*) '        alex      knot         number    fraction    homf'
    write(*,*) 
    do a=1, amax
       do h=1, hmax
          if(alex_list(h).eq.alex_found(a)) then
             write(*,('(F15.4,A4,A10,I8,F12.5,A)')) &
             real(alex_list(h))/10000_dbl, '    ', knot(h), nh(h), real(nh(h))/nreduced, trim(homf_list(h))
          endif
       enddo
       write(*,*)
    enddo

    end program find_chiral_knots
