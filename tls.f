!
!   Two-level system: 
!      1) circular graphs
!      2) fixed \beta
!      3) local add/drop updates only 
!      4) uniform distributions in updates
!
!
	use mumbers         ! random # generator
	implicit none

!================ physical param =============================

	real*8  :: t       ! off-diagonal amplitude, in units of \xi=1
	real*8  :: beta    ! inverse temp, same units

!================ configuration =================================
	integer,parameter :: nmnm_max = 512
	integer                        :: nmnm         ! total number of kinks
      INTEGER, DIMENSION(nmnm_max) :: namelist, numnam

	integer :: ass_l(nmnm_max), ass_r(nmnm_max)    ! links
	integer :: state(nmnm_max)                     ! state=1 for |1>  and state=-1 for |2>
	real*8  :: tau(nmnm_max)                       ! Flat interval duration in \tau
	integer :: masha                               ! 1st (leftmost) flat interval

!=============== statistics =========================================
	integer, parameter :: b_n_max=1000    ! Max # of blocks
	integer  :: Z_b                      ! # block size

	integer :: b_n, i_b    ! filled block nbr. & current nbr. of measurements 
	                       ! in the block into the (b_n+1)-th block

	real*8 :: Z            ! MC denominator aka ' partition function'
	real*8 :: ene, ene_stat(b_n_max)   ! energy

!================ flow control =====================================
	integer :: ans1,ans2  ! read or init cnf/stat

! Update manager probabilities
      real*8, dimension(1:5) :: prob

      real*8 :: step, step_p, step_w, step_m, step_t
!     step    counts MC steps
!     step_p  the number of steps between printouts
!     step_w  the number of steps between writing to disk
!     step_m  the number of steps between measurements
!     step_t  the number of steps for thermolization

      real*8 :: i_p, i_w, i_m, i_t
!     Counters for printing, writing, measuring, thermolization, det recalculation

! structure output
      real*8, dimension(0:nmnm_max) :: nmnm_distr    ! kink # distribution

! update addr/acpt counters
	real*8 :: c_a, a_a, c_d,a_d
	character*4 :: status

	logical :: therm_flag

!================== disposable variables ================================
	integer :: j
	real*8 :: r



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++ Done with variables  +++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!+=====================  initialization ==================================

! read parameter file
	open(1,file='par')
	read(1,*)beta
	read(1,*)t
	read(1,*)ans1  ! read/init cnf
	read(1,*)ans2  ! read/init stat
	read(1,*)step_p  ! step for printing
	read(1,*)step_w  ! step for writing
	read(1,*)step_m  ! step for measuring
	read(1,*)step_t  ! step for thermalizations (mln)
	     step_t = step_t*1d6
	close(1)

! init rndm
	call init_rndm(4836,2740)

! update manager probabilities
	prob=0.d0
	prob(1)=0.5d0   ! add 
	prob(2)=1.d0    ! drop

! configuration & thermalization
	if(ans1==0)then; call init_cnf
	else;            call rd_cnf
	endif

	if(step_t>0)call therm

! statistics
	if(ans2==0)then; call init_stat;
	else;            call rd_stat
	endif

      i_p=0.d0; i_w=0.d0; i_m=0.d0; step=0.d0;
	nmnm_distr=0.d0


!============================= Main MC loop ==============================
	do
	
          step=step+1.d0
          i_p=i_p+1.d0; i_w=i_w+1.d0; i_m=i_m+1.d0


!--- Update :)
           r = rndm()
           if     (r<prob(1))then;    call add
           else if(r<prob(2))then;    call drop
           else; print*,'am I nuts or what???'; stop
           endif

!------ this checks configuration; only useful when developing the code
!	     call check
!----------------------------------------------------------------------

           if(nmnm<nmnm_max) nmnm_distr(nmnm) = nmnm_distr(nmnm) + 1.d0

           if (i_m == step_m)  then; i_m=0.d0; call MEASURE; end if
           if (i_p == step_p)  then; i_p=0.d0; call PRNT;    end if
           if (i_w == step_w)  then; i_w=0.d0; call WRT;     end if

	enddo



!******************************************** subroutines below ONLY
	contains

!===================== add a flat interval: =======================
!
!                    name      igor
!         _________...........______________       
!                 |           |
!old=     ........|___________|.............
!                             tb          
!
! [length of an interval is shown at its right-hand end]
!
!               name vova lesha igor
!       __________.._____..._______________
!new=            |  |   |   |
!       .........|__|...|___|..............
!                   t1  t2  tb-(t1+t2)
!
	subroutine add
	integer :: vova,lesha,igor,name,state_new,state_name,j
	real*8  :: ratio,tb,t1,t2,delta_t


	status='add_';	c_a = c_a + 1.d0

	if(nmnm==nmnm_max)then; print*,'add: nmnm > nmnm_max'; stop
	endif

! select an interval
	j = rn(nmnm); name=namelist(j)
	tb=tau(name); 
	state_name = state(name)

! select times
	t1 = tb*rndm()
      t2 = (tb-t1)*rndm()

! state of 'new' is opposite to the state of 'name'
	if(state_name==1)then; state_new=-1
	else; state_new=1
	endif

! ratio
	ratio = t*t*exp(-(state_new-state_name)*t2)   ! conf. weight
	ratio = ratio * tb*(tb-t1)                    ! prob. distribution for \tau-s
	ratio = ratio * nmnm/(nmnm+2)                 ! context factor

! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

! accepted... update configuration

	call GetName(vova); call GetName(lesha)

! links
	igor=ass_r(name); 
	ass_r(lesha)=igor; ass_r(vova)=lesha; ass_r(name)=vova	
	ass_l(igor)=lesha; ass_l(lesha)=vova; ass_l(vova)=name

! tau-s
	tau(name)=t1; tau(vova)=t2; tau(lesha)=tb-t1-t2

! states
	state(vova)=state_new; state(lesha)=state_name

	a_a = a_a + 1.d0

	end subroutine add


!================ drop a flat interval: ========================
!
!               name vova lesha igor
!       __________.._____..._______________
!old=            |  |   |   |
!       .........|__|...|___|..............
!                   t1  t2  tb
!
!                   name      igor
!        _________...........______________
!                |           |
!new     ........|___________|.............
!                            tb+t1+t2
!
	subroutine drop
	real*8 :: tb,t1,t2,ratio
	integer :: lesha,vova,igor,name,state_vova,state_lesha	

	c_d = c_d + 1.d0; status='drop'

	if(nmnm==1)return               ! nothing to drop yet :)

! select an unfortunate one
	j=rn(nmnm); vova=namelist(j)
	lesha=ass_r(vova); name=ass_l(vova)

! tau-s	
	tb = tau(lesha); t2=tau(vova); t1=tau(name)

! states
	state_vova=state(vova); state_lesha=state(lesha)

! ratio
	ratio = t*t*exp(-(state_vova-state_lesha)*t2)
	ratio = ratio * (t2+tb) * (tb+t1+t2)
	ratio = ratio * nmnm/(nmnm+2)
	ratio = 1.d0/ratio

! Metropolis
      if(ratio<1.d0)then; if(rndm()>ratio)return; endif

! accepted... update the configuration

! tau	
	tau(name) = tb + t1 + t2

! links	
	igor = ass_r(lesha)
	ass_r(name)=igor; ass_l(igor)=name

! make sure we've got masha still
	if(vova==masha)masha=name
	if(lesha==masha)masha=name

! drop
	call dropname(vova)
	call dropname(lesha)
	
	a_d = a_d + 1.d0

	end subroutine drop


!================ measure ===========================
	subroutine measure
!
! Measure energy only. Errorbars are estimated via the blocking method.
! 

	Z = Z + 1 
	ene = ene +  esti_ene()       ! energy
	i_b = i_b + 1.d0

!-------------- is a current block full?
	if( i_b == Z_b )then    ! wrap it up
	    b_n = b_n + 1; i_b=0
		ene_stat(b_n) = ene; ene = 0.d0

		if(b_n == b_n_max)then     ! collate blocks: 1000 -> 500 twice bigger blocks
	        call collate( ene_stat, b_n )
			b_n = b_n/2; Z_b = Z_b * 2.d0
		endif
	endif


	end subroutine measure


!--------------- energy estimator ----------------
	real*8 function esti_ene
	integer :: j,name
	real*8  :: sum

	name=masha; sum=0.d0 
	do j=1,nmnm
	   sum = sum + state(name)*tau(name)
	   name=ass_r(name)
	enddo
	
	esti_ene = ( sum - nmnm ) / beta

	end function esti_ene


!-------------- check configuration ----------------
	subroutine check
	integer :: j, name, vova

	name=masha; j=0
	do 
	   if(tau(name)<=0)then; print*,'check: tau ',step,name,tau(name)
				 stop
	   endif
	   vova=ass_r(name); 
           if(ass_l(vova)/=name)then; print*,'check: ',step,name, vova
			              stop
	   endif
	  
	   name=vova; j=j+1
	   if(name==masha)exit

	enddo
	
	if(j/=nmnm)then; print*,'check: windaround',j,nmnm
	stop
	endif

	end subroutine check


! ===================== I/O routines =============================


!----------- print out-------------------------
	subroutine prnt

	print*
	print*,'-----------------------'
	print*,'MC step (mln) = ', step/1d6
	print*,'nmnm now = ', nmnm
	print*

	if(.not.therm_flag)then   ! no statistics while thermalizing

	print*,'doing energy ( should be ',-sqrt(1.d0 + t**2),'):'	 
	call mrg(ene_stat(1:b_n),b_n,Z_b)

	endif

	print*
	print*,'add/drop: ',a_a/(c_a + 1d-6),a_d/(c_d+1d-6)


	end subroutine prnt


!-------- read cnf ---------------
	subroutine rd_cnf
	  print*,'rd_cnf stub'; stop
	end subroutine rd_cnf

!-------- read stat ---------------
	subroutine rd_stat
	  print*,'rd_stat stub'; stop
	end subroutine rd_stat

!-------- write to disk ---------
	subroutine wrt
	integer :: j
	real*8 :: summ
	  
! it's a stub. Writing configuration and statistics should be here.
	print*,'wrt stub...'


	summ = sum(nmnm_distr); if(summ==0.d0)summ=1.d0
	open(1,file='kinks.dat')
	do j=1,nmnm_max
	   write(1,*)j,nmnm_distr(j)/summ
	enddo
	close(1)

	end subroutine wrt

!-------- init stat ---------------------
	subroutine init_stat
	  
	Z_b = 1d3         ! something to start with
	i_b = 0; b_n = 0

	Z=0.d0; ene=0.d0; ene_stat=0.d0

	nmnm_distr=0.d0


	end subroutine init_stat

!-------- init cnf -------------------------
	subroutine init_cnf

! init NameManager data
	do j=1,nmnm_max
          namelist(j)=j; numnam(j)=j
      enddo; nmnm=0;
	
! Starting configuration
	call getname(masha)
	ass_l(masha)=masha; ass_r(masha)=masha
	tau(masha)=beta
	state(masha)=1


	end subroutine init_cnf


!==================== Name manager routines =========================


!----------------------
!---- Getname function for the Name Manager
!----------------------
      SUBROUTINE GetName(nick)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: nick

      IF(nmnm<nmnm_max)THEN;
          nmnm=nmnm+1
        nick=namelist(nmnm); numnam(nick)=nmnm
      ELSE
	PRINT*,'GetName-> list is over!'; stop
      ENDIF

      END SUBROUTINE GetName

!----------------------
!---- DropName function for Name Manager
!----------------------
      SUBROUTINE DropName(nick)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nick
      INTEGER :: nm, lastname

      nm=numnam(nick)
      IF(nm<nmnm)THEN
         lastname=namelist(nmnm); namelist(nm)=lastname
         numnam(lastname)=nm; namelist(nmnm)=nick
         numnam(nick)=nmnm; nmnm=nmnm-1
      ELSE IF(nm==nmnm)THEN
                 nmnm=nmnm-1
      ELSE
	PRINT*,'DropName-> No such name:',nick;  stop
      ENDIF
      RETURN

      END SUBROUTINE DropName


!-------- thermalization -------------
!
!  same MC loop, only without measuring
!
	subroutine therm
	integer :: j

	print*,'thermalization ...'

	therm_flag=.true.

	do  j=1,int(step_t)
	
          step=step+1.d0
          i_p=i_p+1.d0; i_w=i_w+1.d0; i_m=i_m+1.d0

!--- Update :)
           r = rndm()
           if     (r<prob(1))then;    call add
           else if(r<prob(2))then;    call drop
           else; print*,'am I nuts or what???'; stop
           endif

	     call check

           if(nmnm<nmnm_max) nmnm_distr(nmnm) = nmnm_distr(nmnm) + 1.d0

           if (i_p == step_p)  then; i_p=0.d0; call PRNT;    end if
           if (i_w == step_w)  then; i_w=0.d0; call WRT;     end if

	enddo

	therm_flag=.false.

	print*,'thermalization done via ',step_t/1d6, 'mln steps'

	end subroutine therm



!====================== Statistics management routines =======================


!----------- Collate the array ------------------
      subroutine collate(arr,n)
      real*8, dimension(*) :: arr
      integer :: n, Zb        ! array length & # of measurements per array entry

      integer :: i

      do i=1,n/2
          arr(i)=arr(2*i)+arr(2*i-1)
      enddo

      end subroutine collate

!-------------------------------
!--- Analyze block statistics
!-------------------------------
      subroutine bSTAT(arr,n,Zb,av,err)
      integer              :: n, Zb
      real*8, dimension(1:n) :: arr
      real*8               :: av, err

      real*8 :: av2, dif

      av  = sum( arr(1:n) )/Zb/n    
      av2 = sum( arr(1:n)**2 )/Zb/Zb/n

	dif=av2-av*av; if(dif<0)dif=0.d0
      err = sqrt( dif ) / sqrt(1.d0*n)                 


      end subroutine bSTAT


!-------------------------------
!--- Merge blocks & emit av +/- err
!-------------------------------
      subroutine mrg(arr,n,Zb)
      integer, intent(in)              :: n, Zb
      real*8, dimension(1:n), intent(in) :: arr

      real*8  :: av, err, arr1(1:n)
      integer :: i, n1, zb1


      zb1 = zb;       arr1(1:n) = arr(1:n); n1=n            

      print*,'-----------'

      do;

! emit
        call bSTAT(arr1,n1,zb1,av,err)
        print 777, av, err,n1               
 777    format(4x,g12.5,4x,' +/- ',g12.5,8x,I3)

! enough?
        if(n1<3)exit

! merge
        n1=INT(n1/2); zb1=zb1*2
        do i=1,n1
            arr1(i) =  arr1(2*i-1) + arr1(2*i)
        enddo

      enddo

      print*,'------------'; print*;


      end subroutine mrg


	end
