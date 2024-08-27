      module Globalparameters
        
		integer N		
        real*8 pi,aangle,minangle,maxangle,minxvalue,maxxvalue,a
	    real*8 minyvalue,maxyvalue,mintime,maxtime,alambda,tstar
        parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0)
        parameter(N=2,tstar = 300, a=0.015d0)
        parameter(
     +	pi = acos(-1.0),
     +	aangle    = pi/two,
     +	minangle  = aangle - 0.1*aangle,
     +	maxangle  = aangle + 0.1*aangle,
     +	minxvalue = 0.50d0,
     +	maxxvalue = 59.0d0,
     +	minyvalue = 0.0d0,
     +	maxyvalue = 10.0d0,
     +	mintime = 0.0d0,
     +	maxtime = 0.0d0,
     +	alambda = 0.2 ) 

      end module    
C***********************************************************************	  
      SUBROUTINE VEXTERNALDB(LOP,I_ARRAY,NIARRAY,R_ARRAY,NRARRAY)
	  use Globalparameters
      INCLUDE 'vaba_param.inc'
!-----------------------------------------------------------------------
!-----Include additional file for memory management
!-----------------------------------------------------------------------
#include <SMAAspUserSubroutines.hdr>
#include <mpif.h>
!-----------------------------------------------------------------------
!-----Declaration ABAQUS variables
!-----------------------------------------------------------------------
      dimension i_Array(niArray)
      dimension r_Array(nrArray)
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
!     Contents of i_Array
      parameter(i_int_nTotalNodes    = 1,
     *          i_int_nTotalElements = 2,
     *          i_int_kStep          = 3,
     *          i_int_kInc           = 4,
     *          i_int_iStatus        = 5,
     *          i_int_lWriteRestart  = 6)
!     Possible values for i_Array(i_int_iStatus)
      parameter(j_int_Continue          = 0,
     *          j_int_TerminateStep     = 1,
     *          j_int_TerminateAnalysis = 2)
!     Possible values for the lOp argument
      parameter(j_int_StartAnalysis  = 0,
     *          j_int_StartStep      = 1,
     *          j_int_SetupIncrement = 2,
     *          j_int_StartIncrement = 3,
     *          j_int_EndIncrement   = 4,
     *          j_int_EndStep        = 5,
     *          j_int_EndAnalysis    = 6 )
!     Contents of r_Array
      parameter(i_flt_TotalTime = 1,
     *          i_flt_StepTime  = 2,
     *          i_flt_dTime     = 3)
!-----------------------------------------------------------------------
!     Declaration of internal variables
!-----------------------------------------------------------------------
      integer kStep,kInc,kNel,KPROCESSNUM
	  integer myThreadID,check,i_error
      integer LENJOBNAME,LENOUTDIR,NUMPROCESSES
      character*256 JOBNAME,OUTDIR,filename,cwd
      character*1000 line
      REAL*8 x(N),y(N),angle(N),time(N),mindist        ! ARRAYs 
	  real*8 xx(N),yy(N),aaangle(N),xxx(N),aaaangle(N)
      integer ID_x,ID_y,ID_angle,ID_time,ID_inc              ! ID for pointers
      parameter(ID_x=1,ID_y=2,ID_angle=3,ID_time=4)   ! ID for pointers  	  
      pointer(ptr_x, x)         ! pointer link
	  pointer(ptr_y, y)         ! pointer link
	  pointer(ptr_angle, angle) ! pointer link
	  pointer(ptr_time, time)   ! pointer link
	  
	  CALL VGETJOBNAME( JOBNAME, LENJOBNAME )
	  ! cwd = 'P:\AKBAR\Axons Guidance\simulations\final tests\final\Different Growth Rates\New_with_coords\'
      ! filename = trim(cwd)//trim(JOBNAME)//'.txt'
!-----------------------------------------------------------------------
!     Initialization
!-----------------------------------------------------------------------
      kStep = i_Array(i_int_kStep)
      kInc  = i_Array(i_int_kInc)
      kNel  = i_Array(i_int_nTotalElements)
!--------------------------------------------------------------------------
!     Get threadID
!--------------------------------------------------------------------------
      call VGETRANK( KPROCESSNUM ) 
      CALL VGETNUMCPUS( NUMPROCESSES )	  
!-----------------------------------------------------------------------
!     Start of the analysis
!-----------------------------------------------------------------------
      if(lOp.eq.j_int_StartAnalysis)then 
!--------------------------------------------------------------------------
!          CREATE GLOBAL ARRAYS
!--------------------------------------------------------------------------
        ptr_x     = SMAFloatArrayCreate(ID_x, N, 0.0)
		ptr_y     = SMAFloatArrayCreate(ID_y, N, 0.0)
		ptr_angle = SMAFloatArrayCreate(ID_angle, N, 0.0)
		ptr_time  = SMAFloatArrayCreate(ID_time, N, 0.0)	
			
	    IF(kInc==0)THEN	
          if (KPROCESSNUM==0) then	
            !generate random initial points	
            call random_seed()
            call GetRandomBetweenRange(minxvalue, maxxvalue, x(1))
            call GetRandomBetweenRange(minyvalue, maxyvalue, y(1))
            call GetRandomBetweenRange(minangle, maxangle, angle(1))
            call GetRandomBetweenRange(mintime, maxtime, time(1))	
            if (N .gt. 1) then			
              do ii = 2, N
                call GetRandomBetweenRange(minyvalue, maxyvalue, y(ii))
                call GetRandomBetweenRange(minangle, maxangle, angle(ii))
                call GetRandomBetweenRange(mintime, maxtime, time(ii))
                check = 0   
                mindist	= 0.1
                do while (check .ne. 1)
                    call GetRandomBetweenRange(minxvalue, maxxvalue, x(ii))
                    check = 1
                    do jj = 1, N				    
                        if (((abs(x(ii)-x(jj))).lt.mindist) .and. (jj.ne.ii)) then
                            check = 0
                        endif
                    enddo
                enddo
              enddo	
            end if
          endif			  
          call MPI_BCAST(x, N, MPI_REAL8, 0, MPI_COMM_WORLD, i_error)
		  call MPI_BCAST(y, N, MPI_REAL8, 0, MPI_COMM_WORLD, i_error)
		  call MPI_BCAST(angle, N, MPI_REAL8, 0, MPI_COMM_WORLD, i_error)
		  call MPI_BCAST(time, N, MPI_REAL8, 0, MPI_COMM_WORLD, i_error)
	      call MPI_BARRIER( MPI_COMM_WORLD, i_error) 
        ENDIF

!-----------------------------------------------------------------------
!     Setup of the increment
!-----------------------------------------------------------------------
      elseif(lOp.eq.j_int_SetupIncrement)then
!-----------------------------------------------------------------------
!     Start of the increment
!-----------------------------------------------------------------------
      elseif(lOp.eq.j_int_StartIncrement)then  
	    ptr_x     = SMAFloatArrayAccess(ID_x)
	    ptr_y     = SMAFloatArrayAccess(ID_y)
	    ptr_angle = SMAFloatArrayAccess(ID_angle)
	    ptr_time  = SMAFloatArrayAccess(ID_time)
		call MPI_BARRIER( MPI_COMM_WORLD, i_error)
		call MPI_REDUCE(x, xx, N, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, i_error) 
		call MPI_REDUCE(x, xxx, N, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, i_error)
		call MPI_REDUCE(angle, aaangle, N, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, i_error) 
		call MPI_REDUCE(angle, aaaangle, N, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, i_error)		
		if (KPROCESSNUM==0) then
		    x = (xx + xxx)/2
			angle = (aaangle+aaaangle)/2
		end if
		call MPI_BARRIER( MPI_COMM_WORLD, i_error)
		call MPI_BCAST(x, N, MPI_REAL8, 0, MPI_COMM_WORLD, i_error)
		call MPI_BCAST(angle, N, MPI_REAL8, 0, MPI_COMM_WORLD, i_error)
		call MPI_ALLREDUCE(y, yy, N, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, i_error)
		y = yy
		
        call MPI_BARRIER( MPI_COMM_WORLD, i_error)		
		! if (KPROCESSNUM==0) then
		    ! open(unit=4, file=filename, STATUS = 'unknown', ACTION = 'write', POSITION='append')
			! WRITE(4, '(100F10.4, 100F10.4, 100F10.4)') (x(ii), ii=1,N), (y(ii), ii=1,N), (angle(ii), ii=1,N)
		    ! close(4)
		! end if		
		
		
!-----------------------------------------------------------------------
!     End of the increment
!-----------------------------------------------------------------------
      elseif(lOp.eq.j_int_EndIncrement)then	    
!-----------------------------------------------------------------------
!     End of the analysis
!-----------------------------------------------------------------------
      elseif(lOp.eq.j_int_EndAnalysis)then	
      endif
!-----------------------------------------------------------------------
!     End of the subroutine
!-----------------------------------------------------------------------
      return
      end
C***********************************************************************	  
c
c User subroutine VUCHARLENGTH for user-defined element characteristic length
c
      subroutine vucharlength(
c Read only -
     *     nblock, nfieldv, nprops, ncomp, ndim, nnode, nstatev,
     *     kSecPt, kLayer, kIntPt, jElType, jElem,
     *     totalTime, stepTime, dt,
     *     cmname, coordMp, coordNode, direct, T, props, 
     *     field, stateOld,
c Write only -
     *     charLength )
*
      include 'vaba_param.inc'
*     
      dimension jElType(3), jElem(nblock), coordMp(nblock, ndim), 
     *     coordNode(nblock,nnode,ndim),
     *     direct(nblock, 3, 3), T(nblock,3,3), props(nprops),
     *     stateOld(nblock, nstatev), charLength(nblock,ncomp), 
     *     field(nblock, nfieldv)
      character*80 cmname
	  real*8 x1,x2,x3,x4,y1,y2,y3,y4
*      
*T2D2, T3D2, SAX1
      if( jElType(1) .eq. 1 ) then
         do k = 1, nblock
            charLength(k,1) = props(1)
         end do
      end if
* 
*S4R, S4RS, S4RSW, S4, CPS4R, CPE4R, CAX4R, M3D4R, M3D4
      if( jElType(1) .eq. 3 ) then
         do k = 1, nblock
            x1 = coordNode(k,1,1)
			y1 = coordNode(k,1,2)
			x2 = coordNode(k,2,1)
			y2 = coordNode(k,2,2)
			x3 = coordNode(k,3,1)
			y3 = coordNode(k,3,2)
			x4 = coordNode(k,4,1)
			y4 = coordNode(k,4,2)	 
            charLength(k,1)=x1
			charLength(k,2)=x2
			charLength(k,3)=x3
			charLength(k,4)=x4
			charLength(k,5)=y1
			charLength(k,6)=y2
			charLength(k,7)=y3
			charLength(k,8)=y4
         end do
      end if   
*
*C3D8R, C3D8, SC8R
      if( jElType(1) .eq. 6 ) then
         do k = 1, nblock
            diagonal_1 = sqrt( (coordNode(k,1,1)-coordNode(k,7,1) )**2 +
     *           (coordNode(k,1,2)-coordNode(k,7,2))**2 +
     *           (coordNode(k,1,3)-coordNode(k,7,3))**2 )
            diagonal_2 = sqrt( (coordNode(k,2,1)-coordNode(k,8,1) )**2 +
     *           (coordNode(k,2,2)-coordNode(k,8,2))**2 +
     *           (coordNode(k,2,3)-coordNode(k,8,3))**2 )
            diagonal_3 = sqrt( (coordNode(k,3,1)-coordNode(k,5,1) )**2 +
     *           (coordNode(k,3,2)-coordNode(k,5,2))**2 +
     *           (coordNode(k,3,3)-coordNode(k,5,3))**2 )
            diagonal_4 = sqrt( (coordNode(k,4,1)-coordNode(k,6,1) )**2 +
     *           (coordNode(k,4,2)-coordNode(k,6,2))**2 +
     *           (coordNode(k,4,3)-coordNode(k,6,3))**2 )
            charLength(k,1)=props(1)
         end do
      end if   
*
      return
      end
C***********************************************************************
      subroutine vumat (
c Read only -
     +     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     +     stressNew, stateNew, enerInternNew, enerInelasNew )
c
	  use Globalparameters      
	  include 'vaba_param.inc'
#include <SMAAspUserSubroutines.hdr>
#include <mpif.h>
c
      dimension jblock(*), props(nprops),density(*), coordMp(*),
     +     charLength(*), strainInc(*),
     +     relSpinInc(*), tempOld(*),
     +     stretchOld(*),
     +     defgradOld(*),
     +     fieldOld(*), stressOld(*),
     +     stateOld(*), enerInternOld(*),
     +     enerInelasOld(*), tempNew(*),
     +     stretchNew(*),
     +     defgradNew(*),
     +     fieldNew(*),
     +     stressNew(*), stateNew(*),
     +     enerInternNew(*), enerInelasNew(*)
c
      character*80 cmname
      character*256 WHITE,GRAY,SKULL

      parameter (     
     +     i_umt_nblock = 1,
     +     i_umt_npt    = 2,
     +     i_umt_layer  = 3,
     +     i_umt_kspt   = 4,
     +     i_umt_noel   = 5 )
	  
      IF(cmname(1:7) .EQ. 'GRAY') THEN
      !
      ! gray matter 
      !
        call  vumatXtrArg_gray (jblock(i_umt_nblock),
     +     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
     +     jblock(i_umt_noel), jblock(i_umt_npt),
     +     jblock(i_umt_layer), jblock(i_umt_kspt))
	 
	 else if (cmname(1:7) .EQ. 'WHITE') then
      !
      ! white matter 
      !	 
        call  vumatXtrArg_white (jblock(i_umt_nblock),
     +     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
     +     jblock(i_umt_noel), jblock(i_umt_npt),
     +     jblock(i_umt_layer), jblock(i_umt_kspt))	 	 

      end if
      end subroutine vumat

C***********************************************************************

C User subroutine VUMAT
C
      subroutine vumatXtrArg_white (
C Read only -
     *     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     *     stressNew, stateNew, enerInternNew, enerInelasNew )
C
      use Globalparameters
      include 'vaba_param.inc'
#include <SMAAspUserSubroutines.hdr>
#include <mpif.h>
C
      dimension coordMp(nblock,*), charLength(nblock,8), props(nprops),
     1     density(nblock), strainInc(nblock,ndir+nshr),
     2     relSpinInc(nblock,nshr), tempOld(nblock),
     3     stretchOld(nblock,ndir+nshr),
     4     defgradOld(nblock,ndir+nshr+nshr),
     5     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6     stateOld(nblock,nstatev), enerInternOld(nblock),
     7     enerInelasOld(nblock), tempNew(nblock),
     8     stretchNew(nblock,ndir+nshr),
     9     defgradNew(nblock,ndir+nshr+nshr),eigVec(3,3),
     1     fieldNew(nblock,nfieldv), eigVal(nblock*3), 
     2     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3     enerInternNew(nblock), enerInelasNew(nblock)
C     
      character*256 cmname, filename, OUTDIR, JOBNAME
	  integer LENJOBNAME,LENOUTDIR,KPROCESSNUM,myThreadID,r

      integer i,j,ii,jj,kk,km,damage,process_rank,kblock,check
	  integer ABA_COMM_WORLD,counting,i_error
      real*8 Iden(3,3),F_tau(3,3),U_tau(3,3),Be_tau(3,3),a
      real*8 T_tau(3,3),R_tau(3,3),U_inv(3,3),detF,trace,Gsub
      real*8 pwrinct,stress_power,Je,JU,detU,Jg,b0(3,3),a0(3,3)
      real*8 c10,d1,xI1,range,n1(2),nA(3),stress(ndir+nshr)
	  real*8 coeff1,trbmat,xkirch(3,3),coeff3,xener,randi,tempRand
	  real*8 x(N),y(N),angle(N),time(N),minx,maxx,miny,maxy,noise
	  real*8 gradient,xm,ym,gradCx,gradCy,dangle,c10g,thetag_t,thetag_tau
	  real*8 twomu,alamda,xnu,K0,prindir,tetha,rotation_matrix(2,2)
	  real*8 x11,x22,x33,x44,y11,y22,y33,y44,R_A(3),growth_rate
	  real*8 Trotation_matrix(2,2),principlestress(2,2),stressTensor(2,2)
	  real*8 sigmax,sigmin,deltaUx,deltaUy,xx(N),yy(N),aaangle(N)
	  real*8 omega(3),norm_omega,normal_omega(3),cross_omegaTonA(3)
	  real*8 Fg_tau(3,3),Fe_tau(3,3),Fginv(3,3)
      !-----Declare global variables
      integer ID_x,ID_y,ID_angle,ID_time          ! ID for pointers
      parameter(ID_x=1,ID_y=2,ID_angle=3,ID_time=4)   ! ID for pointers     
      pointer(ptr_x, x)         ! pointer link
	  pointer(ptr_y, y)         ! pointer link
	  pointer(ptr_angle, angle) ! pointer link
	  pointer(ptr_time, time)   ! pointer link	   
	  
      ! Obtain material properties 
      !
      c10         = props(1)
      d1          = props(2)
	  Gsub        = props(3)
  
      twomu  = four*c10
	  K0 = 2/d1
	  xnu = (6*K0/twomu - 2)/(12*K0/twomu + 2)
      alamda = twomu * xnu / ( one - two * xnu )	  

      do kblock = 1, nblock
         do k = 1, ndir+nshr
            stressNew(kblock,k) = zero
         end do
         enerInternNew(kblock) = zero
         enerInelasNew(kblock) = zero
      end do
	  
      !If stepTime equals to zero, assume the material pure elastic 
      !and use initial elastic modulus
      if ( stepTime .eq. zero ) then  
        do k = 1, nblock
          !Trial stress
          trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
          stressNew(k,1) = stressOld(k,1) 
     *         + twomu * strainInc(k,1) + alamda * trace
          stressNew(k,2) = stressOld(k,2) 
     *         + twomu * strainInc(k,2) + alamda * trace
          stressNew(k,3) = stressOld(k,3) 
     *         + twomu * strainInc(k,3) + alamda * trace
          stressNew(k,4)=stressOld(k,4) + twomu * strainInc(k,4)
		  stateOld(k,1) = zero
		  stateOld(k,2) = c10
          stateOld(k,3) = coordMp(k,1)
          stateOld(k,4) = coordMp(k,2)
          stateOld(k,5) = one		  
        end do		
      else
	 				
       !
       ! START LOOP OVER MATERIAL POINTS:	   
       do km=1,nblock
   
        c10   = props(1)

        !initialize state variables
        stateNew(km,1) = stateOld(km,1)   !!damage
		stateNew(km,2) = stateOld(km,2)   !!gradient
		
        !reading x , y coordinates of axon's previous position		
		damage       =    stateOld(km,1)
		c10g         =    stateOld(km,2)

		!obtain middle point coordinates of the current material point
        xm = coordMp(km,1)
        ym = coordMp(km,2)
        stateNew(km,3) = xm
        stateNew(km,4) = ym		
		
        !obtain nodal coordinates of the elements
        x11 = charLength(km,1)
		x22 = charLength(km,2)
		x33 = charLength(km,3)
		x44 = charLength(km,4)
		y11 = charLength(km,5)
		y22 = charLength(km,6)
		y33 = charLength(km,7)
		y44 = charLength(km,8)
		
		!find the minimum and maximum values of the element
		minx = min(x11,x22,x33,x44)
		maxx = max(x11,x22,x33,x44)
		miny = min(y11,y22,y33,y44)
		maxy = max(y11,y22,y33,y44)	  
		
        !check if the initial point is inside the current element. If so, active it.
		!
	    ptr_x     = SMAFloatArrayAccess(ID_x)
	    ptr_y     = SMAFloatArrayAccess(ID_y)
	    ptr_angle = SMAFloatArrayAccess(ID_angle)
	    ptr_time  = SMAFloatArrayAccess(ID_time)

		do jj = 1, N
            if (totalTime.ge.time(jj)) then	
		        ! generate random numbers
                call random_seed()
	            call random_number(randi)
	            noise = two*randi-one
	            call point_in_quad(x(jj), y(jj), r, x11, y11, x22, y22, x33, y33, x44, y44)
				if (r .eq. 1) then	
		            damage = 1
			        c10g = c10*2

		            !!!calculating priciple values and directions
                    do 10 i=1,ndir
                        stress(i) = stressOld(km,i)
   10               continue
                    do 20 i = ndir + 1, nshr + ndir
                        stress(i) = stressOld(km,i)            
   20               continue
   
                    call sprindir(stress,principlestress,tetha)
                    ! find maximum principal stress and corresponding direction
	                index = 1
                    max_prin_SIG = principlestress(index,index)
                    if (principlestress(2,2) .gt. max_prin_SIG) then
			            max_prin_SIG = principlestress(2,2)
				        index = 2
			        endif					

					if (max_prin_SIG .gt. zero) then
					    prindir = max_prin_SIG
					    if (index .eq. 1) then
						    n1(1) = cos(tetha)
						    n1(2) = sin(tetha)
							n1(3) = zero 
							theta_n1nA = acos(dot_product(n1,nA)/(norm2(n1)*norm2(nA)))
							if (theta_n1nA .gt. pi/2) then
							    n1 = -n1
							endif
						else
						    tetha = tetha + pi/2
					        n1(1) = cos(tetha)
						    n1(2) = sin(tetha)
							n1(3) = zero  	
							theta_n1nA = acos(dot_product(n1,nA)/(norm2(n1)*norm2(nA)))
							if (theta_n1nA .gt. pi/2) then
							    n1 = -n1
							endif						
						endif
				    else 
					    tetha = zero
					endif
					
                    dangle = 0.125*noise
			        angle(jj) = angle(jj) + 0.1*dangle
				
			        if ((tetha.ne.zero).and.(tetha.ne.angle(jj)).and.(abs(tetha-angle(jj)).gt.0.1)) then
			            nA(1) = cos(angle(jj))
			            nA(2) = sin(angle(jj))
						nA(3) = 0						
						omega(1) = zero
						omega(2) = zero
						omega(3) = (pi/(2*tstar))*(nA(1)*n1(2)-nA(2)*n1(1))
						norm_omega = abs(omega(3))
						normal_omega = omega/abs(omega(3))
						cross_omegaTonA(1)= -normal_omega(3)*nA(2)
						cross_omegaTonA(2)= normal_omega(3)*nA(1)
						cross_omegaTonA(3) = zero
						R_A = cos(0.1d0*norm_omega)*nA+sin(0.1d0*norm_omega)*cross_omegaTonA
						R_A = R_A + (1-cos(0.1d0*norm_omega))*(dot_product(normal_omega,nA))*normal_omega		
		                angle(jj) = atan2(R_A(2),R_A(1))
			        end if
                    
					deltaUx = stateNew(km,3)-stateOld(km,3)
					deltaUy = stateNew(km,4)-stateOld(km,4)
					
					growth_rate = 1.0 + abs(prindir)*a*1000000.0d0
			        x(jj) = x(jj) + growth_rate*dt*cos(angle(jj)) + deltaUx
			        y(jj) = y(jj) + growth_rate*dt*sin(angle(jj)) + deltaUy
		
		        end if
			end if
		end do
		
        !!update state variables
        stateNew(km,1) = damage
		stateNew(km,2) = c10g
		c10            = c10g
		
        ! Copy old and new deformation gradients
        !
        F_tau(1,1) = defgradNew(km,1)
        F_tau(2,2) = defgradNew(km,2)
        F_tau(3,3) = defgradNew(km,3)
        F_tau(1,2) = defgradNew(km,4)
        U_tau(1,1) = stretchNew(km,1)
        U_tau(2,2) = stretchNew(km,2)
        U_tau(3,3) = stretchNew(km,3)
        U_tau(1,2) = stretchNew(km,4)
        if(nshr .lt. 2) then
          ! 2D case
          F_tau(2,1) = defgradNew(km,5)
          F_tau(1,3) = zero
          F_tau(2,3) = zero
          F_tau(3,1) = zero
          F_tau(3,2) = zero
          U_tau(2,1) = U_tau(1,2)
          U_tau(1,3) = zero
          U_tau(2,3) = zero
          U_tau(3,1) = zero
          U_tau(3,2) = zero
        else
          ! 3D case
          F_tau(2,3) = defgradNew(km,5)
          F_tau(3,1) = defgradNew(km,6)
          F_tau(2,1) = defgradNew(km,7)
          F_tau(3,2) = defgradNew(km,8)
          F_tau(1,3) = defgradNew(km,9)
          U_tau(2,3) = stretchNew(km,5)
          U_tau(3,1) = stretchNew(km,6)
          U_tau(2,1) = U_tau(1,2)
          U_tau(3,2) = U_tau(2,3)
          U_tau(1,3) = U_tau(3,1)
        end if

        ! Identity matrix for later use.
        !
        call onem(Iden)
		
        ! Compute the relative volume change
        !
        call mdet(F_tau,detF)
		
        ! Compute the determinant of U_tau
        !
        call mdet(U_tau,detU)

        ! Compute the inverse of the U_tau
        !
	    U_inv(1,1) = (U_tau(2,2)*U_tau(3,3) - U_tau(2,3)*U_tau(3,2))/detU
	    U_inv(1,2) = (U_tau(1,3)*U_tau(3,2) - U_tau(1,2)*U_tau(3,3))/detU
	    U_inv(1,3) = (U_tau(1,2)*U_tau(2,3) - U_tau(1,3)*U_tau(2,2))/detU
	    U_inv(2,1) = (U_tau(2,3)*U_tau(3,1) - U_tau(2,1)*U_tau(3,3))/detU
	    U_inv(2,2) = (U_tau(1,1)*U_tau(3,3) - U_tau(1,3)*U_tau(3,1))/detU
	    U_inv(2,3) = (U_tau(1,3)*U_tau(2,1) - U_tau(1,1)*U_tau(2,3))/detU
	    U_inv(3,1) = (U_tau(2,1)*U_tau(3,2) - U_tau(2,2)*U_tau(3,1))/detU
	    U_inv(3,2) = (U_tau(1,2)*U_tau(3,1) - U_tau(1,1)*U_tau(3,2))/detU
	    U_inv(3,3) = (U_tau(1,1)*U_tau(2,2) - U_tau(1,2)*U_tau(2,1))/detU		

        ! compute the R=FU^-1
        R_tau = matmul(F_tau,U_inv)
		
		! Read old state variables         
        thetag_t = stateOld(km,5) ! growth parameter at time t
        thetag_tau = thetag_t  
        thetag_tau = thetag_tau + (Gsub)*dt 

        if(stepTime.eq.zero) then		
		  Fg_tau = Iden
		else
          Fg_tau = ((thetag_tau)**third)*Iden 
		end if
		
        ! Jacobian of the Fg
        ! 
	    Jg = Fg_tau(1,1) * Fg_tau(2,2) * Fg_tau(3,3) -
     *	    Fg_tau(1,2) * Fg_tau(2,1) * Fg_tau(3,3) +
     *	    Fg_tau(1,2) * Fg_tau(2,3) * Fg_tau(3,1) +
     *	    Fg_tau(1,3) * Fg_tau(3,2) * Fg_tau(2,1) -
     *	    Fg_tau(1,3) * Fg_tau(3,1) * Fg_tau(2,2) -
     *	    Fg_tau(2,3) * Fg_tau(3,2) * Fg_tau(1,1)		
		
        ! Jacobian of Fe
        Je = detF/Jg
		
        ! inverse of the growth Fg
        !
	    Fginv(1,1) = (Fg_tau(2,2)*Fg_tau(3,3) - Fg_tau(2,3)*Fg_tau(3,2))/Jg
	    Fginv(1,2) = (Fg_tau(1,3)*Fg_tau(3,2) - Fg_tau(1,2)*Fg_tau(3,3))/Jg
	    Fginv(1,3) = (Fg_tau(1,2)*Fg_tau(2,3) - Fg_tau(1,3)*Fg_tau(2,2))/Jg
	    Fginv(2,1) = (Fg_tau(2,3)*Fg_tau(3,1) - Fg_tau(2,1)*Fg_tau(3,3))/Jg
	    Fginv(2,2) = (Fg_tau(1,1)*Fg_tau(3,3) - Fg_tau(1,3)*Fg_tau(3,1))/Jg
	    Fginv(2,3) = (Fg_tau(1,3)*Fg_tau(2,1) - Fg_tau(1,1)*Fg_tau(2,3))/Jg
	    Fginv(3,1) = (Fg_tau(2,1)*Fg_tau(3,2) - Fg_tau(2,2)*Fg_tau(3,1))/Jg
	    Fginv(3,2) = (Fg_tau(1,2)*Fg_tau(3,1) - Fg_tau(1,1)*Fg_tau(3,2))/Jg
	    Fginv(3,3) = (Fg_tau(1,1)*Fg_tau(2,2) - Fg_tau(1,2)*Fg_tau(2,1))/Jg		
		
        ! elastic Fe
        ! 
        Fe_tau = matmul(F_tau,Fginv)
		
        ! Left Cauchy Green tensor  
        ! 
        Be_tau = matmul(Fe_tau,transpose(Fe_tau))
		
		xI1 = Be_tau(1,1)+Be_tau(2,2)+Be_tau(3,3)
		
		coeff1 = two*c10/(Je**(two/three))
		
		trbmat= coeff1*xI1/three
		
C     KIRCHHOFF STRESS PART 1   
        do I=1,3
          do J=1,3		
            T_tau(I,J) = Be_tau(I,J) * coeff1			
          end do
        end do
	  
C     SUBTRACT THE PART 2   
        do I = 1,3	  
          T_tau(I,I) = T_tau(I,I) - trbmat		
        end do
	  
C     FORM PART 3      
        if (d1 .eq. zero) then
           coeff3 = zero
        else		   
           coeff3 = two*(Je-one)*Je/d1
		end if
	  
C     ADD TO THE PREVIOUS PARTS      
        do I = 1,3  
          T_tau(I,I) = T_tau(I,I) + coeff3	
        end do	
	
        ! ABAQUS/Explicit uses stress measure (transpose(R) T R)
        !
		T_tau = T_tau/Je
        T_tau = matmul(transpose(R_tau),matmul(T_tau,R_tau))

        do i=1,ndir
          stressNew(km,i) = T_tau(i,i)
        end do
        if(nshr.ne.0) then
          stressNew(km,ndir+1) = T_tau(1,2)
          if(nshr.ne.1) then
            stressNew(km, ndir+2) = T_tau(2,3)
            if(nshr.ne.2) then
              stressNew(km,ndir+3) = T_tau(1,3)
            end if
          end if
        end if
		
		stateNew(km,5) = thetag_tau ! growth parameter at time tau

        ! Update the specific internal energy
        !
        stress_power = 0.d0
        do i = 1,ndir
          stress_power = stress_power +
     +           0.5*((StressOld(km,i)+StressNew(km,i))*
     +           StrainInc(km,i))
        end do
         
        select case (nshr)
         case(1)
            stress_power = stress_power + 
     +           0.5*((StressOld(km,ndir+1)+StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1))
         case(3)
            stress_power = stress_power + 
     +           0.5*(((StressOld(km,ndir+1) + StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1)) +
     +           ((StressOld(km,ndir+2)+ StressNew(km,ndir+2)) *
     +           StrainInc(km,ndir+2))+
     +           ((StressOld(km,ndir+3) + StressNew(km,ndir+3))*
     +           StrainInc(km,ndir+3)))
        end select

		if (d1 .eq. zero) then
		   xener=c10*(xI1/(Je**(two/three))-three)
		else
           xener=c10*(xI1/(Je**(two/three))-three)+(one/d1)*((Je-one)**two)
		end if
        enerInternNew(km) = xener/density(km)
	  
C        enerInternNew(km) = enerInternOld(km) + 
C     +        stress_power/density(km)
           
        enerInelasNew(km) = enerInelasOld(km) + 
     +        pwrinct/density(km)
           
       end do
	  end if
	  return
	  end
 	  
C****************************************************************************
C***********************************************************************
      subroutine vumatXtrArg_gray (
c Read only -
     +     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
c Read only extra arguments -
     +     nElement, nMatPoint, nLayer, nSecPoint )

      use Globalparameters
      include 'vaba_param.inc'

      dimension props(nprops), density(nblock), coordMp(nblock,*),
     +     charLength(nblock), strainInc(nblock,ndir+nshr),
     +     relSpinInc(nblock,nshr), tempOld(nblock),
     +     stretchOld(nblock,ndir+nshr),
     +     defgradOld(nblock,ndir+nshr+nshr),
     +     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     +     stateOld(nblock,nstatev), enerInternOld(nblock),
     +     enerInelasOld(nblock), tempNew(nblock),
     +     stretchNew(nblock,ndir+nshr),
     +     defgradNew(nblock,ndir+nshr+nshr),
     +     fieldNew(nblock,nfieldv),
     +     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     +     enerInternNew(nblock), enerInelasNew(nblock)

      dimension nElement(nblock)

      character*80 cmname

      integer i,j,l,i1,j1,ii,jj,kk,ll,km,ifail

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),U_t(3,3),U_tau(3,3),Fp_t(3,3)
      real*8 Fp_tau(3,3),Me_t(3,3),Me_tau(3,3),nuP_t,nuP_tau,Y_t,Y_tau
      real*8 gBarP_t,gBarP_tau,T_tau(3,3),R_tau(3,3),U_inv(3,3),detF
      real*8 Fp_inv(3,3),Ee_tau(3,3),Re_tau(3,3),Ue_tau(3,3)
      real*8 pnu0,damage_t,damage_tau,mag_Dp_tau,pwrinct,stress_power
      real*8 nu1,nu3,nu5, Je,JU,detU,Jg,xini,yini,zini,rini
      real*8 lamg_t,lamg_tau,matProps(nprops)
      real*8 E_tau,E_t,thetag_t,thetag_tau,mu_g,lambda_g,Gctx
      real*8 coordx,coordy,coordz,b0(3,3),a0(3,3)
      real*8 curvature,NormCurve,G_total,Gcur,c10,d1,xI1
	  real*8 Fg_tau(3,3),Fginv(3,3),Be_tau(3,3),Fe_tau(3,3)
	  real*8 coeff1,trbmat,xkirch(3,3),coeff3,xener
	  real*8 twomu,K0,xnu,alamda
	  
      ! Obtain material properties 
      !
      c10         = props(1)
      d1          = props(2)
      Gctx        = props(3)

      twomu  = four*c10
	  K0 = 2/d1
	  xnu = (6*K0/twomu - 2)/(12*K0/twomu + 2)
      alamda = twomu * xnu / ( one - two * xnu )	  

      do kblock = 1, nblock
         do k = 1, ndir+nshr
            stressNew(kblock,k) = zero
         end do
         enerInternNew(kblock) = zero
         enerInelasNew(kblock) = zero
      end do

      !If stepTime equals to zero, assume the material pure elastic 
      !and use initial elastic modulus
      if ( stepTime .eq. zero ) then  
        do k = 1, nblock
          !Trial stress
          trace = strainInc(k,1) + strainInc(k,2) + strainInc(k,3)
          stressNew(k,1) = stressOld(k,1) 
     *         + twomu * strainInc(k,1) + alamda * trace
          stressNew(k,2) = stressOld(k,2) 
     *         + twomu * strainInc(k,2) + alamda * trace
          stressNew(k,3) = stressOld(k,3) 
     *         + twomu * strainInc(k,3) + alamda * trace
          stressNew(k,4)=stressOld(k,4) + twomu * strainInc(k,4)
		  stateOld(k,1) = one	  
        end do		
      else

       do km=1,nblock
 
        ! Copy old and new deformation gradients
        !
        F_tau(1,1) = defgradNew(km,1)
        F_tau(2,2) = defgradNew(km,2)
        F_tau(3,3) = defgradNew(km,3)
        F_tau(1,2) = defgradNew(km,4)
        U_tau(1,1) = stretchNew(km,1)
        U_tau(2,2) = stretchNew(km,2)
        U_tau(3,3) = stretchNew(km,3)
        U_tau(1,2) = stretchNew(km,4)
        if(nshr .lt. 2) then
          ! 2D case
          F_tau(2,1) = defgradNew(km,5)
          F_tau(1,3) = zero
          F_tau(2,3) = zero
          F_tau(3,1) = zero
          F_tau(3,2) = zero
          U_tau(2,1) = U_tau(1,2)
          U_tau(1,3) = zero
          U_tau(2,3) = zero
          U_tau(3,1) = zero
          U_tau(3,2) = zero
        else
          ! 3D case
          F_tau(2,3) = defgradNew(km,5)
          F_tau(3,1) = defgradNew(km,6)
          F_tau(2,1) = defgradNew(km,7)
          F_tau(3,2) = defgradNew(km,8)
          F_tau(1,3) = defgradNew(km,9)
          U_tau(2,3) = stretchNew(km,5)
          U_tau(3,1) = stretchNew(km,6)
          U_tau(2,1) = U_tau(1,2)
          U_tau(3,2) = U_tau(2,3)
          U_tau(1,3) = U_tau(3,1)
        end if

        ! Identity matrix for later use.
        !
        call onem(Iden)
		
        ! Compute the relative volume change
        !
        call mdet(F_tau,detF)
		
        ! Compute the determinant of U_tau
        !
        call mdet(U_tau,detU)

        ! Compute the inverse of the U_tau
        !
	    U_inv(1,1) = (U_tau(2,2)*U_tau(3,3) - U_tau(2,3)*U_tau(3,2))/detU
	    U_inv(1,2) = (U_tau(1,3)*U_tau(3,2) - U_tau(1,2)*U_tau(3,3))/detU
	    U_inv(1,3) = (U_tau(1,2)*U_tau(2,3) - U_tau(1,3)*U_tau(2,2))/detU
	    U_inv(2,1) = (U_tau(2,3)*U_tau(3,1) - U_tau(2,1)*U_tau(3,3))/detU
	    U_inv(2,2) = (U_tau(1,1)*U_tau(3,3) - U_tau(1,3)*U_tau(3,1))/detU
	    U_inv(2,3) = (U_tau(1,3)*U_tau(2,1) - U_tau(1,1)*U_tau(2,3))/detU
	    U_inv(3,1) = (U_tau(2,1)*U_tau(3,2) - U_tau(2,2)*U_tau(3,1))/detU
	    U_inv(3,2) = (U_tau(1,2)*U_tau(3,1) - U_tau(1,1)*U_tau(3,2))/detU
	    U_inv(3,3) = (U_tau(1,1)*U_tau(2,2) - U_tau(1,2)*U_tau(2,1))/detU		

        ! compute the R=FU^-1
        R_tau = matmul(F_tau,U_inv)
		
        ! Read old state variables         
        thetag_t = stateOld(km,1) ! growth parameter at time t
        thetag_tau = thetag_t  
        thetag_tau = thetag_tau + (Gctx)*dt   
		
        ! update  kinematics 
        ! area growth 	
        a0(1,1) = 0
        a0(2,1) = 0
        a0(3,1) = 1

        b0(1,1) = a0(1,1)*a0(1,1)
        b0(1,2) = a0(1,1)*a0(2,1)
		b0(1,3) = a0(1,1)*a0(3,1)
		b0(2,1) = a0(2,1)*a0(1,1)
		b0(2,2) = a0(2,1)*a0(2,1)
		b0(2,3) = a0(2,1)*a0(3,1)
		b0(3,1) = a0(3,1)*a0(1,1)
		b0(3,2) = a0(3,1)*a0(2,1)
		b0(3,3) = a0(3,1)*a0(3,1)	        
		
        if(stepTime.eq.zero) then		
		  Fg_tau = Iden
		else
          Fg_tau = dsqrt(thetag_tau)*Iden 
     +          +(one - dsqrt(thetag_tau))*b0
		end if
 
        ! Jacobian of the Fg
        ! 
	    Jg = Fg_tau(1,1) * Fg_tau(2,2) * Fg_tau(3,3) -
     *	    Fg_tau(1,2) * Fg_tau(2,1) * Fg_tau(3,3) +
     *	    Fg_tau(1,2) * Fg_tau(2,3) * Fg_tau(3,1) +
     *	    Fg_tau(1,3) * Fg_tau(3,2) * Fg_tau(2,1) -
     *	    Fg_tau(1,3) * Fg_tau(3,1) * Fg_tau(2,2) -
     *	    Fg_tau(2,3) * Fg_tau(3,2) * Fg_tau(1,1)		
		
        ! Jacobian of Fe
        Je = detF/Jg
		
        ! inverse of the growth Fg
        !
	    Fginv(1,1) = (Fg_tau(2,2)*Fg_tau(3,3) - Fg_tau(2,3)*Fg_tau(3,2))/Jg
	    Fginv(1,2) = (Fg_tau(1,3)*Fg_tau(3,2) - Fg_tau(1,2)*Fg_tau(3,3))/Jg
	    Fginv(1,3) = (Fg_tau(1,2)*Fg_tau(2,3) - Fg_tau(1,3)*Fg_tau(2,2))/Jg
	    Fginv(2,1) = (Fg_tau(2,3)*Fg_tau(3,1) - Fg_tau(2,1)*Fg_tau(3,3))/Jg
	    Fginv(2,2) = (Fg_tau(1,1)*Fg_tau(3,3) - Fg_tau(1,3)*Fg_tau(3,1))/Jg
	    Fginv(2,3) = (Fg_tau(1,3)*Fg_tau(2,1) - Fg_tau(1,1)*Fg_tau(2,3))/Jg
	    Fginv(3,1) = (Fg_tau(2,1)*Fg_tau(3,2) - Fg_tau(2,2)*Fg_tau(3,1))/Jg
	    Fginv(3,2) = (Fg_tau(1,2)*Fg_tau(3,1) - Fg_tau(1,1)*Fg_tau(3,2))/Jg
	    Fginv(3,3) = (Fg_tau(1,1)*Fg_tau(2,2) - Fg_tau(1,2)*Fg_tau(2,1))/Jg		
		
        ! elastic Fe
        ! 
        Fe_tau = matmul(F_tau,Fginv)
		
        ! Left Cauchy Green tensor  
        ! 
        Be_tau = matmul(Fe_tau,transpose(Fe_tau))
		
		xI1 = Be_tau(1,1)+Be_tau(2,2)+Be_tau(3,3)
		
		coeff1 = two*c10/(Je**(two/three))
		
		trbmat= coeff1*xI1/three
		
C     KIRCHHOFF STRESS PART 1   
        DO I=1,3
          DO J=1,3		
            T_tau(I,J) = Be_tau(I,J) * coeff1			
          END DO
        END DO
	  
C     SUBTRACT THE PART 2   
        DO I = 1,3	  
          T_tau(I,I) = T_tau(I,I) - trbmat		
        END DO
	  
C     FORM PART 3      
        if (d1 .eq. zero) then
           coeff3 = zero
        else		   
           coeff3 = two*(Je-one)*Je/d1
		end if
	  
C     ADD TO THE PREVIOUS PARTS      
        DO I = 1,3  
          T_tau(I,I) = T_tau(I,I) + coeff3	
        END DO	 

        ! Update state variables
        !
        stateNew(km,1) = thetag_tau ! growth parameter at time tau
        stateNew(km,2) = curvature

        ! ABAQUS/Explicit uses stress measure (transpose(R) T R)
        !
		T_tau = T_tau/Je
        T_tau = matmul(transpose(R_tau),matmul(T_tau,R_tau))

        do i=1,ndir
          stressNew(km,i) = T_tau(i,i)
        end do
        if(nshr.ne.0) then
          stressNew(km,ndir+1) = T_tau(1,2)
          if(nshr.ne.1) then
            stressNew(km, ndir+2) = T_tau(2,3)
            if(nshr.ne.2) then
              stressNew(km,ndir+3) = T_tau(1,3)
            end if
          end if
        end if


        ! Update the specific internal energy
        !
        stress_power = 0.d0
        do i = 1,ndir
          stress_power = stress_power +
     +           0.5*((StressOld(km,i)+StressNew(km,i))*
     +           StrainInc(km,i))
        end do
         
        select case (nshr)
         case(1)
            stress_power = stress_power + 
     +           0.5*((StressOld(km,ndir+1)+StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1))
         case(3)
            stress_power = stress_power + 
     +           0.5*(((StressOld(km,ndir+1) + StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1)) +
     +           ((StressOld(km,ndir+2)+ StressNew(km,ndir+2)) *
     +           StrainInc(km,ndir+2))+
     +           ((StressOld(km,ndir+3) + StressNew(km,ndir+3))*
     +           StrainInc(km,ndir+3)))
        end select

        xener=c10*(xI1/(Je**(two/three))-three)+(one/d1)*((Je-one)**two)
        enerInternNew(km) = xener/density(km)
	  
C        enerInternNew(km) = enerInternOld(km) + 
C     +        stress_power/density(km)
           
        enerInelasNew(km) = enerInelasOld(km) + 
     +        pwrinct/density(km)
           
           
       end do ! end loop over material points
	  endif
	  return
      end

C***********************************************************************
C	THE FOLLOWING SUBROUTINES ARE UTILITY ROUTINES
C**********************************************************************

      SUBROUTINE ONEM(A)

C	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C	3 BY 3 MATRIX [A]
      REAL*8 A(3,3)
      DATA ZERO/0.D0/
      DATA ONE/1.D0/

	  DO 1 I=1,3
	    DO 1 J=1,3
	      IF (I .EQ. J) THEN
              A(I,J) = 1.0
          ELSE
              A(I,J) = 0.0
          end if
1     CONTINUE

	  RETURN
	  END
C**********************************************************************
      SUBROUTINE mdet(A,DET)
 
C 	THIS SUBROUTINE CALCULATES THE DETERMINANT
C 	OF A 3 BY 3 MATRIX [A].

	  REAL*8  A(3,3), DET

	  DET =	  A(1,1)*A(2,2)*A(3,3) 
     +	    + A(1,2)*A(2,3)*A(3,1)
     +	    + A(1,3)*A(2,1)*A(3,2)
     +		- A(3,1)*A(2,2)*A(1,3)
     +		- A(3,2)*A(2,3)*A(1,1)
     +		- A(3,3)*A(2,1)*A(1,2)

	  RETURN
	  END
C**********************************************************************
      SUBROUTINE GetRandomBetweenRange(minValue, maxValue, randomValue)
      REAL*8 minValue, maxValue
      REAL*8 randomValue
      REAL*8 tempRand

      CALL RANDOM_NUMBER(tempRand)
      randomValue = minValue + tempRand * (maxValue - minValue)
	  
	  return
      END
	  
C**********************************************************************
      SUBROUTINE sprindir(stress,eigVal,tetha)
	  use Globalparameters

	  REAL*8 stress(4),eigVal(2,2),tetha,sigmax,sigmin
	  real*8 rotation_matrix(2,2),stressTensor(2,2)
	  real*8 Trotation_matrix(2,2)

        if (stress(1) /= stress(2)) then
            tetha = half*atan2((two*stress(4)/(stress(1)-stress(2))),one)
		else
            tetha = pi/four
        end if
		if (stress(4) .eq. zero) then
		    tetha = zero
		endif
        		
		sigmax = (stress(1)+stress(2))/two +sqrt((half*(stress(1)-stress(2)))**2+stress(4)**2)
		sigmin = (stress(1)+stress(2))/two -sqrt((half*(stress(1)-stress(2)))**2+stress(4)**2)
		! print *, sigmax
		rotation_matrix(1,1) = cos(tetha)
		rotation_matrix(1,2) = sin(tetha)
		rotation_matrix(2,1) = -sin(tetha)
		rotation_matrix(2,2) = cos(tetha)
		Trotation_matrix(1,1) = cos(tetha)
		Trotation_matrix(1,2) = -sin(tetha)
		Trotation_matrix(2,1) = sin(tetha)
		Trotation_matrix(2,2) = cos(tetha)	
        stressTensor(1,1) = stress(1)
        stressTensor(1,2) = stress(4)	
        stressTensor(2,1) = stress(4)
        stressTensor(2,2) = stress(2)						
		eigVal = matmul(rotation_matrix,matmul(stressTensor,Trotation_matrix))
					
	  RETURN
	  END
C**********************************************************************	  
      SUBROUTINE point_in_triangle(px, py, r, x1, y1, x2, y2, x3, y3)
	  use Globalparameters
        real*8 px, py, x1, y1, x2, y2, x3, y3
        real*8 A, A1, A2, A3,isInside
        integer r

        ! Calculate the area of the full triangle
        A = half * abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))

        ! Calculate the area of the sub-triangles with the point
        A1 = half * abs(px*(y2-y3) + x2*(y3-py) + x3*(py-y2))
        A2 = half * abs(x1*(py-y3) + px*(y3-y1) + x3*(y1-py))
        A3 = half * abs(x1*(y2-py) + x2*(py-y1) + px*(y1-y2))

        ! Check if the sum of sub-triangle areas equals the full triangle area
        isInside = abs(A - (A1 + A2 + A3))
		if (isInside .le. 1.0e-6) then
		    r = 1
		else
		    r = 0
		endif
      RETURN
	  end
C**********************************************************************	  
      SUBROUTINE point_in_quad(px, py, r, x1, y1, x2, y2, x3, y3, x4, y4)
        real*8 px, py
        real*8 x1, y1, x2, y2, x3, y3, x4, y4
		integer r, r1, r2, r3, r4

        ! Check if the point is inside either of the two possible triangles
        call point_in_triangle(px, py, r1, x1, y1, x2, y2, x3, y3)
        call point_in_triangle(px, py, r2, x1, y1, x3, y3, x4, y4)
        call point_in_triangle(px, py, r3, x2, y2, x3, y3, x4, y4)
        call point_in_triangle(px, py, r4, x1, y1, x2, y2, x4, y4)
		
		if ((r1 .eq. 1) .or. (r2.eq.1) .or. (r3.eq.1) .or. (r4.eq.1)) then
		    r = 1
		else
		    r = 0
		endif
				   
      return
	  end