!******************************************************************************* 
!                              INTEL CONFIDENTIAL 
!   Copyright(C) 2007-2008 Intel Corporation. All Rights Reserved. 
!   The source code contained  or  described herein and all documents related to 
!   the source code ("Material") are owned by Intel Corporation or its suppliers 
!   or licensors.  Title to the  Material remains with  Intel Corporation or its 
!   suppliers and licensors. The Material contains trade secrets and proprietary 
!   and  confidential  information of  Intel or its suppliers and licensors. The 
!   Material  is  protected  by  worldwide  copyright  and trade secret laws and 
!   treaty  provisions. No part of the Material may be used, copied, reproduced, 
!   modified, published, uploaded, posted, transmitted, distributed or disclosed 
!   in any way without Intel's prior express written permission. 
!   No license  under any  patent, copyright, trade secret or other intellectual 
!   property right is granted to or conferred upon you by disclosure or delivery 
!   of the Materials,  either expressly, by implication, inducement, estoppel or 
!   otherwise.  Any  license  under  such  intellectual property  rights must be 
!   express and approved by Intel in writing.
! 
!*******************************************************************************
!  This example gives the solution of initial value problem for the Van der 
!  Pol equation: 
!
!        y”-1.d6*[(1-y*y)*y’+1.d6*y=0,  0<t<160,    y(0)=2,  y’(0)=0. 
! 
!*******************************************************************************

      PROGRAM ODE_EXAMPLE_F

      IMPLICIT NONE
      
      INTEGER n, ierr, i
! It is higly recommended to declare ipar array of size 128 
! for compatibility with future versions of ODE solvers
      INTEGER kd(2), ipar(128)
      DOUBLE PRECISION t, t_end, h, hm, ep, tr
! As ODE system has size n=2, than the size of dpar array is equal to 
! max{13*n,(7+2*n)*n}=max{26,22}=26. More details on dpar array can be 
! found in the Manual
      DOUBLE PRECISION y(2), dpar(26)
      EXTERNAL rhs_v_d_p, jacmat_v_d_p
      REAL time_begin, time_end

! global parameter settings suitable for all 6 dodesol routines  
! minimal step size for the methods
         hm=1.d-12
! relative tolerance. The code cannot guarantee the requested accuracy for ep<1.d-9
         ep=1.d-6
! absolute tolerance
         tr=1.d-3

c****************************** dodesol ********************************
! Please don't forget to initialize ipar array with zeros before the first call 
! to dodesol routines
         DO i=1,128
            ipar(i)=0
         END DO  
         
         t=0.d0
         h=1.d-7       
         
! setting size of the system n, end of integration interval t_end, and initial 
! value y at t=0
         CALL example_v_d_p(n,t_end,y)

         CALL CPU_TIME(time_begin)
! universal solver
         CALL dodesol(ipar,n,t,t_end,y,rhs_v_d_p,jacmat_v_d_p,
     &                h,hm,ep,tr,dpar,kd,ierr)

         CALL CPU_TIME(time_end)

         IF(ierr.ne.0) THEN
            PRINT*,'========================'
            PRINT*,'DODESOL FORTRAN example FAILED'
            PRINT*,'dodesol routine exited with error code',ierr  
            STOP 1
         END IF
   
         PRINT*
         PRINT*, 'dodesol results'
         PRINT*
         PRINT*, 'ipar(2)=',ipar(2),', ipar(4)=',ipar(4)
         PRINT*, 't=',t
         PRINT*, 'Solution','  y1=',y(1),'  y2=',y(2)
         PRINT*, '-----------------------------------------------------'
         PRINT*, 'CPU time=',time_end-time_begin,' seconds'
         PRINT*, '====================================================='
         PRINT*
         IF(dabs(y(1)-1.878d0)+dabs(y(2)+0.7436d0).gt.1.d-2) THEN
            PRINT*,'Solution seems to be inaccurate. Probably, ',
     &                                               'example FAILED...'
            STOP 1
      	 END IF
		 
c****************************** dodesol_rkm9st *****************************
! Please don't forget to initialize ipar array with zeros before the first call 
! to dodesol routines
         DO i=1,128
            ipar(i)=0
         END DO
         
         t=0.d0
         h=1.d-7
         
! setting size of the system n, end of integration interval t_end, and initial 
! value y at t=0
         CALL example_v_d_p(n,t_end,y)

         CALL CPU_TIME(time_begin)
! explicit solver
         CALL dodesol_rkm9st(ipar,n,t,t_end,y,rhs_v_d_p,h,hm,ep,tr,
     &                       dpar,ierr)
         
         CALL CPU_TIME(time_end)
       
         IF(ierr.ne.0) THEN
            PRINT*,'========================'
            PRINT*,'DODESOL FORTRAN example FAILED'
            PRINT*,'dodesol_rkm9st routine exited with error code',ierr
            STOP 1
         END IF
   
         PRINT*
         PRINT*, 'dodesol_rkm9st results'
         PRINT*
         PRINT*, 't=',t
         PRINT*, 'Solution','  y1=',y(1),'  y2=',y(2)
         PRINT*, '-----------------------------------------------------'
         PRINT*, 'CPU time=',time_end-time_begin,' seconds'
         PRINT*, '====================================================='
         PRINT*
         IF(dabs(y(1)-1.878d0)+dabs(y(2)+0.7436d0).gt.1.d-2) THEN
      	    PRINT*,'Solution seems to be inaccurate. Probably, ',
     &                                               'example FAILED...'
            STOP 1
      	 END IF
		 
c****************************** dodesol_mk52lfn ********************************
! Please don't forget to initialize ipar array with zeros before the first call 
! to dodesol routines
         DO i=1,128
            ipar(i)=0
         END DO
         
         t=0.d0
         h=1.d-7
         
! setting size of the system n, end of integration interval t_end, and initial 
! value y at t=0
         CALL example_v_d_p(n,t_end,y)

         CALL CPU_TIME(time_begin)
! implicit solver with automatic numerical Jacobi matrix computations
         CALL dodesol_mk52lfn(ipar,n,t,t_end,y,rhs_v_d_p,h,hm,ep,tr,
     &                        dpar,kd,ierr)
         
         CALL CPU_TIME(time_end)
         
         IF(ierr.ne.0) THEN
            PRINT*,'========================'
            PRINT*,'DODESOL FORTRAN example FAILED'
            PRINT*,'dodesol_mk52lfn routine exited with error code',ierr
            STOP 1
         END IF
   
         PRINT*
         PRINT*, 'dodesol_mk52lfn results'
         PRINT*
         PRINT*, 't=',t
         PRINT*, 'Solution','  y1=',y(1),'  y2=',y(2)
         PRINT*, '-----------------------------------------------------'
         PRINT*, 'CPU time=',time_end-time_begin,' seconds'
         PRINT*, '====================================================='
         PRINT*
         IF(dabs(y(1)-1.878d0)+dabs(y(2)+0.7436d0).gt.1.d-2) THEN
      	    PRINT*,'Solution seems to be inaccurate. Probably, ',
     &                                               'example FAILED...'
            STOP 1
      	 END IF
		 
c****************************** dodesol_mk52lfa ********************************
! Please don't forget to initialize ipar array with zeros before the first call 
! to dodesol routines
         DO i=1,128
            ipar(i)=0
         END DO
         
         t=0.d0
         h=1.d-7
         
! setting size of the system n, end of integration interval t_end, and initial 
! value y at t=0
         CALL example_v_d_p(n,t_end,y)

         CALL CPU_TIME(time_begin)
! implicit solver with user-defined Jacobi matrix computations
         CALL dodesol_mk52lfa(ipar,n,t,t_end,y,rhs_v_d_p,jacmat_v_d_p,
     &                        h,hm,ep,tr,dpar,kd,ierr)

         CALL CPU_TIME(time_end)
         
         IF(ierr.ne.0) THEN
            PRINT*,'========================'
            PRINT*,'DODESOL FORTRAN example FAILED'
            PRINT*,'dodesol_mk52lfa routine exited with error code',ierr
            STOP 1
         END IF
   
         PRINT*
         PRINT*, 'dodesol_mk52lfa results'
         PRINT*
         PRINT*, 't=',t
         PRINT*, 'Solution','  y1=',y(1),'  y2=',y(2)
         PRINT*, '-----------------------------------------------------'
         PRINT*, 'CPU time=',time_end-time_begin,' seconds'
         PRINT*, '====================================================='
         PRINT*
         IF(dabs(y(1)-1.878d0)+dabs(y(2)+0.7436d0).gt.1.d-2) THEN
      	    PRINT*,'Solution seems to be inaccurate. Probably, ',
     &                                               'example FAILED...'
            STOP 1
      	 END IF
		 
c****************************** dodesol_rkm9mkn ********************************
! Please don't forget to initialize ipar array with zeros before the first call 
! to dodesol routines
         DO i=1,128
            ipar(i)=0
         END DO
         
         t=0.d0
         h=1.d-7
         
! setting size of the system n, end of integration interval t_end, and initial 
! value y at t=0
         CALL example_v_d_p(n,t_end,y)

         CALL CPU_TIME(time_begin)
! hybrid solver with automatic numerical Jacobi matrix computations
         CALL dodesol_rkm9mkn(ipar,n,t,t_end,y,rhs_v_d_p,h,hm,ep,tr,
     &                        dpar,kd,ierr)
         
         CALL CPU_TIME(time_end)
        
         IF(ierr.ne.0) THEN
            PRINT*,'========================'
            PRINT*,'DODESOL FORTRAN example FAILED'
            PRINT*,'dodesol_rkm9mkn routine exited with error code',ierr
            STOP 1
         END IF
   
         PRINT*
         PRINT*, 'dodesol_rmk9mkn results'
         PRINT*
         PRINT*, 't=',t
         PRINT*, 'Solution','  y1=',y(1),'  y2=',y(2)
         PRINT*, '-----------------------------------------------------'
         PRINT*, 'CPU time=',time_end-time_begin,' seconds'
         PRINT*, '====================================================='
         PRINT*
         IF(dabs(y(1)-1.878d0)+dabs(y(2)+0.7436d0).gt.1.d-2) THEN
      	    PRINT*,'Solution seems to be inaccurate. Probably, ',
     &                                               'example FAILED...'
            STOP 1
      	 END IF
		 
c****************************** dodesol_rkm9mka ********************************
! Please don't forget to initialize ipar array with zeros before the first call 
! to dodesol routines
         DO i=1,128
            ipar(i)=0
         END DO
         
         t=0.d0
         h=1.d-7
         
! setting size of the system n, end of integration interval t_end, and initial 
! value y at t=0
         CALL example_v_d_p(n,t_end,y)

         CALL CPU_TIME(time_begin)
! hybrid solver with user-defined Jacobi matrix computations
         CALL dodesol_rkm9mka(ipar,n,t,t_end,y,rhs_v_d_p,jacmat_v_d_p,
     &                        h,hm,ep,tr,dpar,kd,ierr)

         CALL CPU_TIME(time_end)
         
         IF(ierr.ne.0) THEN
            PRINT*,'========================'
            PRINT*,'DODESOL FORTRAN example FAILED'
            PRINT*,'dodesol_rkm9mka routine exited with error code',ierr
            STOP 1
         END IF
   
         PRINT*
         PRINT*, 'dodesol_rkm9mka results'
         PRINT*
         PRINT*, 't=',t
         PRINT*, 'Solution','  y1=',y(1),'  y2=',y(2)
         PRINT*, '-----------------------------------------------------'
         PRINT*, 'CPU time=',time_end-time_begin,' seconds'
         PRINT*, '====================================================='
         PRINT*
      	   IF(dabs(y(1)-1.878d0)+dabs(y(2)+0.7436d0).gt.1.d-2) THEN
      	    PRINT*,'Solution seems to be inaccurate. Probably, ',
     &                                               'example FAILED...'
            STOP 1
      	 END IF        
      	 PRINT*, '========================'
         PRINT*, 'DODESOL FORTRAN example successfully PASSED through',
     &                                    ' all steps of computations'
         PRINT*
          
      STOP 0
      END

c********************** Example for Van der Pol equations **************  
      SUBROUTINE example_v_d_p(n,t_end,y)
! The routine initializes the size of the system n, the end of 
! integration interval t_end, and inital data y at t=0.0 
      IMPLICIT NONE
      
      INTEGER n
      DOUBLE PRECISION t_end,y(*)
            
         n=2 
         t_end=160.d0

         y(1)=2.d0
         y(2)=0.d0

      RETURN
      END

c******************* Right hand side of Van der Pol equations ******************  
      SUBROUTINE rhs_v_d_p(n,t,y,f)

      IMPLICIT NONE
      
      INTEGER n
      DOUBLE PRECISION t,y(*),f(*)
 
         f(1)=y(2)
         f(2)=1.d6*((1.d0-y(1)*y(1))*y(2)-y(1))

      RETURN
      END

c************* analytical Jacoby matrix for Van der Pol equations **************
      SUBROUTINE jacmat_v_d_p(n,t,y,a)
      
      IMPLICIT NONE 
      
      INTEGER n
      DOUBLE PRECISION t,y(*),a(n,*)
      
         a(1,1)=0.d0
         a(1,2)=1.d0
         a(2,1)=-1.d6*(1.d0+2.d0*y(1)*y(2))
         a(2,2)= 1.d6*(1.d0-y(1)* y(1))

      RETURN
      END

C************************* End of Fortran code example *************************
