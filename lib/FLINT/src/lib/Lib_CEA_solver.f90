module FLINT_CEA_solver

contains


      SUBROUTINE CEA_solve(p_in, T_in, rhoi_in, T_out, rhoi_out)
!***********************************************************************
! CALCULATE EQUILIBRIUM COMPOSITION AND PROPERTIES.
!***********************************************************************
      use FLINT_CEA_data
      use FLINT_CEA_setup, only: CEA_initialize_local
      use FLINT_Lib_Thermodynamic, only: COMP_THERMO_QUANTS
      IMPLICIT NONE
      REAL*8, DIMENSION(:), INTENT(IN) :: rhoi_in
      REAL*8, INTENT(IN) :: p_in, T_in
      REAL*8, INTENT(OUT) :: T_out
      REAL*8, DIMENSION(:), INTENT(OUT) :: rhoi_out
! LOCAL VARIABLES
      CHARACTER*12 ae,cmp(MAXEL)
      CHARACTER*16 amb
      LOGICAL cpcalc,i2many,newcom,reduce
      INTEGER i,il,ilamb,ilamb1,inc,ipr,iq2,iter,ix,ixsing,iz,j,ja,jb, &
              jbx,jc,jcondi,jcons,jdelg,jex,jj,jkg,jneg,jsw,k,kc,kg,kk,&
              kneg,l,lc,lcs(MAXEL),le,lelim,lk,ll,lncvg,ls,lsing, &
              lz,maxitn,ncvg,njc,nn,numb,kmat
      INTEGER IABS
      REAL*8 aa,ambda,ambda1,bigen,bigneg,delg,dlnt,dpie,ensol,esize, &
             gap,gasfrc,pie,pisave(MAXMAT-2),siz9,sizeg,smalno,smnol, &
             sum,sum1,szgj,tem,tmelt,tsize,ween,xi,xln,xsize,xx(MAXMAT)
      REAL*8 DABS,DEXP,DLOG,DMAX1, Tln, Tm
      REAL*8 Vv, pp, tt, Hsub0
      REAL*8 B0p(MAXEL, 2)
      REAL*8 Enln(MAXNGC), Enn, En(MAXNGC, NCOL), Ennl, Sumn, B0(MAXEL), G(MAXMAT, MAXMAT+1)
      REAL(8) :: Cp(MAXNGC), H0(MAXNGC), Mu(MAXNGC), S(MAXNGC)
      LOGICAL Convg, Gonly

      DATA smalno/1.E-6/,smnol/ - 13.815511/

      call CEA_initialize_local(T_in, rhoi_in, Vv, B0p, Hsub0)
      B0 = B0p(:,2)

      Pp = p_in
      Tt = T_in

      ! INITIAL ESTIMATES
      Npr = 0
      Enn = .1D0
      Ennl = -2.3025851
      Sumn = Enn
      xi = dble(Ng)
      if ( xi==0.d0 ) xi = 1.
      xi = Enn/xi
      xln = DLOG(xi)
      do i = 1,Nc
        j = Ng + i
        En(j,1) = 0.D0
        Enln(j) = 0.D0
      enddo
      do j = 1,Ng
        En(j,1) = xi
        Enln(j) = xln
      enddo

      ixsing = 0
      lsing = 0
      jsw = 0
      jdelg = 0
      maxitn = 50
      ncvg = 0
      lncvg = 3*Nlm
      reduce = .FALSE.
      siz9 = Size - 9.2103404D0
      tsize = Size
      xsize = Size + 6.90775528D0
      IF ( Trace.NE.0. ) THEN
        maxitn = maxitn + Ngc/2
        xsize = -DLOG(Trace)
        IF ( xsize.LT.Size ) xsize = Size + .1
      ENDIF
      IF ( xsize.GT.80. ) xsize = 80.D0
      esize = MIN(80.D0,xsize+6.90775528D0)
      jcons = 0
      pie = 0.
      i2many = .FALSE.
      Pderiv = .FALSE.
      Convg = .FALSE.
      numb = 0
      cpcalc = .TRUE.
      GOTO 400
      k = 1
 100  j = Jcond(k)
      jc = j - Ng
      kg = -Ifz(jc)
      DO i = 1,9
        kg = kg + 1
        kc = jc + kg
        IF ( Tt.LE.Temp(2,kc) ) THEN
          IF ( kg.NE.0 ) THEN
            Jcond(k) = j + kg
            En(j+kg,Npt) = En(j,Npt)
            En(j,Npt) = 0.
          ENDIF
          GOTO 300
        ELSEIF ( kc.GE.Nc.OR.Ifz(kc+1).LE.Ifz(kc) ) THEN
          GOTO 200
        ENDIF
      ENDDO
 200  IF ( .NOT.Tp ) THEN
        Tt = Temp(2,kc) - 10.D0
        k = 1
        GOTO 100
      ENDIF
      !WRITE (*,99028) Prod(j)
      En(j,Npt) = 0.D0
      Enln(j) = 0.D0
      Deln(j) = 0.D0
      DO i = k,Npr
        Jcond(i) = Jcond(i+1)
      ENDDO
      Npr = Npr - 1
 300  k = k + 1
      IF ( k.LE.Npr ) GOTO 100
 400  Tln = DLOG(Tt)
      IF ( Vol ) Pp = Rr*Enn*Tt/Vv
      CALL COMP_THERMO_QUANTS(Tt,Ng,Cp,H0,S)
      Tm = DLOG(Pp/Enn)
      le = Nlm
      IF ( Lsave.NE.0.AND.Nlm.NE.Lsave ) THEN
        tem = EXP(-tsize)
        DO i = Lsave + 1,Nlm
          DO j = 1,Ng
            IF ( A(i,j).NE.0. ) THEN
              En(j,Npt) = tem
              Enln(j) = -tsize
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      ls = Nlm
      lelim = 0
      lz = ls
      IF ( Ions ) lz = ls - 1
! BEGIN ITERATION
 500  IF ( cpcalc ) THEN
        Cpsum = 0.D0
        DO j = 1,Ng
          Cpsum = Cpsum + En(j,Npt)*Cp(j)
        ENDDO
        IF ( Npr.NE.0 ) THEN
          DO k = 1,Npr
            j = Jcond(k)
            Cpsum = Cpsum + En(j,Npt)*Cp(j)
          ENDDO
          cpcalc = .FALSE.
        ENDIF
      ENDIF
      numb = numb + 1
      CALL MATRIX
      iq2 = Iq1 + 1
      IF ( Convg ) Imat = Imat - 1
      ! IF ( .true. ) THEN
      !   IF ( .true. ) THEN
      !     WRITE (*,99004) numb
      !   ELSE
      !     IF ( .NOT.Pderiv ) WRITE (*,99002)
      !     IF ( Pderiv ) WRITE (*,99003)
      !   ENDIF
      !   kmat = Imat + 1
      !   DO i = 1,Imat
      !     WRITE (*,99006) (G(i,k),k=1,kmat)
      !   ENDDO
      !  if (numb==2) stop
      ! ENDIF
      Msing = 0
      CALL GAUSS
      IF ( Msing.EQ.0 ) THEN
        ! IF ( .true. ) THEN
        ! WRITE (*,*) (cmp(k),k=1,le)
        ! WRITE (*,*) (X(i),i=1,Imat)
        ! ENDIF
        IF ( .NOT.Convg ) THEN
! OBTAIN CORRECTIONS TO THE ESTIMATES
          IF ( Vol ) X(iq2) = X(Iq1)
          IF ( Tp ) X(iq2) = 0.
          dlnt = X(iq2)
          sum = X(Iq1)
          IF ( Vol ) THEN
            X(Iq1) = 0.
            sum = -dlnt
          ENDIF
          DO 520 j = 1,Ng
            IF ( lelim.NE.0 ) THEN
              Deln(j) = 0.
              DO i = lelim,ls
                IF ( A(i,j).NE.0. ) GOTO 520
              ENDDO
            ENDIF
            Deln(j) = -Mu(j) + H0(j)*dlnt + sum
            DO k = 1,Nlm
              Deln(j) = Deln(j) + A(k,j)*X(k)
            ENDDO
            IF ( pie.NE.0. ) Deln(j) = Deln(j) + A(ls,j)*pie
 520      CONTINUE
          IF ( Npr.NE.0 ) THEN
            DO k = 1,Npr
              j = Jcond(k)
              kk = Nlm + k
              Deln(j) = X(kk)
            ENDDO
          ENDIF
! CALCULATE CONTROL FACTOR,AMBDA
          ambda = 1.D0
          ambda1 = 1.D0
          ilamb = 0
          ilamb1 = 0
          sum = DMAX1(DABS(X(Iq1)),DABS(dlnt))
          sum = sum*5.
          DO j = 1,Ng
            IF ( Deln(j).GT.0. ) THEN
              IF ( (Enln(j)-Ennl+Size).LE.0. ) THEN
                sum1 = DABS(Deln(j)-X(Iq1))
                IF ( sum1.GE.siz9 ) THEN
                  sum1 = DABS(-9.2103404D0-Enln(j)+Ennl)/sum1
                  IF ( sum1.LT.ambda1 ) THEN
                    ambda1 = sum1
                    ilamb1 = j
                  ENDIF
                ENDIF
              ELSEIF ( Deln(j).GT.sum ) THEN
                sum = Deln(j)
                ilamb = j
              ENDIF
            ENDIF
          ENDDO
          IF ( sum.GT.2.D0 ) ambda = 2.D0/sum
          IF ( ambda1.LE.ambda ) THEN
            ambda = ambda1
            ilamb = ilamb1
          ENDIF
          ! IF ( .true. ) THEN
          ! ! INTERMEDIATE OUTPUT
          !   WRITE (*,99011) Tt,Enn,Ennl,Pp,Tm,ambda
          !   IF ( ambda.NE.1.D0 ) THEN
          !     amb = 'ENN'
          !     IF ( DABS(X(iq2)).GT.DABS(X(Iq1)) ) amb = 'TEMP'
          !     IF ( ilamb.NE.0 ) amb = Prod(ilamb)
          !     WRITE (*,99012) amb
          !   ENDIF
          !   !IF ( Vol ) WRITE (IOOUT,99013) Vv*.001D0
          !   WRITE (*,99014)
          !   DO j = 1,Ngc
          !    WRITE (*,99015) Prod(j),En(j,Npt),Enln(j),Deln(j),H0(j),S(j),H0(j) - S(j),Mu(j)
          !   ENDDO
          ! ENDIF
! APPLY CORRECTIONS TO ESTIMATES
          Totn(Npt) = 0.D0
          DO j = 1,Ng
            Enln(j) = Enln(j) + ambda*Deln(j)
          ENDDO
          DO 540 j = 1,Ng
            En(j,Npt) = 0.
            IF ( lelim.NE.0 ) THEN
              DO i = lelim,ls
                IF ( A(i,j).NE.0. ) GOTO 540
              ENDDO
            ENDIF
            IF ( (Enln(j)-Ennl+tsize).GT.0. ) THEN
              En(j,Npt) = DEXP(Enln(j))
              Totn(Npt) = Totn(Npt) + En(j,Npt)
            ENDIF
 540      CONTINUE
          IF ( Ions.AND.Elmt(Nlm).EQ.'E' ) THEN
            DO j = 1,Ng
              IF ( A(ls,j).NE.0..AND.En(j,Npt).EQ.0. ) THEN
                IF ( (Enln(j)-Ennl+esize).GT.0. ) THEN
                  En(j,Npt) = DEXP(Enln(j))
                  Totn(Npt) = Totn(Npt) + En(j,Npt)
                ENDIF
              ENDIF
            ENDDO
          ENDIF
          Sumn = Totn(Npt)
          IF ( Npr.NE.0 ) THEN
            DO k = 1,Npr
              j = Jcond(k)
              En(j,Npt) = En(j,Npt) + ambda*Deln(j)
              Totn(Npt) = Totn(Npt) + En(j,Npt)
            ENDDO
          ENDIF
          IF ( .NOT.Tp ) THEN
            Tln = Tln + ambda*dlnt
            Tt = DEXP(Tln)
            cpcalc = .TRUE.
            CALL COMP_THERMO_QUANTS(Tt,Ng,Cp,H0,S)
          ENDIF
          IF ( Vol ) THEN
            Enn = Sumn
            Ennl = DLOG(Enn)
            IF ( Vol ) Pp = Rr*Tt*Enn/Vv
          ELSE
            Ennl = Ennl + ambda*X(Iq1)
            Enn = DEXP(Ennl)
          ENDIF
          Tm = DLOG(Pp/Enn)
          IF ( Elmt(Nlm).EQ.'E' ) THEN
! CHECK ON REMOVING IONS
            DO j = 1,Ngc
              IF ( A(Nlm,j).NE.0. ) THEN
                IF ( En(j,Npt).GT.0. ) GOTO 560
              ENDIF
            ENDDO
            pie = X(Nlm)
            lelim = Nlm
            Nlm = Nlm - 1
            GOTO 500
          ENDIF
! TEST FOR CONVERGENCE
 560      IF ( numb.GT.maxitn ) THEN
            WRITE (*,*) maxitn,Npt
            IF ( Nc.EQ.0.OR.i2many ) GOTO 1500
            i2many = .TRUE.
            IF ( .NOT.Hp.OR.Npt.NE.1.OR.Tt.GT.100. ) THEN
              IF ( Npr.NE.1.OR.Enn.GT.1.E-4 ) GOTO 1500
! HIGH TEMPERATURE, INCLUDED CONDENSED CONDITION
              !WRITE (IOOUT,99020)
              Enn = .1
              Ennl = -2.3025851
              Sumn = Enn
              xi = Ng
              xi = Enn/xi
              xln = DLOG(xi)
              DO j = 1,Ng
                En(j,Npt) = xi
                Enln(j) = xln
              ENDDO
              j = Jcond(1)
              k = 1
              GOTO 1000
            ELSE
             ! WRITE (IOOUT,99008)
              GOTO 1500
            ENDIF
          ELSE
            sum = (X(Iq1)*Enn/Totn(Npt))
            IF ( DABS(sum).GT.0.5E-5 ) GOTO 500
            DO j = 1,Ng
              IF ( DABS(Deln(j))*En(j,Npt)/Totn(Npt).GT.0.5D-5 )GOTO 500
            ENDDO
            IF ( DABS(dlnt).GT.1.D-04 ) GOTO 500
            IF ( Npr.NE.0 ) THEN
              DO k = 1,Npr
                j = Jcond(k)
                IF ( DABS(Deln(j)/Totn(Npt)).GT.0.5D-5 ) GOTO 500
                IF ( En(j,Npt).LT.0. ) GOTO 700
              ENDDO
            ENDIF
            le = Nlm
            DO i = 1,Nlm
              IF ( DABS(B0(i)).GE.1.D-06 ) THEN
                sum = 0.
                DO j = 1,Ngc
                  sum = sum + En(j,Npt)*A(i,j)
                ENDDO
                IF ( DABS(B0(i)-sum).GT.Bcheck ) GOTO 500
              ENDIF
            ENDDO
            IF ( Trace.NE.0. ) THEN
              tsize = xsize
              tem = 1.
              IF ( numb.NE.1 ) THEN
                lk = lz
                IF ( Nlm.LT.lz ) lk = Nlm
                DO i = 1,lk
                  IF ( i.NE.lsing ) THEN
                    tem = 0.
                    IF ( X(i).NE.0. ) THEN
                      tem = DABS((pisave(i)-X(i))/X(i))
                      IF ( tem.GT..001 ) GOTO 565
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
 565          DO i = 1,Nlm
                pisave(i) = X(i)
              ENDDO
              IF ( tem.GT..001 ) GOTO 500
              IF ( Ions ) THEN
! CHECK ON ELECTRON BALANCE
                iter = 1
                IF ( pie.NE.0. ) THEN
                  le = Nlm + 1
                  X(le) = pie
                ENDIF
 566            sum1 = 0.
                sum = 0.
                pie = X(le)
                DO j = 1,Ng
                  IF ( A(ls,j).NE.0. ) THEN
                    En(j,Npt) = 0.
                    tem = 0.
                    IF ( Enln(j).GT.-87. ) tem = DEXP(Enln(j))
                    IF ( (Enln(j)-Ennl+tsize).GT.0..AND.Elmt(Nlm).EQ.'E' ) THEN
                      pie = 0.
                      En(j,Npt) = tem
                    ENDIF
                    aa = A(ls,j)*tem
                    sum = sum + aa
                    sum1 = sum1 + aa*A(ls,j)
                  ENDIF
                ENDDO
                IF ( sum1.NE.0. ) THEN
                  dpie = -sum/sum1
                  DO j = 1,Ng
                    IF ( A(ls,j).NE.0. ) Enln(j) = Enln(j) + A(ls,j)*dpie
                  ENDDO
                  IF ( DABS(dpie).GT..0001 ) THEN
                    X(le) = X(le) + dpie
                    iter = iter + 1
                    IF ( iter.LE.80 ) GOTO 566
                    !WRITE (IOOUT,99017)
                    GOTO 1500
                  ELSEIF ( Elmt(Nlm).EQ.'E'.AND.pie.NE.0. ) THEN
                    Nlm = Nlm - 1
                    newcom = .TRUE.
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ELSEIF ( .NOT.Pderiv ) THEN
! TEMPERATURE DERIVATIVES--CONVG=T, PDERIV=F
          Dlvtp(Npt) = 1. - X(Iq1)
          Cpr(Npt) = G(iq2,iq2)
          DO j = 1,Iq1
            Cpr(Npt) = Cpr(Npt) - G(iq2,j)*X(j)
          ENDDO
! PRESSURE DERIVATIVE--CONVG=T, PDERIV=T
          Pderiv = .TRUE.
          GOTO 500
        ELSE
          Dlvpt(Npt) = -1. + X(Iq1)
          IF ( Jliq/=0 ) THEN
            En(Jsol,Npt) = ensol
            Hsum(Npt) = Hsum(Npt) + En(Jliq,Npt)*(H0(Jliq)-H0(Jsol))
            Npr = Npr + 1
            Jcond(Npr) = Jliq
          ENDIF
          GOTO 1400
        ENDIF
! SINGULAR MATRIX
      ELSE
        IF ( Convg ) THEN
          !WRITE (IOOUT,99007)
          Dlvpt(Npt) = -1.
          Dlvtp(Npt) = 1.
          Cpr(Npt) = Cpsum
          GOTO 1400
        ELSE
          !WRITE (IOOUT,99009) numb,Msing
          lsing = Msing
          ixsing = ixsing + 1
          IF ( ixsing.LE.8 ) THEN
            xsize = 80.
            tsize = xsize
            IF ( Msing.GT.Nlm.AND.numb.LT.1.AND.Npr.GT.1.AND.jdelg.GT.0 ) THEN
              ween = 1000.
              j = 0
              DO 570 i = 1,Npr
                jcondi = Jcond(i)
                IF ( jcondi.NE.jdelg ) THEN
                  DO ll = 1,Nlm
                    IF ( A(ll,jdelg).NE.0.AND.A(ll,jcondi).NE.0. ) THEN
                      IF ( En(jcondi,Npt).LE.ween ) THEN
                        ween = En(jcondi,Npt)
                        j = jcondi
                        k = i
                      ENDIF
                      GOTO 570
                    ENDIF
                  ENDDO
                ENDIF
 570          CONTINUE
              IF ( j.GT.0 ) THEN
                !WRITE (IOOUT,99020)
                GOTO 1000
              ENDIF
            ELSEIF ( .NOT.Hp.OR.Npt.NE.1.OR.Nc.EQ.0.OR.Tt.GT.100. ) THEN
              IF ( ixsing.GE.3 ) THEN
                IF ( Msing.LT.Iq1 ) THEN
                  IF ( reduce.AND.Msing.LE.Nlm ) THEN
                    IF ( Nlm.LT.lelim ) GOTO 1500
                    !WRITE (IOOUT,99010) Npt,Elmt(Nlm)
                    Nlm = Nlm - 1
                    GOTO 500
                  ELSEIF ( Msing.LE.Nlm ) THEN
! FIND NEW COMPONENTS
                    IF ( .NOT.Ions ) GOTO 1100
                    IF ( Elmt(Nlm).NE.'E' ) GOTO 1100
                    DO j = 1,Ng
                      IF ( A(Nlm,j).NE.0. ) En(j,Npt) = 0.D0
                    ENDDO
                    pie = X(Nlm)
                    Nlm = Nlm - 1
                    IF ( Msing.GT.Nlm ) GOTO 500
                    GOTO 1100
                  ELSE
! REMOVE CONDENSED SPECIES TO CORRECT SINGULARITY
                    k = Msing - Nlm
                    j = Jcond(k)
                    IF ( j.NE.jcons ) THEN
                      jcons = j
                      GOTO 1000
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
              DO 575 jj = 1,Ng
                IF ( Ions ) THEN
                  IF ( Elmt(Nlm).NE.'E' ) THEN
                    IF ( A(ls,jj).NE.0. ) GOTO 575
                  ENDIF
                ENDIF
                IF ( En(jj,Npt).EQ.0. ) THEN
                  En(jj,Npt) = smalno
                  Enln(jj) = smnol
                ENDIF
 575          CONTINUE
              GOTO 500
            ELSE
              !WRITE (IOOUT,99008)
            ENDIF
          ENDIF
        ENDIF
        GOTO 1500
      ENDIF
! CALCULATE ENTROPY, CHECK ON DELTA S FOR SP PROBLEMS
 600  Ssum(Npt) = 0.
      DO j = 1,Ng
        Ssum(Npt) = Ssum(Npt) + En(j,Npt)*(S(j)-Enln(j)-Tm)
      ENDDO
      IF ( Npr.GT.0 ) THEN
        DO k = 1,Npr
          j = Jcond(k)
          Ssum(Npt) = Ssum(Npt) + En(j,Npt)*S(j)
        ENDDO
      ENDIF
      Convg = .TRUE.
! CONVERGENCE TESTS ARE SATISFIED, TEST CONDENSED SPECIES.
 700  ncvg = ncvg + 1
      IF ( ncvg.GT.lncvg ) THEN
! ERROR, SET TT=0
        !WRITE (IOOUT,99034) lncvg
        GOTO 1500
      ELSE
        DO il = 1,le
          xx(il) = X(il)
        ENDDO
        newcom = .FALSE.
        IF ( Npr.NE.0 ) THEN
          bigneg = 0.
          jneg = 0
          DO k = 1,Npr
            j = Jcond(k)
            IF ( En(j,Npt)*Cp(j).LE.bigneg ) THEN
              bigneg = En(j,Npt)*Cp(j)
              jneg = j
              kneg = k
            ENDIF
          ENDDO
          IF ( jneg.NE.0 ) THEN
            j = jneg
            k = kneg
            IF ( j.EQ.Jsol.OR.j.EQ.Jliq ) THEN
              Jsol = 0
              Jliq = 0
            ENDIF
            GOTO 1000
          ENDIF
        ENDIF
        IF ( Ngc.NE.Ng.OR.Tp ) THEN
          Ng = Ngc
          CALL COMP_THERMO_QUANTS(Tt,Ng,Cp,H0,S)
          Ng = Ngp1 - 1
          cpcalc = .TRUE.
          IF ( Ngc.EQ.Ng ) GOTO 750
          !CALL ALLCON
          IF ( Npr.NE.0.AND..NOT.Tp ) THEN
            gap = 50.
            DO 710 ipr = 1,Npr
              j = Jcond(ipr)
              IF ( j.NE.Jsol.AND.j.NE.Jliq ) THEN
                inc = j - Ng
                kg = -Ifz(inc)
                DO iz = 1,20
                  kg = kg + 1
                  kc = inc + kg
                  IF ( Tt.LE.Temp(2,kc) ) THEN
                    IF ( kg.NE.0 ) THEN
                      jkg = j + kg
                      IF ( IABS(kg).GT.1.OR.Prod(j).EQ.Prod(jkg) )GOTO 740
                      IF ( jkg.EQ.jsw ) GOTO 720
                      IF ( Tt.LT.Temp(1,inc)-gap.OR.Tt.GT.Temp(2,inc)+gap ) GOTO 740
                      GOTO 720
                    ENDIF
                    GOTO 710
                  ELSEIF ( Ifz(kc+1).LE.Ifz(kc) ) THEN
                    GOTO 710
                  ENDIF
                ENDDO
                IF ( Tt.GT.Temp(2,kc)*1.2D0 ) GOTO 1000
              ENDIF
 710        CONTINUE
          ENDIF
          sizeg = 0.
          szgj = 0.
          DO inc = 1,Nc
            j = inc + Ng
          ENDDO
          IF ( sizeg.EQ.0..AND.szgj.EQ.0. ) GOTO 750
          IF ( sizeg.NE.0. ) THEN
            j = jdelg
            GOTO 800
          ELSE
            !WRITE (IOOUT,99026) Prod(jcons)
            GOTO 1500
          ENDIF
 720      kk = MAX(0,kg)
          tmelt = Temp(kk+1,inc)
          Tt = tmelt
          Tln = DLOG(Tt)
          Jsol = MIN(j,jkg)
          Jliq = Jsol + 1
          En(jkg,Npt) = .5D0*En(j,Npt)
          En(j,Npt) = En(jkg,Npt)
          j = jkg
          GOTO 800
! WRONG PHASE INCLUDED FOR T INTERVAL, SWITCH EN
 740      En(jkg,Npt) = En(j,Npt)
          Jcond(ipr) = jkg
          En(j,Npt) = 0.
          jsw = j
          !IF ( Prod(j).NE.Prod(jkg).AND..NOT.Short ) WRITE (IOOUT,99023)Prod(j),Prod(jkg)
          j = jkg
          GOTO 900
        ENDIF
! CONVERGED WITH NO CONDENSED CHANGES.  IF BOTH SOLID & LIQ PRESENT,
! TEMPORARILY REMOVE LIQ TO PREVENT SINGULAR DERIVATIVE MATRIX.
 750    Sumn = Enn
        IF ( Jsol.NE.0 ) THEN
          ensol = En(Jsol,Npt)
          En(Jsol,Npt) = En(Jsol,Npt) + En(Jliq,Npt)
          Dlvtp(Npt) = 0.
          Cpr(Npt) = 0.
          Pderiv = .TRUE.
          DO k = 1,Npr
            IF ( Jcond(k).EQ.Jliq ) GOTO 760
          ENDDO
 760      DO i = k,Npr
            Jcond(i) = Jcond(i+1)
          ENDDO
          Npr = Npr - 1
        ENDIF
        GOTO 500
      ENDIF
! ADD CONDENSED SPECIES
 800  Npr = Npr + 1
      i = Npr
      DO ix = 2,Npr
        Jcond(i) = Jcond(i-1)
        i = i - 1
      ENDDO
      Jcond(1) = j
      !IF ( .NOT.Short ) WRITE (IOOUT,99027) Prod(j)
 900  inc = j - Ng
      Convg = .FALSE.
      IF ( Tp ) cpcalc = .FALSE.
      numb = -1
      GOTO 500
! REMOVE CONDENSED SPECIES
 1000 En(j,Npt) = 0.D0
      Deln(j) = 0.D0
      Enln(j) = 0.D0
      DO i = k,Npr
        Jcond(i) = Jcond(i+1)
      ENDDO
      !IF ( .NOT.Short ) WRITE (IOOUT,99028) Prod(j)
      Npr = Npr - 1
      DO i = 1,Nlm
        IF ( cmp(i).EQ.Prod(j) ) THEN
          numb = -1
          Convg = .FALSE.
          IF ( Tp ) cpcalc = .FALSE.
          GOTO 1100
        ENDIF
      ENDDO
      GOTO 900
 1100 newcom = .FALSE.
      nn = Nlm
      IF ( Elmt(Nlm).EQ.'E' ) nn = Nlm - 1
! FIND ORDER OF SPECIES FOR COMPONENTS - BIGGEST TO SMALLEST
      njc = 0
      DO lc = 1,nn
        lcs(lc) = 0
      ENDDO
 1200 bigen = -1.D-35
      DO j = 1,Ng
        IF ( En(j,Npt).GT.bigen ) THEN
          IF ( .NOT.Ions.OR.A(ls,j).EQ.0. ) THEN
            bigen = En(j,Npt)
            jbx = j
          ENDIF
        ENDIF
      ENDDO
      IF ( bigen.GT.0. ) THEN
        DO 1250 lc = 1,nn
          IF ( jbx.EQ.0 ) jbx = Jx(lc)
          IF ( A(lc,jbx).GT.smalno ) THEN
            IF ( njc.NE.0 ) THEN
              DO 1205 i = 1,njc
                l = lcs(i)
                IF ( l.EQ.lc ) GOTO 1250
                IF ( l.EQ.0 ) GOTO 1210
                j = Jcm(l)
                DO l = 1,nn
                  IF ( A(l,jbx).NE.A(l,j) ) GOTO 1205
                ENDDO
                GOTO 1250
 1205         CONTINUE
            ENDIF
 1210       DO i = 1,nn
              IF ( i.NE.lc ) THEN
                jex = Jx(i)
                IF ( DABS(A(lc,jbx)*A(i,jex)-A(lc,jex)*A(i,jbx)).LE.smalno ) GOTO 1250
              ENDIF
            ENDDO
            njc = njc + 1
            IF ( jbx.NE.Jcm(lc) ) newcom = .TRUE.
            Jcm(lc) = jbx
            lcs(njc) = lc
            GOTO 1300
          ENDIF
 1250   CONTINUE
 1300   En(jbx,Npt) = -En(jbx,Npt)
        IF ( njc.LT.nn ) GOTO 1200
      ENDIF
      DO j = 1,Ng
        En(j,Npt) = DABS(En(j,Npt))
      ENDDO
      IF ( newcom ) THEN
! SWITCH COMPONENTS
        DO lc = 1,nn
          jb = Jcm(lc)
          IF ( A(lc,jb).EQ.0. ) THEN
            jb = Jx(lc)
            Jcm(lc) = jb
          ENDIF
          tem = A(lc,jb)
          IF ( tem.NE.0. ) THEN
            pisave(lc) = H0(jb) - S(jb)
            IF ( jb.LE.Ng ) pisave(lc) = pisave(lc) + Enln(jb) + Tm
            cmp(lc) = Prod(jb)
! CALCULATE NEW COEFFICIENTS
            IF ( tem.NE.1. ) THEN
              B0(lc) = B0(lc)/tem
              B0p(lc,1) = B0p(lc,1)/tem
              B0p(lc,2) = B0p(lc,2)/tem
              DO j = 1,Nspx
                A(lc,j) = A(lc,j)/tem
              ENDDO
            ENDIF
            DO i = 1,nn
              IF ( A(i,jb).NE.0..AND.i.NE.lc ) THEN
                tem = A(i,jb)
                DO j = 1,Nspx
                  A(i,j) = A(i,j) - A(lc,j)*tem
                  IF ( DABS(A(i,j)).LT.1.E-5 ) A(i,j) = 0.
                ENDDO
                B0(i) = B0(i) - B0(lc)*tem
                B0p(i,1) = B0p(i,1) - B0p(lc,1)*tem
                B0p(i,2) = B0p(i,2) - B0p(lc,2)*tem
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      IF ( Msing.NE.0 ) THEN
! SWITCH ORDER OF MSING AND NLM COMPONENTS
        reduce = .TRUE.
        lelim = Nlm
        lsing = Nlm
        IF ( Msing.NE.Nlm ) THEN
          DO j = 1,Nspx
            aa = A(Msing,j)
            A(Msing,j) = A(Nlm,j)
            A(Nlm,j) = aa
          ENDDO
          ja = Jcm(Msing)
          Jcm(Msing) = Jcm(Nlm)
          Jcm(Nlm) = ja
          ae = cmp(Msing)
          cmp(Msing) = cmp(Nlm)
          cmp(Nlm) = ae
          ae = Elmt(Msing)
          Elmt(Msing) = Elmt(Nlm)
          Elmt(Nlm) = ae
          ja = Jx(Msing)
          Jx(Msing) = Jx(Nlm)
          Jx(Nlm) = ja
          aa = Atwt(Msing)
          Atwt(Msing) = Atwt(Nlm)
          Atwt(Nlm) = aa
          aa = B0(Msing)
          B0(Msing) = B0(Nlm)
          B0(Nlm) = aa
          aa = pisave(Msing)
          pisave(Msing) = pisave(Nlm)
          pisave(Nlm) = aa
          DO i = 1,2
            aa = B0p(Msing,i)
            B0p(Msing,i) = B0p(Nlm,i)
            B0p(Nlm,i) = aa
          ENDDO
        ENDIF
      ELSEIF ( .NOT.newcom.AND.Trace.EQ.0. ) THEN
        GOTO 600
      ENDIF
      Msing = 0
      tsize = xsize
      GOTO 500
 1400 Ttt(Npt) = Tt
      Ppp(Npt) = Pp
      Hsum(Npt) = Hsum(Npt)*Tt
      Wm(Npt) = 1./Enn
      gasfrc = Enn/Totn(Npt)
      t_out = Tt
      rhoi_out = en(1:Nreac,1)*Mw(1:Nreac)
      IF ( Trace.NE.0. ) THEN
        DO 1450 j = 1,Ng
          IF ( lelim.NE.0 ) THEN
            DO i = lelim,ls
              IF ( A(i,j).NE.0. ) GOTO 1450
            ENDDO
          ENDIF
          IF ( Enln(j).GT.-87. ) En(j,Npt) = DEXP(Enln(j))
 1450   CONTINUE
      ENDIF
      Npt = Npt + 1
 1500 Tt = 0.
      Npt = Npt - 1
      !WRITE (IOOUT,99035) Npt
 1600 Lsave = Nlm
      Nlm = ls
      IF ( Npr.GT.0 ) Gonly = .FALSE.

contains

      SUBROUTINE MATRIX
!***********************************************************************
! SET UP ITERATION OR DERIVATIVE MATRIX.
!***********************************************************************
      use FLINT_CEA_data
      IMPLICIT NONE
! LOCAL VARIABLES
      INTEGER i,iq,iq2,iq3,isym,j,k,kk,kmat
      REAL*8 energyl,f,h,ss,sss,term,term1
      SAVE energyl,f,h,i,iq,iq2,iq3,isym,j,k,kk,kmat,ss,sss,term,term1

      iq = Nlm + Npr
      Iq1 = iq + 1
      iq2 = Iq1 + 1
      iq3 = iq2 + 1
      kmat = iq3
      IF ( .NOT.Convg.AND.Tp ) kmat = iq2
      Imat = kmat - 1
! CLEAR MATRIX STORAGES TO ZERO
      DO i = 1,Imat
        DO k = 1,kmat
          G(i,k) = 0.0D0
        ENDDO
      ENDDO
      G(iq2,Iq1) = 0.D0
      sss = 0.D0
      Hsum(Npt) = 0.D0
! BEGIN SET-UP OF ITERATION OR DERIVATIVE MATRIX
      DO j = 1,Ng
        Mu(j) = H0(j) - S(j) + Enln(j) + Tm
        IF ( En(j,Npt).NE.0.D0 ) THEN
          h = H0(j)*En(j,Npt)
          f = Mu(j)*En(j,Npt)
          ss = h - f
          term1 = h
          IF ( kmat.EQ.iq2 ) term1 = f
          DO i = 1,Nlm
            IF ( A(i,j).NE.0. ) THEN
              term = A(i,j)*En(j,Npt)
              DO k = i,Nlm
                G(i,k) = G(i,k) + A(k,j)*term
              ENDDO
              G(i,Iq1) = G(i,Iq1) + term
              G(i,iq2) = G(i,iq2) + A(i,j)*term1
              IF ( .NOT.(Convg.OR.Tp) ) THEN
                G(i,iq3) = G(i,iq3) + A(i,j)*f
              ENDIF
            ENDIF
          ENDDO
          IF ( kmat.NE.iq2 ) THEN
            IF ( Convg.OR.Hp ) THEN
              G(iq2,iq2) = G(iq2,iq2) + H0(j)*h
              IF ( .NOT.Convg ) THEN
                G(iq2,iq3) = G(iq2,iq3) + H0(j)*f
                G(Iq1,iq3) = G(Iq1,iq3) + f
              ENDIF
            ELSE
              G(iq2,Iq1) = G(iq2,Iq1) + ss
              G(iq2,iq2) = G(iq2,iq2) + H0(j)*ss
              G(iq2,iq3) = G(iq2,iq3) + Mu(j)*ss
              G(Iq1,iq3) = G(Iq1,iq3) + f
            ENDIF
          ENDIF
          G(Iq1,iq2) = G(Iq1,iq2) + term1
        ENDIF
      ENDDO
! CONDENSED SPECIES
      IF ( Npr.NE.0 ) THEN
        DO k = 1,Npr
          j = Jcond(k)
          kk = Nlm + k
          Mu(j) = H0(j) - S(j)
          DO i = 1,Nlm
            G(i,kk) = A(i,j)
            G(i,kmat) = G(i,kmat) - A(i,j)*En(j,Npt)
          ENDDO
          G(kk,iq2) = H0(j)
          G(kk,kmat) = Mu(j)
          Hsum(Npt) = Hsum(Npt) + H0(j)*En(j,Npt)
        ENDDO
      ENDIF
      sss = sss + G(iq2,Iq1)
      Hsum(Npt) = Hsum(Npt) + G(Iq1,iq2)
      G(Iq1,Iq1) = Sumn - Enn
! REFLECT SYMMETRIC PORTIONS OF THE MATRIX
      isym = Iq1
      IF ( Hp.OR.Convg ) isym = iq2
      DO i = 1,isym
!CDIR$ IVDEP
        DO j = i,isym
          G(j,i) = G(i,j)
        ENDDO
      ENDDO
! COMPLETE THE RIGHT HAND SIDE
      IF ( .NOT.Convg ) THEN
        DO i = 1,Nlm
          G(i,kmat) = G(i,kmat) + B0(i) - G(i,Iq1)
        ENDDO
        G(Iq1,kmat) = G(Iq1,kmat) + Enn - Sumn
! COMPLETE ENERGY ROW AND TEMPERATURE COLUMN
        IF ( kmat.NE.iq2 ) THEN
          IF ( Hp ) energyl = Hsub0/Tt - Hsum(Npt)
          G(iq2,iq3) = G(iq2,iq3) + energyl
          G(iq2,iq2) = G(iq2,iq2) + Cpsum
        ENDIF
      ELSE
        IF ( Pderiv ) THEN
! PDERIV = .TRUE.-- SET UP MATRIX TO SOLVE FOR DLVPT
          G(Iq1,iq2) = Enn
          DO i = 1,iq
            G(i,iq2) = G(i,Iq1)
          ENDDO
        ENDIF
        G(iq2,iq2) = G(iq2,iq2) + Cpsum
      ENDIF
      IF ( Vol.AND..NOT.Convg ) THEN
! CONSTANT VOLUME MATRIX
        IF ( kmat.EQ.iq2 ) THEN
          DO i = 1,iq
            G(i,Iq1) = G(i,iq2)
          ENDDO
        ELSE
!CDIR$ IVDEP
          DO i = 1,iq
            G(Iq1,i) = G(iq2,i) - G(Iq1,i)
            G(i,Iq1) = G(i,iq2) - G(i,Iq1)
            G(i,iq2) = G(i,iq3)
          ENDDO
          G(Iq1,Iq1) = G(iq2,iq2) - G(Iq1,iq2) - G(iq2,Iq1)
          G(Iq1,iq2) = G(iq2,iq3) - G(Iq1,iq3)
          IF ( Hp ) G(Iq1,iq2) = G(Iq1,iq2) + Enn
        ENDIF
        kmat = Imat
        Imat = Imat - 1
      ENDIF
      END

      SUBROUTINE GAUSS
!***********************************************************************
! SOLVE ANY LINEAR SET OF UP TO MAXMAT EQUATIONS
! NUMBER OF EQUATIONS = IMAT
!***********************************************************************
      use FLINT_CEA_data
      IMPLICIT NONE
! LOCAL VARIABLES
      INTEGER i,imatp1,j,k,nn,nnp1
      REAL*8 bigno,coefx(50),tmp
      REAL*8 DABS,DMAX1
      SAVE coefx,i,imatp1,j,k,nn,nnp1,tmp

      DATA bigno/1.E+25/
! BEGIN ELIMINATION OF NNTH VARIABLE
      imatp1 = Imat + 1
      DO nn = 1,Imat
        IF ( nn.NE.Imat ) THEN
! SEARCH FOR MAXIMUM COEFFICIENT IN EACH ROW
          nnp1 = nn + 1
          DO i = nn,Imat
            coefx(i) = bigno
            IF ( G(i,nn).NE.0. ) THEN
              coefx(i) = 0.
              DO j = nnp1,imatp1
                coefx(i) = DMAX1(coefx(i),DABS(G(i,j)))
              ENDDO
              tmp = DABS(G(i,nn))
              IF ( bigno*tmp.GT.coefx(i) ) THEN
                coefx(i) = coefx(i)/tmp
              ELSE
                coefx(i) = bigno
              ENDIF
            ENDIF
          ENDDO
! LOCATE ROW WITH SMALLEST MAXIMUM COEFFICIENT
          tmp = bigno
          i = 0
          DO j = nn,Imat
            IF ( coefx(j).LT.tmp ) THEN
              tmp = coefx(j)
              i = j
            ENDIF
          ENDDO
          IF ( i.EQ.0 ) THEN
            Msing = nn
            GOTO 99999
! INDEX I LOCATES EQUATION TO BE USED FOR ELIMINATING THE NTH
! VARIABLE FROM THE REMAINING EQUATIONS
! INTERCHANGE EQUATIONS I AND NN
          ELSEIF ( nn.NE.i ) THEN
            DO j = nn,imatp1
              tmp = G(i,j)
              G(i,j) = G(nn,j)
              G(nn,j) = tmp
            ENDDO
          ENDIF
        ELSEIF ( G(nn,nn).EQ.0 ) THEN
          Msing = nn
          GOTO 99999
        ENDIF
! DIVIDE NTH ROW BY NTH DIAGONAL ELEMENT AND ELIMINATE THE NTH
! VARIABLE FROM THE REMAINING EQUATIONS
        k = nn + 1
        tmp = G(nn,nn)
        IF ( tmp.EQ.0. ) THEN
          Msing = nn
          GOTO 99999
        ELSE
          DO j = k,imatp1
            G(nn,j) = G(nn,j)/tmp
          ENDDO
          IF ( k.NE.imatp1 ) THEN
            DO i = k,Imat
              DO j = k,imatp1
                G(i,j) = G(i,j) - G(i,nn)*G(nn,j)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
! BACKSOLVE FOR THE VARIABLES
      k = Imat
 100  j = k + 1
      X(k) = 0.0D0
      tmp = 0.0
      IF ( Imat.GE.j ) THEN
        DO i = j,Imat
          tmp = tmp + G(k,i)*X(i)
        ENDDO
      ENDIF
      X(k) = G(k,imatp1) - tmp
      k = k - 1
      IF ( k.NE.0 ) GOTO 100
99999 END

      END subroutine CEA_solve

end module FLINT_CEA_solver