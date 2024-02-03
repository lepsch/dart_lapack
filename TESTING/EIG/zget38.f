      SUBROUTINE ZGET38( RMAX, LMAX, NINFO, KNT, NIN )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, NIN;
      // ..
      // .. Array Arguments ..
      int                LMAX( 3 ), NINFO( 3 );
      double             RMAX( 3 );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                LDT, LWORK;
      const              LDT = 20, LWORK = 2*LDT*( 10+LDT ) ;
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 ;
      double             EPSIN;
      const              EPSIN = 5.9605D-8 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, ISCL, ISRT, ITMP, J, KMIN, M, N, NDIM;
      double             BIGNUM, EPS, S, SEP, SEPIN, SEPTMP, SIN, SMLNUM, STMP, TNRM, TOL, TOLIN, V, VMAX, VMIN, VMUL;
      // ..
      // .. Local Arrays ..
      bool               SELECT( LDT );
      int                IPNT( LDT ), ISELEC( LDT );
      double             RESULT( 2 ), RWORK( LDT ), VAL( 3 ), WSRT( LDT );
      COMPLEX*16         Q( LDT, LDT ), QSAV( LDT, LDT ), QTMP( LDT, LDT ), T( LDT, LDT ), TMP( LDT, LDT ), TSAV( LDT, LDT ), TSAV1( LDT, LDT ), TTMP( LDT, LDT ), W( LDT ), WORK( LWORK ), WTMP( LDT );
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL, ZGEHRD, ZHSEQR, ZHST01, ZLACPY, ZTRSEN, ZUNGHR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DIMAG, MAX, SQRT
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      // EPSIN = 2**(-24) = precision to which input data computed

      EPS = MAX( EPS, EPSIN )
      RMAX( 1 ) = ZERO
      RMAX( 2 ) = ZERO
      RMAX( 3 ) = ZERO
      LMAX( 1 ) = 0
      LMAX( 2 ) = 0
      LMAX( 3 ) = 0
      KNT = 0
      NINFO( 1 ) = 0
      NINFO( 2 ) = 0
      NINFO( 3 ) = 0
      VAL( 1 ) = SQRT( SMLNUM )
      VAL( 2 ) = ONE
      VAL( 3 ) = SQRT( SQRT( BIGNUM ) )

      // Read input data until N=0.  Assume input eigenvalues are sorted
      // lexicographically (increasing by real part, then decreasing by
      // imaginary part)

      } // 10
      READ( NIN, FMT = * )N, NDIM, ISRT
      if (N == 0) RETURN;
      READ( NIN, FMT = * )( ISELEC( I ), I = 1, NDIM )
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
      } // 20
      READ( NIN, FMT = * )SIN, SEPIN

      TNRM = ZLANGE( 'M', N, N, TMP, LDT, RWORK )
      for (ISCL = 1; ISCL <= 3; ISCL++) { // 200

         // Scale input matrix

         KNT = KNT + 1
         zlacpy('F', N, N, TMP, LDT, T, LDT );
         VMUL = VAL( ISCL )
         for (I = 1; I <= N; I++) { // 30
            zdscal(N, VMUL, T( 1, I ), 1 );
         } // 30
         if (TNRM == ZERO) VMUL = ONE;
         zlacpy('F', N, N, T, LDT, TSAV, LDT );

         // Compute Schur form

         zgehrd(N, 1, N, T, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 1 ) = KNT
            NINFO( 1 ) = NINFO( 1 ) + 1
            GO TO 200
         }

         // Generate unitary matrix

         zlacpy('L', N, N, T, LDT, Q, LDT );
         zunghr(N, 1, N, Q, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, INFO );

         // Compute Schur form

         for (J = 1; J <= N - 2; J++) { // 50
            for (I = J + 2; I <= N; I++) { // 40
               T( I, J ) = CZERO
            } // 40
         } // 50
         zhseqr('S', 'V', N, 1, N, T, LDT, W, Q, LDT, WORK, LWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 2 ) = KNT
            NINFO( 2 ) = NINFO( 2 ) + 1
            GO TO 200
         }

         // Sort, select eigenvalues

         for (I = 1; I <= N; I++) { // 60
            IPNT( I ) = I
            SELECT( I ) = false;
         } // 60
         if ( ISRT == 0 ) {
            for (I = 1; I <= N; I++) { // 70
               WSRT( I ) = DBLE( W( I ) )
            } // 70
         } else {
            for (I = 1; I <= N; I++) { // 80
               WSRT( I ) = DIMAG( W( I ) )
            } // 80
         }
         for (I = 1; I <= N - 1; I++) { // 100
            KMIN = I
            VMIN = WSRT( I )
            for (J = I + 1; J <= N; J++) { // 90
               if ( WSRT( J ).LT.VMIN ) {
                  KMIN = J
                  VMIN = WSRT( J )
               }
            } // 90
            WSRT( KMIN ) = WSRT( I )
            WSRT( I ) = VMIN
            ITMP = IPNT( I )
            IPNT( I ) = IPNT( KMIN )
            IPNT( KMIN ) = ITMP
         } // 100
         for (I = 1; I <= NDIM; I++) { // 110
            SELECT( IPNT( ISELEC( I ) ) ) = true;
         } // 110

         // Compute condition numbers

         zlacpy('F', N, N, Q, LDT, QSAV, LDT );
         zlacpy('F', N, N, T, LDT, TSAV1, LDT );
         ztrsen('B', 'V', SELECT, N, T, LDT, Q, LDT, WTMP, M, S, SEP, WORK, LWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 200
         }
         SEPTMP = SEP / VMUL
         STMP = S

         // Compute residuals

         zhst01(N, 1, N, TSAV, LDT, T, LDT, Q, LDT, WORK, LWORK, RWORK, RESULT );
         VMAX = MAX( RESULT( 1 ), RESULT( 2 ) )
         if ( VMAX.GT.RMAX( 1 ) ) {
            RMAX( 1 ) = VMAX
            IF( NINFO( 1 ) == 0 ) LMAX( 1 ) = KNT
         }

         // Compare condition number for eigenvalue cluster
         // taking its condition number into account

         V = MAX( TWO*DBLE( N )*EPS*TNRM, SMLNUM )
         if (TNRM == ZERO) V = ONE;
         if ( V.GT.SEPTMP ) {
            TOL = ONE
         } else {
            TOL = V / SEPTMP
         }
         if ( V.GT.SEPIN ) {
            TOLIN = ONE
         } else {
            TOLIN = V / SEPIN
         }
         TOL = MAX( TOL, SMLNUM / EPS )
         TOLIN = MAX( TOLIN, SMLNUM / EPS )
         if ( EPS*( SIN-TOLIN ).GT.STMP+TOL ) {
            VMAX = ONE / EPS
         } else if ( SIN-TOLIN.GT.STMP+TOL ) {
            VMAX = ( SIN-TOLIN ) / ( STMP+TOL )
         } else if ( SIN+TOLIN.LT.EPS*( STMP-TOL ) ) {
            VMAX = ONE / EPS
         } else if ( SIN+TOLIN.LT.STMP-TOL ) {
            VMAX = ( STMP-TOL ) / ( SIN+TOLIN )
         } else {
            VMAX = ONE
         }
         if ( VMAX.GT.RMAX( 2 ) ) {
            RMAX( 2 ) = VMAX
            IF( NINFO( 2 ) == 0 ) LMAX( 2 ) = KNT
         }

         // Compare condition numbers for invariant subspace
         // taking its condition number into account

         if ( V.GT.SEPTMP*STMP ) {
            TOL = SEPTMP
         } else {
            TOL = V / STMP
         }
         if ( V.GT.SEPIN*SIN ) {
            TOLIN = SEPIN
         } else {
            TOLIN = V / SIN
         }
         TOL = MAX( TOL, SMLNUM / EPS )
         TOLIN = MAX( TOLIN, SMLNUM / EPS )
         if ( EPS*( SEPIN-TOLIN ).GT.SEPTMP+TOL ) {
            VMAX = ONE / EPS
         } else if ( SEPIN-TOLIN.GT.SEPTMP+TOL ) {
            VMAX = ( SEPIN-TOLIN ) / ( SEPTMP+TOL )
         } else if ( SEPIN+TOLIN.LT.EPS*( SEPTMP-TOL ) ) {
            VMAX = ONE / EPS
         } else if ( SEPIN+TOLIN.LT.SEPTMP-TOL ) {
            VMAX = ( SEPTMP-TOL ) / ( SEPIN+TOLIN )
         } else {
            VMAX = ONE
         }
         if ( VMAX.GT.RMAX( 2 ) ) {
            RMAX( 2 ) = VMAX
            IF( NINFO( 2 ) == 0 ) LMAX( 2 ) = KNT
         }

         // Compare condition number for eigenvalue cluster
         // without taking its condition number into account

         if ( SIN.LE.DBLE( 2*N )*EPS .AND. STMP.LE.DBLE( 2*N )*EPS ) {
            VMAX = ONE
         } else if ( EPS*SIN.GT.STMP ) {
            VMAX = ONE / EPS
         } else if ( SIN.GT.STMP ) {
            VMAX = SIN / STMP
         } else if ( SIN.LT.EPS*STMP ) {
            VMAX = ONE / EPS
         } else if ( SIN.LT.STMP ) {
            VMAX = STMP / SIN
         } else {
            VMAX = ONE
         }
         if ( VMAX.GT.RMAX( 3 ) ) {
            RMAX( 3 ) = VMAX
            IF( NINFO( 3 ) == 0 ) LMAX( 3 ) = KNT
         }

         // Compare condition numbers for invariant subspace
         // without taking its condition number into account

         if ( SEPIN.LE.V .AND. SEPTMP.LE.V ) {
            VMAX = ONE
         } else if ( EPS*SEPIN.GT.SEPTMP ) {
            VMAX = ONE / EPS
         } else if ( SEPIN.GT.SEPTMP ) {
            VMAX = SEPIN / SEPTMP
         } else if ( SEPIN.LT.EPS*SEPTMP ) {
            VMAX = ONE / EPS
         } else if ( SEPIN.LT.SEPTMP ) {
            VMAX = SEPTMP / SEPIN
         } else {
            VMAX = ONE
         }
         if ( VMAX.GT.RMAX( 3 ) ) {
            RMAX( 3 ) = VMAX
            IF( NINFO( 3 ) == 0 ) LMAX( 3 ) = KNT
         }

         // Compute eigenvalue condition number only and compare
         // Update Q

         VMAX = ZERO
         zlacpy('F', N, N, TSAV1, LDT, TTMP, LDT );
         zlacpy('F', N, N, QSAV, LDT, QTMP, LDT );
         SEPTMP = -ONE
         STMP = -ONE
         ztrsen('E', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WTMP, M, STMP, SEPTMP, WORK, LWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 200
         }
         if (S.NE.STMP) VMAX = ONE / EPS          IF( -ONE.NE.SEPTMP ) VMAX = ONE / EPS;
         for (I = 1; I <= N; I++) { // 130
            for (J = 1; J <= N; J++) { // 120
               IF( TTMP( I, J ).NE.T( I, J ) ) VMAX = ONE / EPS                IF( QTMP( I, J ).NE.Q( I, J ) ) VMAX = ONE / EPS
            } // 120
         } // 130

         // Compute invariant subspace condition number only and compare
         // Update Q

         zlacpy('F', N, N, TSAV1, LDT, TTMP, LDT );
         zlacpy('F', N, N, QSAV, LDT, QTMP, LDT );
         SEPTMP = -ONE
         STMP = -ONE
         ztrsen('V', 'V', SELECT, N, TTMP, LDT, QTMP, LDT, WTMP, M, STMP, SEPTMP, WORK, LWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 200
         }
         if (-ONE.NE.STMP) VMAX = ONE / EPS          IF( SEP.NE.SEPTMP ) VMAX = ONE / EPS;
         for (I = 1; I <= N; I++) { // 150
            for (J = 1; J <= N; J++) { // 140
               IF( TTMP( I, J ).NE.T( I, J ) ) VMAX = ONE / EPS                IF( QTMP( I, J ).NE.Q( I, J ) ) VMAX = ONE / EPS
            } // 140
         } // 150

         // Compute eigenvalue condition number only and compare
         // Do not update Q

         zlacpy('F', N, N, TSAV1, LDT, TTMP, LDT );
         zlacpy('F', N, N, QSAV, LDT, QTMP, LDT );
         SEPTMP = -ONE
         STMP = -ONE
         ztrsen('E', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WTMP, M, STMP, SEPTMP, WORK, LWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 200
         }
         if (S.NE.STMP) VMAX = ONE / EPS          IF( -ONE.NE.SEPTMP ) VMAX = ONE / EPS;
         for (I = 1; I <= N; I++) { // 170
            for (J = 1; J <= N; J++) { // 160
               IF( TTMP( I, J ).NE.T( I, J ) ) VMAX = ONE / EPS                IF( QTMP( I, J ).NE.QSAV( I, J ) ) VMAX = ONE / EPS
            } // 160
         } // 170

         // Compute invariant subspace condition number only and compare
         // Do not update Q

         zlacpy('F', N, N, TSAV1, LDT, TTMP, LDT );
         zlacpy('F', N, N, QSAV, LDT, QTMP, LDT );
         SEPTMP = -ONE
         STMP = -ONE
         ztrsen('V', 'N', SELECT, N, TTMP, LDT, QTMP, LDT, WTMP, M, STMP, SEPTMP, WORK, LWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 200
         }
         if (-ONE.NE.STMP) VMAX = ONE / EPS          IF( SEP.NE.SEPTMP ) VMAX = ONE / EPS;
         for (I = 1; I <= N; I++) { // 190
            for (J = 1; J <= N; J++) { // 180
               IF( TTMP( I, J ).NE.T( I, J ) ) VMAX = ONE / EPS                IF( QTMP( I, J ).NE.QSAV( I, J ) ) VMAX = ONE / EPS
            } // 180
         } // 190
         if ( VMAX.GT.RMAX( 1 ) ) {
            RMAX( 1 ) = VMAX
            IF( NINFO( 1 ) == 0 ) LMAX( 1 ) = KNT
         }
      } // 200
      GO TO 10

      // End of ZGET38

      }
