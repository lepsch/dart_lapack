      SUBROUTINE ZGET37( RMAX, LMAX, NINFO, KNT, NIN )

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
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 ;
      double             EPSIN;
      const              EPSIN = 5.9605D-8 ;
      int                LDT, LWORK;
      const              LDT = 20, LWORK = 2*LDT*( 10+LDT ) ;
      // ..
      // .. Local Scalars ..
      int                I, ICMP, INFO, ISCL, ISRT, J, KMIN, M, N;
      double             BIGNUM, EPS, SMLNUM, TNRM, TOL, TOLIN, V, VCMIN, VMAX, VMIN, VMUL;
      // ..
      // .. Local Arrays ..
      bool               SELECT( LDT );
      int                LCMP( 3 );
      double             DUM( 1 ), RWORK( 2*LDT ), S( LDT ), SEP( LDT ), SEPIN( LDT ), SEPTMP( LDT ), SIN( LDT ), STMP( LDT ), VAL( 3 ), WIIN( LDT ), WRIN( LDT ), WSRT( LDT );
      COMPLEX*16         CDUM( 1 ), LE( LDT, LDT ), RE( LDT, LDT ), T( LDT, LDT ), TMP( LDT, LDT ), W( LDT ), WORK( LWORK ), WTMP( LDT )
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSCAL, ZCOPY, ZDSCAL, ZGEHRD, ZHSEQR, ZLACPY, ZTREVC, ZTRSNA
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
      VAL( 3 ) = SQRT( BIGNUM )

      // Read input data until N=0.  Assume input eigenvalues are sorted
      // lexicographically (increasing by real part if ISRT = 0,
      // increasing by imaginary part if ISRT = 1)

      } // 10
      READ( NIN, FMT = * )N, ISRT
      IF( N.EQ.0 ) RETURN
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
      } // 20
      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )WRIN( I ), WIIN( I ), SIN( I ), SEPIN( I )
      } // 30
      TNRM = ZLANGE( 'M', N, N, TMP, LDT, RWORK )
      for (ISCL = 1; ISCL <= 3; ISCL++) { // 260

         // Scale input matrix

         KNT = KNT + 1
         zlacpy('F', N, N, TMP, LDT, T, LDT );
         VMUL = VAL( ISCL )
         for (I = 1; I <= N; I++) { // 40
            zdscal(N, VMUL, T( 1, I ), 1 );
         } // 40
         IF( TNRM.EQ.ZERO ) VMUL = ONE

         // Compute eigenvalues and eigenvectors

         zgehrd(N, 1, N, T, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 1 ) = KNT
            NINFO( 1 ) = NINFO( 1 ) + 1
            GO TO 260
         }
         DO 60 J = 1, N - 2
            DO 50 I = J + 2, N
               T( I, J ) = ZERO
            } // 50
         } // 60

         // Compute Schur form

         zhseqr('S', 'N', N, 1, N, T, LDT, W, CDUM, 1, WORK, LWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 2 ) = KNT
            NINFO( 2 ) = NINFO( 2 ) + 1
            GO TO 260
         }

         // Compute eigenvectors

         for (I = 1; I <= N; I++) { // 70
            SELECT( I ) = .TRUE.
         } // 70
         ztrevc('B', 'A', SELECT, N, T, LDT, LE, LDT, RE, LDT, N, M, WORK, RWORK, INFO );

         // Compute condition numbers

         ztrsna('B', 'A', SELECT, N, T, LDT, LE, LDT, RE, LDT, S, SEP, N, M, WORK, N, RWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 260
         }

         // Sort eigenvalues and condition numbers lexicographically
         // to compare with inputs

         zcopy(N, W, 1, WTMP, 1 );
         if ( ISRT.EQ.0 ) {

            // Sort by increasing real part

            for (I = 1; I <= N; I++) { // 80
               WSRT( I ) = DBLE( W( I ) )
            } // 80
         } else {

            // Sort by increasing imaginary part

            for (I = 1; I <= N; I++) { // 90
               WSRT( I ) = DIMAG( W( I ) )
            } // 90
         }
         dcopy(N, S, 1, STMP, 1 );
         dcopy(N, SEP, 1, SEPTMP, 1 );
         dscal(N, ONE / VMUL, SEPTMP, 1 );
         DO 110 I = 1, N - 1
            KMIN = I
            VMIN = WSRT( I )
            DO 100 J = I + 1, N
               if ( WSRT( J ).LT.VMIN ) {
                  KMIN = J
                  VMIN = WSRT( J )
               }
            } // 100
            WSRT( KMIN ) = WSRT( I )
            WSRT( I ) = VMIN
            VCMIN = DBLE( WTMP( I ) )
            WTMP( I ) = W( KMIN )
            WTMP( KMIN ) = VCMIN
            VMIN = STMP( KMIN )
            STMP( KMIN ) = STMP( I )
            STMP( I ) = VMIN
            VMIN = SEPTMP( KMIN )
            SEPTMP( KMIN ) = SEPTMP( I )
            SEPTMP( I ) = VMIN
         } // 110

         // Compare condition numbers for eigenvalues
         // taking their condition numbers into account

         V = MAX( TWO*DBLE( N )*EPS*TNRM, SMLNUM )
         IF( TNRM.EQ.ZERO ) V = ONE
         for (I = 1; I <= N; I++) { // 120
            if ( V.GT.SEPTMP( I ) ) {
               TOL = ONE
            } else {
               TOL = V / SEPTMP( I )
            }
            if ( V.GT.SEPIN( I ) ) {
               TOLIN = ONE
            } else {
               TOLIN = V / SEPIN( I )
            }
            TOL = MAX( TOL, SMLNUM / EPS )
            TOLIN = MAX( TOLIN, SMLNUM / EPS )
            if ( EPS*( SIN( I )-TOLIN ).GT.STMP( I )+TOL ) {
               VMAX = ONE / EPS
            } else if ( SIN( I )-TOLIN.GT.STMP( I )+TOL ) {
               VMAX = ( SIN( I )-TOLIN ) / ( STMP( I )+TOL )
            } else if ( SIN( I )+TOLIN.LT.EPS*( STMP( I )-TOL ) ) {
               VMAX = ONE / EPS
            } else if ( SIN( I )+TOLIN.LT.STMP( I )-TOL ) {
               VMAX = ( STMP( I )-TOL ) / ( SIN( I )+TOLIN )
            } else {
               VMAX = ONE
            }
            if ( VMAX.GT.RMAX( 2 ) ) {
               RMAX( 2 ) = VMAX
               IF( NINFO( 2 ).EQ.0 ) LMAX( 2 ) = KNT
            }
         } // 120

         // Compare condition numbers for eigenvectors
         // taking their condition numbers into account

         for (I = 1; I <= N; I++) { // 130
            if ( V.GT.SEPTMP( I )*STMP( I ) ) {
               TOL = SEPTMP( I )
            } else {
               TOL = V / STMP( I )
            }
            if ( V.GT.SEPIN( I )*SIN( I ) ) {
               TOLIN = SEPIN( I )
            } else {
               TOLIN = V / SIN( I )
            }
            TOL = MAX( TOL, SMLNUM / EPS )
            TOLIN = MAX( TOLIN, SMLNUM / EPS )
            if ( EPS*( SEPIN( I )-TOLIN ).GT.SEPTMP( I )+TOL ) {
               VMAX = ONE / EPS
            } else if ( SEPIN( I )-TOLIN.GT.SEPTMP( I )+TOL ) {
               VMAX = ( SEPIN( I )-TOLIN ) / ( SEPTMP( I )+TOL )
            } else if ( SEPIN( I )+TOLIN.LT.EPS*( SEPTMP( I )-TOL ) ) {
               VMAX = ONE / EPS
            } else if ( SEPIN( I )+TOLIN.LT.SEPTMP( I )-TOL ) {
               VMAX = ( SEPTMP( I )-TOL ) / ( SEPIN( I )+TOLIN )
            } else {
               VMAX = ONE
            }
            if ( VMAX.GT.RMAX( 2 ) ) {
               RMAX( 2 ) = VMAX
               IF( NINFO( 2 ).EQ.0 ) LMAX( 2 ) = KNT
            }
         } // 130

         // Compare condition numbers for eigenvalues
         // without taking their condition numbers into account

         for (I = 1; I <= N; I++) { // 140
            if ( SIN( I ).LE.DBLE( 2*N )*EPS .AND. STMP( I ).LE. DBLE( 2*N )*EPS ) {
               VMAX = ONE
            } else if ( EPS*SIN( I ).GT.STMP( I ) ) {
               VMAX = ONE / EPS
            } else if ( SIN( I ).GT.STMP( I ) ) {
               VMAX = SIN( I ) / STMP( I )
            } else if ( SIN( I ).LT.EPS*STMP( I ) ) {
               VMAX = ONE / EPS
            } else if ( SIN( I ).LT.STMP( I ) ) {
               VMAX = STMP( I ) / SIN( I )
            } else {
               VMAX = ONE
            }
            if ( VMAX.GT.RMAX( 3 ) ) {
               RMAX( 3 ) = VMAX
               IF( NINFO( 3 ).EQ.0 ) LMAX( 3 ) = KNT
            }
         } // 140

         // Compare condition numbers for eigenvectors
         // without taking their condition numbers into account

         for (I = 1; I <= N; I++) { // 150
            if ( SEPIN( I ).LE.V .AND. SEPTMP( I ).LE.V ) {
               VMAX = ONE
            } else if ( EPS*SEPIN( I ).GT.SEPTMP( I ) ) {
               VMAX = ONE / EPS
            } else if ( SEPIN( I ).GT.SEPTMP( I ) ) {
               VMAX = SEPIN( I ) / SEPTMP( I )
            } else if ( SEPIN( I ).LT.EPS*SEPTMP( I ) ) {
               VMAX = ONE / EPS
            } else if ( SEPIN( I ).LT.SEPTMP( I ) ) {
               VMAX = SEPTMP( I ) / SEPIN( I )
            } else {
               VMAX = ONE
            }
            if ( VMAX.GT.RMAX( 3 ) ) {
               RMAX( 3 ) = VMAX
               IF( NINFO( 3 ).EQ.0 ) LMAX( 3 ) = KNT
            }
         } // 150

         // Compute eigenvalue condition numbers only and compare

         VMAX = ZERO
         DUM( 1 ) = -ONE
         dcopy(N, DUM, 0, STMP, 1 );
         dcopy(N, DUM, 0, SEPTMP, 1 );
         ztrsna('E', 'A', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, RWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 260
         }
         for (I = 1; I <= N; I++) { // 160
            IF( STMP( I ).NE.S( I ) ) VMAX = ONE / EPS             IF( SEPTMP( I ).NE.DUM( 1 ) ) VMAX = ONE / EPS
         } // 160

         // Compute eigenvector condition numbers only and compare

         dcopy(N, DUM, 0, STMP, 1 );
         dcopy(N, DUM, 0, SEPTMP, 1 );
         ztrsna('V', 'A', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, RWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 260
         }
         for (I = 1; I <= N; I++) { // 170
            IF( STMP( I ).NE.DUM( 1 ) ) VMAX = ONE / EPS             IF( SEPTMP( I ).NE.SEP( I ) ) VMAX = ONE / EPS
         } // 170

         // Compute all condition numbers using SELECT and compare

         for (I = 1; I <= N; I++) { // 180
            SELECT( I ) = .TRUE.
         } // 180
         dcopy(N, DUM, 0, STMP, 1 );
         dcopy(N, DUM, 0, SEPTMP, 1 );
         ztrsna('B', 'S', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, RWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 260
         }
         for (I = 1; I <= N; I++) { // 190
            IF( SEPTMP( I ).NE.SEP( I ) ) VMAX = ONE / EPS             IF( STMP( I ).NE.S( I ) ) VMAX = ONE / EPS
         } // 190

         // Compute eigenvalue condition numbers using SELECT and compare

         dcopy(N, DUM, 0, STMP, 1 );
         dcopy(N, DUM, 0, SEPTMP, 1 );
         ztrsna('E', 'S', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, RWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 260
         }
         for (I = 1; I <= N; I++) { // 200
            IF( STMP( I ).NE.S( I ) ) VMAX = ONE / EPS             IF( SEPTMP( I ).NE.DUM( 1 ) ) VMAX = ONE / EPS
         } // 200

         // Compute eigenvector condition numbers using SELECT and compare

         dcopy(N, DUM, 0, STMP, 1 );
         dcopy(N, DUM, 0, SEPTMP, 1 );
         ztrsna('V', 'S', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, RWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 260
         }
         for (I = 1; I <= N; I++) { // 210
            IF( STMP( I ).NE.DUM( 1 ) ) VMAX = ONE / EPS             IF( SEPTMP( I ).NE.SEP( I ) ) VMAX = ONE / EPS
         } // 210
         if ( VMAX.GT.RMAX( 1 ) ) {
            RMAX( 1 ) = VMAX
            IF( NINFO( 1 ).EQ.0 ) LMAX( 1 ) = KNT
         }

         // Select second and next to last eigenvalues

         for (I = 1; I <= N; I++) { // 220
            SELECT( I ) = .FALSE.
         } // 220
         ICMP = 0
         if ( N.GT.1 ) {
            ICMP = 1
            LCMP( 1 ) = 2
            SELECT( 2 ) = .TRUE.
            zcopy(N, RE( 1, 2 ), 1, RE( 1, 1 ), 1 );
            zcopy(N, LE( 1, 2 ), 1, LE( 1, 1 ), 1 );
         }
         if ( N.GT.3 ) {
            ICMP = 2
            LCMP( 2 ) = N - 1
            SELECT( N-1 ) = .TRUE.
            zcopy(N, RE( 1, N-1 ), 1, RE( 1, 2 ), 1 );
            zcopy(N, LE( 1, N-1 ), 1, LE( 1, 2 ), 1 );
         }

         // Compute all selected condition numbers

         dcopy(ICMP, DUM, 0, STMP, 1 );
         dcopy(ICMP, DUM, 0, SEPTMP, 1 );
         ztrsna('B', 'S', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, RWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 260
         }
         for (I = 1; I <= ICMP; I++) { // 230
            J = LCMP( I )
            IF( SEPTMP( I ).NE.SEP( J ) ) VMAX = ONE / EPS             IF( STMP( I ).NE.S( J ) ) VMAX = ONE / EPS
         } // 230

         // Compute selected eigenvalue condition numbers

         dcopy(ICMP, DUM, 0, STMP, 1 );
         dcopy(ICMP, DUM, 0, SEPTMP, 1 );
         ztrsna('E', 'S', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, RWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 260
         }
         for (I = 1; I <= ICMP; I++) { // 240
            J = LCMP( I )
            IF( STMP( I ).NE.S( J ) ) VMAX = ONE / EPS             IF( SEPTMP( I ).NE.DUM( 1 ) ) VMAX = ONE / EPS
         } // 240

         // Compute selected eigenvector condition numbers

         dcopy(ICMP, DUM, 0, STMP, 1 );
         dcopy(ICMP, DUM, 0, SEPTMP, 1 );
         ztrsna('V', 'S', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, RWORK, INFO );
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 260
         }
         for (I = 1; I <= ICMP; I++) { // 250
            J = LCMP( I )
            IF( STMP( I ).NE.DUM( 1 ) ) VMAX = ONE / EPS             IF( SEPTMP( I ).NE.SEP( J ) ) VMAX = ONE / EPS
         } // 250
         if ( VMAX.GT.RMAX( 1 ) ) {
            RMAX( 1 ) = VMAX
            IF( NINFO( 1 ).EQ.0 ) LMAX( 1 ) = KNT
         }
      } // 260
      GO TO 10

      // End of ZGET37

      }
