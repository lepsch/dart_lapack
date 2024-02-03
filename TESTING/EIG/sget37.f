      SUBROUTINE SGET37( RMAX, LMAX, NINFO, KNT, NIN )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KNT, NIN;
      // ..
      // .. Array Arguments ..
      int                LMAX( 3 ), NINFO( 3 );
      REAL               RMAX( 3 )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 ;
      REAL               EPSIN
      const              EPSIN = 5.9605E-8 ;
      int                LDT, LWORK;
      const              LDT = 20, LWORK = 2*LDT*( 10+LDT ) ;
      // ..
      // .. Local Scalars ..
      int                I, ICMP, IFND, INFO, ISCL, J, KMIN, M, N;
      REAL               BIGNUM, EPS, SMLNUM, TNRM, TOL, TOLIN, V, VIMIN, VMAX, VMUL, VRMIN
      // ..
      // .. Local Arrays ..
      bool               SELECT( LDT );
      int                IWORK( 2*LDT ), LCMP( 3 );
      REAL               DUM( 1 ), LE( LDT, LDT ), RE( LDT, LDT ), S( LDT ), SEP( LDT ), SEPIN( LDT ), SEPTMP( LDT ), SIN( LDT ), STMP( LDT ), T( LDT, LDT ), TMP( LDT, LDT ), VAL( 3 ), WI( LDT ), WIIN( LDT ), WITMP( LDT ), WORK( LWORK ), WR( LDT ), WRIN( LDT ), WRTMP( LDT )
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEHRD, SHSEQR, SLACPY, SSCAL, STREVC, STRSNA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL, SQRT
      // ..
      // .. Executable Statements ..

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
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
      // lexicographically (increasing by real part, then decreasing by
      // imaginary part)

      } // 10
      READ( NIN, FMT = * )N
      if (N == 0) RETURN;
      for (I = 1; I <= N; I++) { // 20
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
      } // 20
      for (I = 1; I <= N; I++) { // 30
         READ( NIN, FMT = * )WRIN( I ), WIIN( I ), SIN( I ), SEPIN( I )
      } // 30
      TNRM = SLANGE( 'M', N, N, TMP, LDT, WORK )

      // Begin test

      for (ISCL = 1; ISCL <= 3; ISCL++) { // 240

         // Scale input matrix

         KNT = KNT + 1
         slacpy('F', N, N, TMP, LDT, T, LDT );
         VMUL = VAL( ISCL )
         for (I = 1; I <= N; I++) { // 40
            sscal(N, VMUL, T( 1, I ), 1 );
         } // 40
         if (TNRM == ZERO) VMUL = ONE;

         // Compute eigenvalues and eigenvectors

         sgehrd(N, 1, N, T, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, INFO );
         if ( INFO != 0 ) {
            LMAX( 1 ) = KNT
            NINFO( 1 ) = NINFO( 1 ) + 1
            GO TO 240
         }
         for (J = 1; J <= N - 2; J++) { // 60
            for (I = J + 2; I <= N; I++) { // 50
               T( I, J ) = ZERO
            } // 50
         } // 60

         // Compute Schur form

         shseqr('S', 'N', N, 1, N, T, LDT, WR, WI, DUM, 1, WORK, LWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 2 ) = KNT
            NINFO( 2 ) = NINFO( 2 ) + 1
            GO TO 240
         }

         // Compute eigenvectors

         strevc('Both', 'All', SELECT, N, T, LDT, LE, LDT, RE, LDT, N, M, WORK, INFO );

         // Compute condition numbers

         strsna('Both', 'All', SELECT, N, T, LDT, LE, LDT, RE, LDT, S, SEP, N, M, WORK, N, IWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }

         // Sort eigenvalues and condition numbers lexicographically
         // to compare with inputs

         scopy(N, WR, 1, WRTMP, 1 );
         scopy(N, WI, 1, WITMP, 1 );
         scopy(N, S, 1, STMP, 1 );
         scopy(N, SEP, 1, SEPTMP, 1 );
         sscal(N, ONE / VMUL, SEPTMP, 1 );
         for (I = 1; I <= N - 1; I++) { // 80
            KMIN = I
            VRMIN = WRTMP( I )
            VIMIN = WITMP( I )
            for (J = I + 1; J <= N; J++) { // 70
               if ( WRTMP( J ) < VRMIN ) {
                  KMIN = J
                  VRMIN = WRTMP( J )
                  VIMIN = WITMP( J )
               }
            } // 70
            WRTMP( KMIN ) = WRTMP( I )
            WITMP( KMIN ) = WITMP( I )
            WRTMP( I ) = VRMIN
            WITMP( I ) = VIMIN
            VRMIN = STMP( KMIN )
            STMP( KMIN ) = STMP( I )
            STMP( I ) = VRMIN
            VRMIN = SEPTMP( KMIN )
            SEPTMP( KMIN ) = SEPTMP( I )
            SEPTMP( I ) = VRMIN
         } // 80

         // Compare condition numbers for eigenvalues
         // taking their condition numbers into account

         V = MAX( TWO*REAL( N )*EPS*TNRM, SMLNUM )
         if (TNRM == ZERO) V = ONE;
         for (I = 1; I <= N; I++) { // 90
            if ( V > SEPTMP( I ) ) {
               TOL = ONE
            } else {
               TOL = V / SEPTMP( I )
            }
            if ( V > SEPIN( I ) ) {
               TOLIN = ONE
            } else {
               TOLIN = V / SEPIN( I )
            }
            TOL = MAX( TOL, SMLNUM / EPS )
            TOLIN = MAX( TOLIN, SMLNUM / EPS )
            if ( EPS*( SIN( I )-TOLIN ) > STMP( I )+TOL ) {
               VMAX = ONE / EPS
            } else if ( SIN( I )-TOLIN > STMP( I )+TOL ) {
               VMAX = ( SIN( I )-TOLIN ) / ( STMP( I )+TOL )
            } else if ( SIN( I )+TOLIN < EPS*( STMP( I )-TOL ) ) {
               VMAX = ONE / EPS
            } else if ( SIN( I )+TOLIN < STMP( I )-TOL ) {
               VMAX = ( STMP( I )-TOL ) / ( SIN( I )+TOLIN )
            } else {
               VMAX = ONE
            }
            if ( VMAX > RMAX( 2 ) ) {
               RMAX( 2 ) = VMAX
               IF( NINFO( 2 ) == 0 ) LMAX( 2 ) = KNT
            }
         } // 90

         // Compare condition numbers for eigenvectors
         // taking their condition numbers into account

         for (I = 1; I <= N; I++) { // 100
            if ( V > SEPTMP( I )*STMP( I ) ) {
               TOL = SEPTMP( I )
            } else {
               TOL = V / STMP( I )
            }
            if ( V > SEPIN( I )*SIN( I ) ) {
               TOLIN = SEPIN( I )
            } else {
               TOLIN = V / SIN( I )
            }
            TOL = MAX( TOL, SMLNUM / EPS )
            TOLIN = MAX( TOLIN, SMLNUM / EPS )
            if ( EPS*( SEPIN( I )-TOLIN ) > SEPTMP( I )+TOL ) {
               VMAX = ONE / EPS
            } else if ( SEPIN( I )-TOLIN > SEPTMP( I )+TOL ) {
               VMAX = ( SEPIN( I )-TOLIN ) / ( SEPTMP( I )+TOL )
            } else if ( SEPIN( I )+TOLIN < EPS*( SEPTMP( I )-TOL ) ) {
               VMAX = ONE / EPS
            } else if ( SEPIN( I )+TOLIN < SEPTMP( I )-TOL ) {
               VMAX = ( SEPTMP( I )-TOL ) / ( SEPIN( I )+TOLIN )
            } else {
               VMAX = ONE
            }
            if ( VMAX > RMAX( 2 ) ) {
               RMAX( 2 ) = VMAX
               IF( NINFO( 2 ) == 0 ) LMAX( 2 ) = KNT
            }
         } // 100

         // Compare condition numbers for eigenvalues
         // without taking their condition numbers into account

         for (I = 1; I <= N; I++) { // 110
            if ( SIN( I ) <= REAL( 2*N )*EPS && STMP( I ) <= REAL( 2*N )*EPS ) {
               VMAX = ONE
            } else if ( EPS*SIN( I ) > STMP( I ) ) {
               VMAX = ONE / EPS
            } else if ( SIN( I ) > STMP( I ) ) {
               VMAX = SIN( I ) / STMP( I )
            } else if ( SIN( I ) < EPS*STMP( I ) ) {
               VMAX = ONE / EPS
            } else if ( SIN( I ) < STMP( I ) ) {
               VMAX = STMP( I ) / SIN( I )
            } else {
               VMAX = ONE
            }
            if ( VMAX > RMAX( 3 ) ) {
               RMAX( 3 ) = VMAX
               IF( NINFO( 3 ) == 0 ) LMAX( 3 ) = KNT
            }
         } // 110

         // Compare condition numbers for eigenvectors
         // without taking their condition numbers into account

         for (I = 1; I <= N; I++) { // 120
            if ( SEPIN( I ) <= V && SEPTMP( I ) <= V ) {
               VMAX = ONE
            } else if ( EPS*SEPIN( I ) > SEPTMP( I ) ) {
               VMAX = ONE / EPS
            } else if ( SEPIN( I ) > SEPTMP( I ) ) {
               VMAX = SEPIN( I ) / SEPTMP( I )
            } else if ( SEPIN( I ) < EPS*SEPTMP( I ) ) {
               VMAX = ONE / EPS
            } else if ( SEPIN( I ) < SEPTMP( I ) ) {
               VMAX = SEPTMP( I ) / SEPIN( I )
            } else {
               VMAX = ONE
            }
            if ( VMAX > RMAX( 3 ) ) {
               RMAX( 3 ) = VMAX
               IF( NINFO( 3 ) == 0 ) LMAX( 3 ) = KNT
            }
         } // 120

         // Compute eigenvalue condition numbers only and compare

         VMAX = ZERO
         DUM( 1 ) = -ONE
         scopy(N, DUM, 0, STMP, 1 );
         scopy(N, DUM, 0, SEPTMP, 1 );
         strsna('Eigcond', 'All', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         for (I = 1; I <= N; I++) { // 130
            IF( STMP( I ) != S( I ) ) VMAX = ONE / EPS             IF( SEPTMP( I ) != DUM( 1 ) ) VMAX = ONE / EPS
         } // 130

         // Compute eigenvector condition numbers only and compare

         scopy(N, DUM, 0, STMP, 1 );
         scopy(N, DUM, 0, SEPTMP, 1 );
         strsna('Veccond', 'All', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         for (I = 1; I <= N; I++) { // 140
            IF( STMP( I ) != DUM( 1 ) ) VMAX = ONE / EPS             IF( SEPTMP( I ) != SEP( I ) ) VMAX = ONE / EPS
         } // 140

         // Compute all condition numbers using SELECT and compare

         for (I = 1; I <= N; I++) { // 150
            SELECT( I ) = true;
         } // 150
         scopy(N, DUM, 0, STMP, 1 );
         scopy(N, DUM, 0, SEPTMP, 1 );
         strsna('Bothcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         for (I = 1; I <= N; I++) { // 160
            IF( SEPTMP( I ) != SEP( I ) ) VMAX = ONE / EPS             IF( STMP( I ) != S( I ) ) VMAX = ONE / EPS
         } // 160

         // Compute eigenvalue condition numbers using SELECT and compare

         scopy(N, DUM, 0, STMP, 1 );
         scopy(N, DUM, 0, SEPTMP, 1 );
         strsna('Eigcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         for (I = 1; I <= N; I++) { // 170
            IF( STMP( I ) != S( I ) ) VMAX = ONE / EPS             IF( SEPTMP( I ) != DUM( 1 ) ) VMAX = ONE / EPS
         } // 170

         // Compute eigenvector condition numbers using SELECT and compare

         scopy(N, DUM, 0, STMP, 1 );
         scopy(N, DUM, 0, SEPTMP, 1 );
         strsna('Veccond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         for (I = 1; I <= N; I++) { // 180
            IF( STMP( I ) != DUM( 1 ) ) VMAX = ONE / EPS             IF( SEPTMP( I ) != SEP( I ) ) VMAX = ONE / EPS
         } // 180
         if ( VMAX > RMAX( 1 ) ) {
            RMAX( 1 ) = VMAX
            IF( NINFO( 1 ) == 0 ) LMAX( 1 ) = KNT
         }

         // Select first real and first complex eigenvalue

         if ( WI( 1 ) == ZERO ) {
            LCMP( 1 ) = 1
            IFND = 0
            for (I = 2; I <= N; I++) { // 190
               if ( IFND == 1 || WI( I ) == ZERO ) {
                  SELECT( I ) = false;
               } else {
                  IFND = 1
                  LCMP( 2 ) = I
                  LCMP( 3 ) = I + 1
                  scopy(N, RE( 1, I ), 1, RE( 1, 2 ), 1 );
                  scopy(N, RE( 1, I+1 ), 1, RE( 1, 3 ), 1 );
                  scopy(N, LE( 1, I ), 1, LE( 1, 2 ), 1 );
                  scopy(N, LE( 1, I+1 ), 1, LE( 1, 3 ), 1 );
               }
            } // 190
            if ( IFND == 0 ) {
               ICMP = 1
            } else {
               ICMP = 3
            }
         } else {
            LCMP( 1 ) = 1
            LCMP( 2 ) = 2
            IFND = 0
            for (I = 3; I <= N; I++) { // 200
               if ( IFND == 1 || WI( I ) != ZERO ) {
                  SELECT( I ) = false;
               } else {
                  LCMP( 3 ) = I
                  IFND = 1
                  scopy(N, RE( 1, I ), 1, RE( 1, 3 ), 1 );
                  scopy(N, LE( 1, I ), 1, LE( 1, 3 ), 1 );
               }
            } // 200
            if ( IFND == 0 ) {
               ICMP = 2
            } else {
               ICMP = 3
            }
         }

         // Compute all selected condition numbers

         scopy(ICMP, DUM, 0, STMP, 1 );
         scopy(ICMP, DUM, 0, SEPTMP, 1 );
         strsna('Bothcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         for (I = 1; I <= ICMP; I++) { // 210
            J = LCMP( I )
            IF( SEPTMP( I ) != SEP( J ) ) VMAX = ONE / EPS             IF( STMP( I ) != S( J ) ) VMAX = ONE / EPS
         } // 210

         // Compute selected eigenvalue condition numbers

         scopy(ICMP, DUM, 0, STMP, 1 );
         scopy(ICMP, DUM, 0, SEPTMP, 1 );
         strsna('Eigcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         for (I = 1; I <= ICMP; I++) { // 220
            J = LCMP( I )
            IF( STMP( I ) != S( J ) ) VMAX = ONE / EPS             IF( SEPTMP( I ) != DUM( 1 ) ) VMAX = ONE / EPS
         } // 220

         // Compute selected eigenvector condition numbers

         scopy(ICMP, DUM, 0, STMP, 1 );
         scopy(ICMP, DUM, 0, SEPTMP, 1 );
         strsna('Veccond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO );
         if ( INFO != 0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         for (I = 1; I <= ICMP; I++) { // 230
            J = LCMP( I )
            IF( STMP( I ) != DUM( 1 ) ) VMAX = ONE / EPS             IF( SEPTMP( I ) != SEP( J ) ) VMAX = ONE / EPS
         } // 230
         if ( VMAX > RMAX( 1 ) ) {
            RMAX( 1 ) = VMAX
            IF( NINFO( 1 ) == 0 ) LMAX( 1 ) = KNT
         }
      } // 240
      GO TO 10

      // End of SGET37

      }
