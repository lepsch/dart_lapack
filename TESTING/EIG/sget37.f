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

   10 CONTINUE
      READ( NIN, FMT = * )N
      IF( N.EQ.0 ) RETURN
      DO 20 I = 1, N
         READ( NIN, FMT = * )( TMP( I, J ), J = 1, N )
   20 CONTINUE
      DO 30 I = 1, N
         READ( NIN, FMT = * )WRIN( I ), WIIN( I ), SIN( I ), SEPIN( I )
   30 CONTINUE
      TNRM = SLANGE( 'M', N, N, TMP, LDT, WORK )

      // Begin test

      DO 240 ISCL = 1, 3

         // Scale input matrix

         KNT = KNT + 1
         CALL SLACPY( 'F', N, N, TMP, LDT, T, LDT )
         VMUL = VAL( ISCL )
         DO 40 I = 1, N
            CALL SSCAL( N, VMUL, T( 1, I ), 1 )
   40    CONTINUE
         IF( TNRM.EQ.ZERO ) VMUL = ONE

         // Compute eigenvalues and eigenvectors

         CALL SGEHRD( N, 1, N, T, LDT, WORK( 1 ), WORK( N+1 ), LWORK-N, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 1 ) = KNT
            NINFO( 1 ) = NINFO( 1 ) + 1
            GO TO 240
         }
         DO 60 J = 1, N - 2
            DO 50 I = J + 2, N
               T( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE

         // Compute Schur form

         CALL SHSEQR( 'S', 'N', N, 1, N, T, LDT, WR, WI, DUM, 1, WORK, LWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 2 ) = KNT
            NINFO( 2 ) = NINFO( 2 ) + 1
            GO TO 240
         }

         // Compute eigenvectors

         CALL STREVC( 'Both', 'All', SELECT, N, T, LDT, LE, LDT, RE, LDT, N, M, WORK, INFO )

         // Compute condition numbers

         CALL STRSNA( 'Both', 'All', SELECT, N, T, LDT, LE, LDT, RE, LDT, S, SEP, N, M, WORK, N, IWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }

         // Sort eigenvalues and condition numbers lexicographically
         // to compare with inputs

         CALL SCOPY( N, WR, 1, WRTMP, 1 )
         CALL SCOPY( N, WI, 1, WITMP, 1 )
         CALL SCOPY( N, S, 1, STMP, 1 )
         CALL SCOPY( N, SEP, 1, SEPTMP, 1 )
         CALL SSCAL( N, ONE / VMUL, SEPTMP, 1 )
         DO 80 I = 1, N - 1
            KMIN = I
            VRMIN = WRTMP( I )
            VIMIN = WITMP( I )
            DO 70 J = I + 1, N
               if ( WRTMP( J ).LT.VRMIN ) {
                  KMIN = J
                  VRMIN = WRTMP( J )
                  VIMIN = WITMP( J )
               }
   70       CONTINUE
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
   80    CONTINUE

         // Compare condition numbers for eigenvalues
         // taking their condition numbers into account

         V = MAX( TWO*REAL( N )*EPS*TNRM, SMLNUM )
         IF( TNRM.EQ.ZERO ) V = ONE
         DO 90 I = 1, N
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
   90    CONTINUE

         // Compare condition numbers for eigenvectors
         // taking their condition numbers into account

         DO 100 I = 1, N
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
  100    CONTINUE

         // Compare condition numbers for eigenvalues
         // without taking their condition numbers into account

         DO 110 I = 1, N
            if ( SIN( I ).LE.REAL( 2*N )*EPS .AND. STMP( I ).LE. REAL( 2*N )*EPS ) {
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
  110    CONTINUE

         // Compare condition numbers for eigenvectors
         // without taking their condition numbers into account

         DO 120 I = 1, N
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
  120    CONTINUE

         // Compute eigenvalue condition numbers only and compare

         VMAX = ZERO
         DUM( 1 ) = -ONE
         CALL SCOPY( N, DUM, 0, STMP, 1 )
         CALL SCOPY( N, DUM, 0, SEPTMP, 1 )
         CALL STRSNA( 'Eigcond', 'All', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         DO 130 I = 1, N
            IF( STMP( I ).NE.S( I ) ) VMAX = ONE / EPS             IF( SEPTMP( I ).NE.DUM( 1 ) ) VMAX = ONE / EPS
  130    CONTINUE

         // Compute eigenvector condition numbers only and compare

         CALL SCOPY( N, DUM, 0, STMP, 1 )
         CALL SCOPY( N, DUM, 0, SEPTMP, 1 )
         CALL STRSNA( 'Veccond', 'All', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         DO 140 I = 1, N
            IF( STMP( I ).NE.DUM( 1 ) ) VMAX = ONE / EPS             IF( SEPTMP( I ).NE.SEP( I ) ) VMAX = ONE / EPS
  140    CONTINUE

         // Compute all condition numbers using SELECT and compare

         DO 150 I = 1, N
            SELECT( I ) = .TRUE.
  150    CONTINUE
         CALL SCOPY( N, DUM, 0, STMP, 1 )
         CALL SCOPY( N, DUM, 0, SEPTMP, 1 )
         CALL STRSNA( 'Bothcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         DO 160 I = 1, N
            IF( SEPTMP( I ).NE.SEP( I ) ) VMAX = ONE / EPS             IF( STMP( I ).NE.S( I ) ) VMAX = ONE / EPS
  160    CONTINUE

         // Compute eigenvalue condition numbers using SELECT and compare

         CALL SCOPY( N, DUM, 0, STMP, 1 )
         CALL SCOPY( N, DUM, 0, SEPTMP, 1 )
         CALL STRSNA( 'Eigcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         DO 170 I = 1, N
            IF( STMP( I ).NE.S( I ) ) VMAX = ONE / EPS             IF( SEPTMP( I ).NE.DUM( 1 ) ) VMAX = ONE / EPS
  170    CONTINUE

         // Compute eigenvector condition numbers using SELECT and compare

         CALL SCOPY( N, DUM, 0, STMP, 1 )
         CALL SCOPY( N, DUM, 0, SEPTMP, 1 )
         CALL STRSNA( 'Veccond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         DO 180 I = 1, N
            IF( STMP( I ).NE.DUM( 1 ) ) VMAX = ONE / EPS             IF( SEPTMP( I ).NE.SEP( I ) ) VMAX = ONE / EPS
  180    CONTINUE
         if ( VMAX.GT.RMAX( 1 ) ) {
            RMAX( 1 ) = VMAX
            IF( NINFO( 1 ).EQ.0 ) LMAX( 1 ) = KNT
         }

         // Select first real and first complex eigenvalue

         if ( WI( 1 ).EQ.ZERO ) {
            LCMP( 1 ) = 1
            IFND = 0
            DO 190 I = 2, N
               if ( IFND.EQ.1 .OR. WI( I ).EQ.ZERO ) {
                  SELECT( I ) = .FALSE.
               } else {
                  IFND = 1
                  LCMP( 2 ) = I
                  LCMP( 3 ) = I + 1
                  CALL SCOPY( N, RE( 1, I ), 1, RE( 1, 2 ), 1 )
                  CALL SCOPY( N, RE( 1, I+1 ), 1, RE( 1, 3 ), 1 )
                  CALL SCOPY( N, LE( 1, I ), 1, LE( 1, 2 ), 1 )
                  CALL SCOPY( N, LE( 1, I+1 ), 1, LE( 1, 3 ), 1 )
               }
  190       CONTINUE
            if ( IFND.EQ.0 ) {
               ICMP = 1
            } else {
               ICMP = 3
            }
         } else {
            LCMP( 1 ) = 1
            LCMP( 2 ) = 2
            IFND = 0
            DO 200 I = 3, N
               if ( IFND.EQ.1 .OR. WI( I ).NE.ZERO ) {
                  SELECT( I ) = .FALSE.
               } else {
                  LCMP( 3 ) = I
                  IFND = 1
                  CALL SCOPY( N, RE( 1, I ), 1, RE( 1, 3 ), 1 )
                  CALL SCOPY( N, LE( 1, I ), 1, LE( 1, 3 ), 1 )
               }
  200       CONTINUE
            if ( IFND.EQ.0 ) {
               ICMP = 2
            } else {
               ICMP = 3
            }
         }

         // Compute all selected condition numbers

         CALL SCOPY( ICMP, DUM, 0, STMP, 1 )
         CALL SCOPY( ICMP, DUM, 0, SEPTMP, 1 )
         CALL STRSNA( 'Bothcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         DO 210 I = 1, ICMP
            J = LCMP( I )
            IF( SEPTMP( I ).NE.SEP( J ) ) VMAX = ONE / EPS             IF( STMP( I ).NE.S( J ) ) VMAX = ONE / EPS
  210    CONTINUE

         // Compute selected eigenvalue condition numbers

         CALL SCOPY( ICMP, DUM, 0, STMP, 1 )
         CALL SCOPY( ICMP, DUM, 0, SEPTMP, 1 )
         CALL STRSNA( 'Eigcond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         DO 220 I = 1, ICMP
            J = LCMP( I )
            IF( STMP( I ).NE.S( J ) ) VMAX = ONE / EPS             IF( SEPTMP( I ).NE.DUM( 1 ) ) VMAX = ONE / EPS
  220    CONTINUE

         // Compute selected eigenvector condition numbers

         CALL SCOPY( ICMP, DUM, 0, STMP, 1 )
         CALL SCOPY( ICMP, DUM, 0, SEPTMP, 1 )
         CALL STRSNA( 'Veccond', 'Some', SELECT, N, T, LDT, LE, LDT, RE, LDT, STMP, SEPTMP, N, M, WORK, N, IWORK, INFO )
         if ( INFO.NE.0 ) {
            LMAX( 3 ) = KNT
            NINFO( 3 ) = NINFO( 3 ) + 1
            GO TO 240
         }
         DO 230 I = 1, ICMP
            J = LCMP( I )
            IF( STMP( I ).NE.DUM( 1 ) ) VMAX = ONE / EPS             IF( SEPTMP( I ).NE.SEP( J ) ) VMAX = ONE / EPS
  230    CONTINUE
         if ( VMAX.GT.RMAX( 1 ) ) {
            RMAX( 1 ) = VMAX
            IF( NINFO( 1 ).EQ.0 ) LMAX( 1 ) = KNT
         }
  240 CONTINUE
      GO TO 10

      // End of SGET37

      }
