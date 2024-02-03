      SUBROUTINE SLATMR( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RSIGN, GRADE, DL, MODEL, CONDL, DR, MODER, CONDR, PIVTNG, IPIVOT, KL, KU, SPARSE, ANORM, PACK, A, LDA, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIST, GRADE, PACK, PIVTNG, RSIGN, SYM;
      int                INFO, KL, KU, LDA, M, MODE, MODEL, MODER, N;
      REAL               ANORM, COND, CONDL, CONDR, DMAX, SPARSE
      // ..
      // .. Array Arguments ..
      int                IPIVOT( * ), ISEED( 4 ), IWORK( * );
      REAL               A( LDA, * ), D( * ), DL( * ), DR( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      REAL               ONE
      const              ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               BADPVT, DZERO, FULBND;
      int                I, IDIST, IGRADE, IISUB, IPACK, IPVTNG, IRSIGN, ISUB, ISYM, J, JJSUB, JSUB, K, KLL, KUU, MNMIN, MNSUB, MXSUB, NPVTS;
      REAL               ALPHA, ONORM, TEMP
      // ..
      // .. Local Arrays ..
      REAL               TEMPA( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLANGB, SLANGE, SLANSB, SLANSP, SLANSY, SLATM2, SLATM3       EXTERNAL           LSAME, SLANGB, SLANGE, SLANSB, SLANSP, SLANSY, SLATM2, SLATM3
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLATM1, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, MOD
      // ..
      // .. Executable Statements ..

      // 1)      Decode and Test the input parameters.
              // Initialize flags & seed.

      INFO = 0

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      // Decode DIST

      if ( LSAME( DIST, 'U' ) ) {
         IDIST = 1
      } else if ( LSAME( DIST, 'S' ) ) {
         IDIST = 2
      } else if ( LSAME( DIST, 'N' ) ) {
         IDIST = 3
      } else {
         IDIST = -1
      }

      // Decode SYM

      if ( LSAME( SYM, 'S' ) ) {
         ISYM = 0
      } else if ( LSAME( SYM, 'N' ) ) {
         ISYM = 1
      } else if ( LSAME( SYM, 'H' ) ) {
         ISYM = 0
      } else {
         ISYM = -1
      }

      // Decode RSIGN

      if ( LSAME( RSIGN, 'F' ) ) {
         IRSIGN = 0
      } else if ( LSAME( RSIGN, 'T' ) ) {
         IRSIGN = 1
      } else {
         IRSIGN = -1
      }

      // Decode PIVTNG

      if ( LSAME( PIVTNG, 'N' ) ) {
         IPVTNG = 0
      } else if ( LSAME( PIVTNG, ' ' ) ) {
         IPVTNG = 0
      } else if ( LSAME( PIVTNG, 'L' ) ) {
         IPVTNG = 1
         NPVTS = M
      } else if ( LSAME( PIVTNG, 'R' ) ) {
         IPVTNG = 2
         NPVTS = N
      } else if ( LSAME( PIVTNG, 'B' ) ) {
         IPVTNG = 3
         NPVTS = MIN( N, M )
      } else if ( LSAME( PIVTNG, 'F' ) ) {
         IPVTNG = 3
         NPVTS = MIN( N, M )
      } else {
         IPVTNG = -1
      }

      // Decode GRADE

      if ( LSAME( GRADE, 'N' ) ) {
         IGRADE = 0
      } else if ( LSAME( GRADE, 'L' ) ) {
         IGRADE = 1
      } else if ( LSAME( GRADE, 'R' ) ) {
         IGRADE = 2
      } else if ( LSAME( GRADE, 'B' ) ) {
         IGRADE = 3
      } else if ( LSAME( GRADE, 'E' ) ) {
         IGRADE = 4
      } else if ( LSAME( GRADE, 'H' ) .OR. LSAME( GRADE, 'S' ) ) {
         IGRADE = 5
      } else {
         IGRADE = -1
      }

      // Decode PACK

      if ( LSAME( PACK, 'N' ) ) {
         IPACK = 0
      } else if ( LSAME( PACK, 'U' ) ) {
         IPACK = 1
      } else if ( LSAME( PACK, 'L' ) ) {
         IPACK = 2
      } else if ( LSAME( PACK, 'C' ) ) {
         IPACK = 3
      } else if ( LSAME( PACK, 'R' ) ) {
         IPACK = 4
      } else if ( LSAME( PACK, 'B' ) ) {
         IPACK = 5
      } else if ( LSAME( PACK, 'Q' ) ) {
         IPACK = 6
      } else if ( LSAME( PACK, 'Z' ) ) {
         IPACK = 7
      } else {
         IPACK = -1
      }

      // Set certain internal parameters

      MNMIN = MIN( M, N )
      KLL = MIN( KL, M-1 )
      KUU = MIN( KU, N-1 )

      // If inv(DL) is used, check to see if DL has a zero entry.

      DZERO = .FALSE.
      if ( IGRADE.EQ.4 .AND. MODEL.EQ.0 ) {
         DO 10 I = 1, M
            IF( DL( I ).EQ.ZERO ) DZERO = .TRUE.
   10    CONTINUE
      }

      // Check values in IPIVOT

      BADPVT = .FALSE.
      if ( IPVTNG.GT.0 ) {
         DO 20 J = 1, NPVTS
            IF( IPIVOT( J ).LE.0 .OR. IPIVOT( J ).GT.NPVTS ) BADPVT = .TRUE.
   20    CONTINUE
      }

      // Set INFO if an error

      if ( M.LT.0 ) {
         INFO = -1
      } else if ( M.NE.N .AND. ISYM.EQ.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( IDIST.EQ.-1 ) {
         INFO = -3
      } else if ( ISYM.EQ.-1 ) {
         INFO = -5
      } else if ( MODE.LT.-6 .OR. MODE.GT.6 ) {
         INFO = -7
      } else if ( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. COND.LT.ONE ) {
         INFO = -8
      } else if ( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. IRSIGN.EQ.-1 ) {
         INFO = -10
      } else if ( IGRADE.EQ.-1 .OR. ( IGRADE.EQ.4 .AND. M.NE.N ) .OR. ( ( IGRADE.GE.1 .AND. IGRADE.LE.4 ) .AND. ISYM.EQ.0 ) ) {
         INFO = -11
      } else if ( IGRADE.EQ.4 .AND. DZERO ) {
         INFO = -12
      } else if ( ( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ.5 ) .AND. ( MODEL.LT.-6 .OR. MODEL.GT.6 ) ) {
         INFO = -13
      } else if ( ( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ.5 ) .AND. ( MODEL.NE.-6 .AND. MODEL.NE.0 .AND. MODEL.NE.6 ) .AND. CONDL.LT.ONE ) {
         INFO = -14
      } else if ( ( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) .AND. ( MODER.LT.-6 .OR. MODER.GT.6 ) ) {
         INFO = -16
      } else if ( ( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) .AND. ( MODER.NE.-6 .AND. MODER.NE.0 .AND. MODER.NE.6 ) .AND. CONDR.LT.ONE ) {
         INFO = -17
      } else if ( IPVTNG.EQ.-1 .OR. ( IPVTNG.EQ.3 .AND. M.NE.N ) .OR. ( ( IPVTNG.EQ.1 .OR. IPVTNG.EQ.2 ) .AND. ISYM.EQ.0 ) ) {
         INFO = -18
      } else if ( IPVTNG.NE.0 .AND. BADPVT ) {
         INFO = -19
      } else if ( KL.LT.0 ) {
         INFO = -20
      } else if ( KU.LT.0 .OR. ( ISYM.EQ.0 .AND. KL.NE.KU ) ) {
         INFO = -21
      } else if ( SPARSE.LT.ZERO .OR. SPARSE.GT.ONE ) {
         INFO = -22
      } else if ( IPACK.EQ.-1 .OR. ( ( IPACK.EQ.1 .OR. IPACK.EQ.2 .OR. IPACK.EQ.5 .OR. IPACK.EQ.6 ) .AND. ISYM.EQ.1 ) .OR. ( IPACK.EQ.3 .AND. ISYM.EQ.1 .AND. ( KL.NE.0 .OR. M.NE. N ) ) .OR. ( IPACK.EQ.4 .AND. ISYM.EQ.1 .AND. ( KU.NE. 0 .OR. M.NE.N ) ) ) {
         INFO = -24
      } else if ( ( ( IPACK.EQ.0 .OR. IPACK.EQ.1 .OR. IPACK.EQ.2 ) .AND. LDA.LT.MAX( 1, M ) ) .OR. ( ( IPACK.EQ.3 .OR. IPACK.EQ. 4 ) .AND. LDA.LT.1 ) .OR. ( ( IPACK.EQ.5 .OR. IPACK.EQ. 6 ) .AND. LDA.LT.KUU+1 ) .OR. ( IPACK.EQ.7 .AND. LDA.LT.KLL+KUU+1 ) ) {
         INFO = -26
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SLATMR', -INFO )
         RETURN
      }

      // Decide if we can pivot consistently

      FULBND = .FALSE.
      IF( KUU.EQ.N-1 .AND. KLL.EQ.M-1 ) FULBND = .TRUE.

      // Initialize random number generator

      DO 30 I = 1, 4
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
   30 CONTINUE

      ISEED( 4 ) = 2*( ISEED( 4 ) / 2 ) + 1

      // 2)      Set up D, DL, and DR, if indicated.

              // Compute D according to COND and MODE

      CALL SLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, INFO )
      if ( INFO.NE.0 ) {
         INFO = 1
         RETURN
      }
      if ( MODE.NE.0 .AND. MODE.NE.-6 .AND. MODE.NE.6 ) {

         // Scale by DMAX

         TEMP = ABS( D( 1 ) )
         DO 40 I = 2, MNMIN
            TEMP = MAX( TEMP, ABS( D( I ) ) )
   40    CONTINUE
         if ( TEMP.EQ.ZERO .AND. DMAX.NE.ZERO ) {
            INFO = 2
            RETURN
         }
         if ( TEMP.NE.ZERO ) {
            ALPHA = DMAX / TEMP
         } else {
            ALPHA = ONE
         }
         DO 50 I = 1, MNMIN
            D( I ) = ALPHA*D( I )
   50    CONTINUE

      }

      // Compute DL if grading set

      if ( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ. 5 ) {
         CALL SLATM1( MODEL, CONDL, 0, IDIST, ISEED, DL, M, INFO )
         if ( INFO.NE.0 ) {
            INFO = 3
            RETURN
         }
      }

      // Compute DR if grading set

      if ( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) {
         CALL SLATM1( MODER, CONDR, 0, IDIST, ISEED, DR, N, INFO )
         if ( INFO.NE.0 ) {
            INFO = 4
            RETURN
         }
      }

      // 3)     Generate IWORK if pivoting

      if ( IPVTNG.GT.0 ) {
         DO 60 I = 1, NPVTS
            IWORK( I ) = I
   60    CONTINUE
         if ( FULBND ) {
            DO 70 I = 1, NPVTS
               K = IPIVOT( I )
               J = IWORK( I )
               IWORK( I ) = IWORK( K )
               IWORK( K ) = J
   70       CONTINUE
         } else {
            DO 80 I = NPVTS, 1, -1
               K = IPIVOT( I )
               J = IWORK( I )
               IWORK( I ) = IWORK( K )
               IWORK( K ) = J
   80       CONTINUE
         }
      }

      // 4)      Generate matrices for each kind of PACKing
              // Always sweep matrix columnwise (if symmetric, upper
              // half only) so that matrix generated does not depend
              // on PACK

      if ( FULBND ) {

         // Use SLATM3 so matrices generated with differing PIVOTing only
         // differ only in the order of their rows and/or columns.

         if ( IPACK.EQ.0 ) {
            if ( ISYM.EQ.0 ) {
               DO 100 J = 1, N
                  DO 90 I = 1, J
                     TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = TEMP
                     A( JSUB, ISUB ) = TEMP
   90             CONTINUE
  100          CONTINUE
            } else if ( ISYM.EQ.1 ) {
               DO 120 J = 1, N
                  DO 110 I = 1, M
                     TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = TEMP
  110             CONTINUE
  120          CONTINUE
            }

         } else if ( IPACK.EQ.1 ) {

            DO 140 J = 1, N
               DO 130 I = 1, J
                  TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  A( MNSUB, MXSUB ) = TEMP
                  IF( MNSUB.NE.MXSUB ) A( MXSUB, MNSUB ) = ZERO
  130          CONTINUE
  140       CONTINUE

         } else if ( IPACK.EQ.2 ) {

            DO 160 J = 1, N
               DO 150 I = 1, J
                  TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  A( MXSUB, MNSUB ) = TEMP
                  IF( MNSUB.NE.MXSUB ) A( MNSUB, MXSUB ) = ZERO
  150          CONTINUE
  160       CONTINUE

         } else if ( IPACK.EQ.3 ) {

            DO 180 J = 1, N
               DO 170 I = 1, J
                  TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )

                  // Compute K = location of (ISUB,JSUB) entry in packed
                  // array

                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  K = MXSUB*( MXSUB-1 ) / 2 + MNSUB

                  // Convert K to (IISUB,JJSUB) location

                  JJSUB = ( K-1 ) / LDA + 1
                  IISUB = K - LDA*( JJSUB-1 )

                  A( IISUB, JJSUB ) = TEMP
  170          CONTINUE
  180       CONTINUE

         } else if ( IPACK.EQ.4 ) {

            DO 200 J = 1, N
               DO 190 I = 1, J
                  TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )

                  // Compute K = location of (I,J) entry in packed array

                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  if ( MNSUB.EQ.1 ) {
                     K = MXSUB
                  } else {
                     K = N*( N+1 ) / 2 - ( N-MNSUB+1 )*( N-MNSUB+2 ) / 2 + MXSUB - MNSUB + 1
                  }

                  // Convert K to (IISUB,JJSUB) location

                  JJSUB = ( K-1 ) / LDA + 1
                  IISUB = K - LDA*( JJSUB-1 )

                  A( IISUB, JJSUB ) = TEMP
  190          CONTINUE
  200       CONTINUE

         } else if ( IPACK.EQ.5 ) {

            DO 220 J = 1, N
               DO 210 I = J - KUU, J
                  if ( I.LT.1 ) {
                     A( J-I+1, I+N ) = ZERO
                  } else {
                     TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     MNSUB = MIN( ISUB, JSUB )
                     MXSUB = MAX( ISUB, JSUB )
                     A( MXSUB-MNSUB+1, MNSUB ) = TEMP
                  }
  210          CONTINUE
  220       CONTINUE

         } else if ( IPACK.EQ.6 ) {

            DO 240 J = 1, N
               DO 230 I = J - KUU, J
                  TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  A( MNSUB-MXSUB+KUU+1, MXSUB ) = TEMP
  230          CONTINUE
  240       CONTINUE

         } else if ( IPACK.EQ.7 ) {

            if ( ISYM.EQ.0 ) {
               DO 260 J = 1, N
                  DO 250 I = J - KUU, J
                     TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     MNSUB = MIN( ISUB, JSUB )
                     MXSUB = MAX( ISUB, JSUB )
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = TEMP
                     IF( I.LT.1 ) A( J-I+1+KUU, I+N ) = ZERO                      IF( I.GE.1 .AND. MNSUB.NE.MXSUB ) A( MXSUB-MNSUB+1+KUU, MNSUB ) = TEMP
  250             CONTINUE
  260          CONTINUE
            } else if ( ISYM.EQ.1 ) {
               DO 280 J = 1, N
                  DO 270 I = J - KUU, J + KLL
                     TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB-JSUB+KUU+1, JSUB ) = TEMP
  270             CONTINUE
  280          CONTINUE
            }

         }

      } else {

         // Use SLATM2

         if ( IPACK.EQ.0 ) {
            if ( ISYM.EQ.0 ) {
               DO 300 J = 1, N
                  DO 290 I = 1, J
                     A( I, J ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( J, I ) = A( I, J )
  290             CONTINUE
  300          CONTINUE
            } else if ( ISYM.EQ.1 ) {
               DO 320 J = 1, N
                  DO 310 I = 1, M
                     A( I, J ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  310             CONTINUE
  320          CONTINUE
            }

         } else if ( IPACK.EQ.1 ) {

            DO 340 J = 1, N
               DO 330 I = 1, J
                  A( I, J ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )                   IF( I.NE.J ) A( J, I ) = ZERO
  330          CONTINUE
  340       CONTINUE

         } else if ( IPACK.EQ.2 ) {

            DO 360 J = 1, N
               DO 350 I = 1, J
                  A( J, I ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )                   IF( I.NE.J ) A( I, J ) = ZERO
  350          CONTINUE
  360       CONTINUE

         } else if ( IPACK.EQ.3 ) {

            ISUB = 0
            JSUB = 1
            DO 380 J = 1, N
               DO 370 I = 1, J
                  ISUB = ISUB + 1
                  if ( ISUB.GT.LDA ) {
                     ISUB = 1
                     JSUB = JSUB + 1
                  }
                  A( ISUB, JSUB ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  370          CONTINUE
  380       CONTINUE

         } else if ( IPACK.EQ.4 ) {

            if ( ISYM.EQ.0 ) {
               DO 400 J = 1, N
                  DO 390 I = 1, J

                     // Compute K = location of (I,J) entry in packed array

                     if ( I.EQ.1 ) {
                        K = J
                     } else {
                        K = N*( N+1 ) / 2 - ( N-I+1 )*( N-I+2 ) / 2 + J - I + 1
                     }

                     // Convert K to (ISUB,JSUB) location

                     JSUB = ( K-1 ) / LDA + 1
                     ISUB = K - LDA*( JSUB-1 )

                     A( ISUB, JSUB ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  390             CONTINUE
  400          CONTINUE
            } else {
               ISUB = 0
               JSUB = 1
               DO 420 J = 1, N
                  DO 410 I = J, M
                     ISUB = ISUB + 1
                     if ( ISUB.GT.LDA ) {
                        ISUB = 1
                        JSUB = JSUB + 1
                     }
                     A( ISUB, JSUB ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  410             CONTINUE
  420          CONTINUE
            }

         } else if ( IPACK.EQ.5 ) {

            DO 440 J = 1, N
               DO 430 I = J - KUU, J
                  if ( I.LT.1 ) {
                     A( J-I+1, I+N ) = ZERO
                  } else {
                     A( J-I+1, I ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  }
  430          CONTINUE
  440       CONTINUE

         } else if ( IPACK.EQ.6 ) {

            DO 460 J = 1, N
               DO 450 I = J - KUU, J
                  A( I-J+KUU+1, J ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  450          CONTINUE
  460       CONTINUE

         } else if ( IPACK.EQ.7 ) {

            if ( ISYM.EQ.0 ) {
               DO 480 J = 1, N
                  DO 470 I = J - KUU, J
                     A( I-J+KUU+1, J ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     IF( I.LT.1 ) A( J-I+1+KUU, I+N ) = ZERO                      IF( I.GE.1 .AND. I.NE.J ) A( J-I+1+KUU, I ) = A( I-J+KUU+1, J )
  470             CONTINUE
  480          CONTINUE
            } else if ( ISYM.EQ.1 ) {
               DO 500 J = 1, N
                  DO 490 I = J - KUU, J + KLL
                     A( I-J+KUU+1, J ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  490             CONTINUE
  500          CONTINUE
            }

         }

      }

      // 5)      Scaling the norm

      if ( IPACK.EQ.0 ) {
         ONORM = SLANGE( 'M', M, N, A, LDA, TEMPA )
      } else if ( IPACK.EQ.1 ) {
         ONORM = SLANSY( 'M', 'U', N, A, LDA, TEMPA )
      } else if ( IPACK.EQ.2 ) {
         ONORM = SLANSY( 'M', 'L', N, A, LDA, TEMPA )
      } else if ( IPACK.EQ.3 ) {
         ONORM = SLANSP( 'M', 'U', N, A, TEMPA )
      } else if ( IPACK.EQ.4 ) {
         ONORM = SLANSP( 'M', 'L', N, A, TEMPA )
      } else if ( IPACK.EQ.5 ) {
         ONORM = SLANSB( 'M', 'L', N, KLL, A, LDA, TEMPA )
      } else if ( IPACK.EQ.6 ) {
         ONORM = SLANSB( 'M', 'U', N, KUU, A, LDA, TEMPA )
      } else if ( IPACK.EQ.7 ) {
         ONORM = SLANGB( 'M', N, KLL, KUU, A, LDA, TEMPA )
      }

      if ( ANORM.GE.ZERO ) {

         if ( ANORM.GT.ZERO .AND. ONORM.EQ.ZERO ) {

            // Desired scaling impossible

            INFO = 5
            RETURN

         } else if ( ( ANORM.GT.ONE .AND. ONORM.LT.ONE ) .OR. ( ANORM.LT.ONE .AND. ONORM.GT.ONE ) ) {

            // Scale carefully to avoid over / underflow

            if ( IPACK.LE.2 ) {
               DO 510 J = 1, N
                  CALL SSCAL( M, ONE / ONORM, A( 1, J ), 1 )
                  CALL SSCAL( M, ANORM, A( 1, J ), 1 )
  510          CONTINUE

            } else if ( IPACK.EQ.3 .OR. IPACK.EQ.4 ) {

               CALL SSCAL( N*( N+1 ) / 2, ONE / ONORM, A, 1 )
               CALL SSCAL( N*( N+1 ) / 2, ANORM, A, 1 )

            } else if ( IPACK.GE.5 ) {

               DO 520 J = 1, N
                  CALL SSCAL( KLL+KUU+1, ONE / ONORM, A( 1, J ), 1 )
                  CALL SSCAL( KLL+KUU+1, ANORM, A( 1, J ), 1 )
  520          CONTINUE

            }

         } else {

            // Scale straightforwardly

            if ( IPACK.LE.2 ) {
               DO 530 J = 1, N
                  CALL SSCAL( M, ANORM / ONORM, A( 1, J ), 1 )
  530          CONTINUE

            } else if ( IPACK.EQ.3 .OR. IPACK.EQ.4 ) {

               CALL SSCAL( N*( N+1 ) / 2, ANORM / ONORM, A, 1 )

            } else if ( IPACK.GE.5 ) {

               DO 540 J = 1, N
                  CALL SSCAL( KLL+KUU+1, ANORM / ONORM, A( 1, J ), 1 )
  540          CONTINUE
            }

         }

      }

      // End of SLATMR

      }
