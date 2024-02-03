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

      if (M.EQ.0 .OR. N.EQ.0) RETURN;

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

      DZERO = false;
      if ( IGRADE.EQ.4 .AND. MODEL.EQ.0 ) {
         for (I = 1; I <= M; I++) { // 10
            IF( DL( I ).EQ.ZERO ) DZERO = true;
         } // 10
      }

      // Check values in IPIVOT

      BADPVT = false;
      if ( IPVTNG.GT.0 ) {
         for (J = 1; J <= NPVTS; J++) { // 20
            IF( IPIVOT( J ).LE.0 .OR. IPIVOT( J ).GT.NPVTS ) BADPVT = true;
         } // 20
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
         xerbla('SLATMR', -INFO );
         RETURN
      }

      // Decide if we can pivot consistently

      FULBND = false;
      if (KUU.EQ.N-1 .AND. KLL.EQ.M-1) FULBND = true ;

      // Initialize random number generator

      for (I = 1; I <= 4; I++) { // 30
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
      } // 30

      ISEED( 4 ) = 2*( ISEED( 4 ) / 2 ) + 1

      // 2)      Set up D, DL, and DR, if indicated.

              // Compute D according to COND and MODE

      slatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, INFO );
      if ( INFO.NE.0 ) {
         INFO = 1
         RETURN
      }
      if ( MODE.NE.0 .AND. MODE.NE.-6 .AND. MODE.NE.6 ) {

         // Scale by DMAX

         TEMP = ABS( D( 1 ) )
         for (I = 2; I <= MNMIN; I++) { // 40
            TEMP = MAX( TEMP, ABS( D( I ) ) )
         } // 40
         if ( TEMP.EQ.ZERO .AND. DMAX.NE.ZERO ) {
            INFO = 2
            RETURN
         }
         if ( TEMP.NE.ZERO ) {
            ALPHA = DMAX / TEMP
         } else {
            ALPHA = ONE
         }
         for (I = 1; I <= MNMIN; I++) { // 50
            D( I ) = ALPHA*D( I )
         } // 50

      }

      // Compute DL if grading set

      if ( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ. 5 ) {
         slatm1(MODEL, CONDL, 0, IDIST, ISEED, DL, M, INFO );
         if ( INFO.NE.0 ) {
            INFO = 3
            RETURN
         }
      }

      // Compute DR if grading set

      if ( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) {
         slatm1(MODER, CONDR, 0, IDIST, ISEED, DR, N, INFO );
         if ( INFO.NE.0 ) {
            INFO = 4
            RETURN
         }
      }

      // 3)     Generate IWORK if pivoting

      if ( IPVTNG.GT.0 ) {
         for (I = 1; I <= NPVTS; I++) { // 60
            IWORK( I ) = I
         } // 60
         if ( FULBND ) {
            for (I = 1; I <= NPVTS; I++) { // 70
               K = IPIVOT( I )
               J = IWORK( I )
               IWORK( I ) = IWORK( K )
               IWORK( K ) = J
            } // 70
         } else {
            DO 80 I = NPVTS, 1, -1
               K = IPIVOT( I )
               J = IWORK( I )
               IWORK( I ) = IWORK( K )
               IWORK( K ) = J
            } // 80
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
               for (J = 1; J <= N; J++) { // 100
                  for (I = 1; I <= J; I++) { // 90
                     TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = TEMP
                     A( JSUB, ISUB ) = TEMP
                  } // 90
               } // 100
            } else if ( ISYM.EQ.1 ) {
               for (J = 1; J <= N; J++) { // 120
                  for (I = 1; I <= M; I++) { // 110
                     TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = TEMP
                  } // 110
               } // 120
            }

         } else if ( IPACK.EQ.1 ) {

            for (J = 1; J <= N; J++) { // 140
               for (I = 1; I <= J; I++) { // 130
                  TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  A( MNSUB, MXSUB ) = TEMP
                  if (MNSUB.NE.MXSUB) A( MXSUB, MNSUB ) = ZERO;
               } // 130
            } // 140

         } else if ( IPACK.EQ.2 ) {

            for (J = 1; J <= N; J++) { // 160
               for (I = 1; I <= J; I++) { // 150
                  TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  A( MXSUB, MNSUB ) = TEMP
                  if (MNSUB.NE.MXSUB) A( MNSUB, MXSUB ) = ZERO;
               } // 150
            } // 160

         } else if ( IPACK.EQ.3 ) {

            for (J = 1; J <= N; J++) { // 180
               for (I = 1; I <= J; I++) { // 170
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
               } // 170
            } // 180

         } else if ( IPACK.EQ.4 ) {

            for (J = 1; J <= N; J++) { // 200
               for (I = 1; I <= J; I++) { // 190
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
               } // 190
            } // 200

         } else if ( IPACK.EQ.5 ) {

            for (J = 1; J <= N; J++) { // 220
               for (I = J - KUU; I <= J; I++) { // 210
                  if ( I.LT.1 ) {
                     A( J-I+1, I+N ) = ZERO
                  } else {
                     TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     MNSUB = MIN( ISUB, JSUB )
                     MXSUB = MAX( ISUB, JSUB )
                     A( MXSUB-MNSUB+1, MNSUB ) = TEMP
                  }
               } // 210
            } // 220

         } else if ( IPACK.EQ.6 ) {

            for (J = 1; J <= N; J++) { // 240
               for (I = J - KUU; I <= J; I++) { // 230
                  TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  A( MNSUB-MXSUB+KUU+1, MXSUB ) = TEMP
               } // 230
            } // 240

         } else if ( IPACK.EQ.7 ) {

            if ( ISYM.EQ.0 ) {
               for (J = 1; J <= N; J++) { // 260
                  for (I = J - KUU; I <= J; I++) { // 250
                     TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     MNSUB = MIN( ISUB, JSUB )
                     MXSUB = MAX( ISUB, JSUB )
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = TEMP
                     if (I.LT.1) A( J-I+1+KUU, I+N ) = ZERO                      IF( I.GE.1 .AND. MNSUB.NE.MXSUB ) A( MXSUB-MNSUB+1+KUU, MNSUB ) = TEMP;
                  } // 250
               } // 260
            } else if ( ISYM.EQ.1 ) {
               for (J = 1; J <= N; J++) { // 280
                  for (I = J - KUU; I <= J + KLL; I++) { // 270
                     TEMP = SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB-JSUB+KUU+1, JSUB ) = TEMP
                  } // 270
               } // 280
            }

         }

      } else {

         // Use SLATM2

         if ( IPACK.EQ.0 ) {
            if ( ISYM.EQ.0 ) {
               for (J = 1; J <= N; J++) { // 300
                  for (I = 1; I <= J; I++) { // 290
                     A( I, J ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( J, I ) = A( I, J )
                  } // 290
               } // 300
            } else if ( ISYM.EQ.1 ) {
               for (J = 1; J <= N; J++) { // 320
                  for (I = 1; I <= M; I++) { // 310
                     A( I, J ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  } // 310
               } // 320
            }

         } else if ( IPACK.EQ.1 ) {

            for (J = 1; J <= N; J++) { // 340
               for (I = 1; I <= J; I++) { // 330
                  A( I, J ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )                   IF( I.NE.J ) A( J, I ) = ZERO
               } // 330
            } // 340

         } else if ( IPACK.EQ.2 ) {

            for (J = 1; J <= N; J++) { // 360
               for (I = 1; I <= J; I++) { // 350
                  A( J, I ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )                   IF( I.NE.J ) A( I, J ) = ZERO
               } // 350
            } // 360

         } else if ( IPACK.EQ.3 ) {

            ISUB = 0
            JSUB = 1
            for (J = 1; J <= N; J++) { // 380
               for (I = 1; I <= J; I++) { // 370
                  ISUB = ISUB + 1
                  if ( ISUB.GT.LDA ) {
                     ISUB = 1
                     JSUB = JSUB + 1
                  }
                  A( ISUB, JSUB ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
               } // 370
            } // 380

         } else if ( IPACK.EQ.4 ) {

            if ( ISYM.EQ.0 ) {
               for (J = 1; J <= N; J++) { // 400
                  for (I = 1; I <= J; I++) { // 390

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
                  } // 390
               } // 400
            } else {
               ISUB = 0
               JSUB = 1
               for (J = 1; J <= N; J++) { // 420
                  for (I = J; I <= M; I++) { // 410
                     ISUB = ISUB + 1
                     if ( ISUB.GT.LDA ) {
                        ISUB = 1
                        JSUB = JSUB + 1
                     }
                     A( ISUB, JSUB ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  } // 410
               } // 420
            }

         } else if ( IPACK.EQ.5 ) {

            for (J = 1; J <= N; J++) { // 440
               for (I = J - KUU; I <= J; I++) { // 430
                  if ( I.LT.1 ) {
                     A( J-I+1, I+N ) = ZERO
                  } else {
                     A( J-I+1, I ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  }
               } // 430
            } // 440

         } else if ( IPACK.EQ.6 ) {

            for (J = 1; J <= N; J++) { // 460
               for (I = J - KUU; I <= J; I++) { // 450
                  A( I-J+KUU+1, J ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
               } // 450
            } // 460

         } else if ( IPACK.EQ.7 ) {

            if ( ISYM.EQ.0 ) {
               for (J = 1; J <= N; J++) { // 480
                  for (I = J - KUU; I <= J; I++) { // 470
                     A( I-J+KUU+1, J ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     if (I.LT.1) A( J-I+1+KUU, I+N ) = ZERO                      IF( I.GE.1 .AND. I.NE.J ) A( J-I+1+KUU, I ) = A( I-J+KUU+1, J );
                  } // 470
               } // 480
            } else if ( ISYM.EQ.1 ) {
               for (J = 1; J <= N; J++) { // 500
                  for (I = J - KUU; I <= J + KLL; I++) { // 490
                     A( I-J+KUU+1, J ) = SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  } // 490
               } // 500
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
               for (J = 1; J <= N; J++) { // 510
                  sscal(M, ONE / ONORM, A( 1, J ), 1 );
                  sscal(M, ANORM, A( 1, J ), 1 );
               } // 510

            } else if ( IPACK.EQ.3 .OR. IPACK.EQ.4 ) {

               sscal(N*( N+1 ) / 2, ONE / ONORM, A, 1 );
               sscal(N*( N+1 ) / 2, ANORM, A, 1 );

            } else if ( IPACK.GE.5 ) {

               for (J = 1; J <= N; J++) { // 520
                  sscal(KLL+KUU+1, ONE / ONORM, A( 1, J ), 1 );
                  sscal(KLL+KUU+1, ANORM, A( 1, J ), 1 );
               } // 520

            }

         } else {

            // Scale straightforwardly

            if ( IPACK.LE.2 ) {
               for (J = 1; J <= N; J++) { // 530
                  sscal(M, ANORM / ONORM, A( 1, J ), 1 );
               } // 530

            } else if ( IPACK.EQ.3 .OR. IPACK.EQ.4 ) {

               sscal(N*( N+1 ) / 2, ANORM / ONORM, A, 1 );

            } else if ( IPACK.GE.5 ) {

               for (J = 1; J <= N; J++) { // 540
                  sscal(KLL+KUU+1, ANORM / ONORM, A( 1, J ), 1 );
               } // 540
            }

         }

      }

      // End of SLATMR

      }
