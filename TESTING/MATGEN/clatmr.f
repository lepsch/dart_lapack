      SUBROUTINE CLATMR( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RSIGN, GRADE, DL, MODEL, CONDL, DR, MODER, CONDR, PIVTNG, IPIVOT, KL, KU, SPARSE, ANORM, PACK, A, LDA, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIST, GRADE, PACK, PIVTNG, RSIGN, SYM;
      int                INFO, KL, KU, LDA, M, MODE, MODEL, MODER, N;
      REAL               ANORM, COND, CONDL, CONDR, SPARSE
      COMPLEX            DMAX
      // ..
      // .. Array Arguments ..
      int                IPIVOT( * ), ISEED( 4 ), IWORK( * );
      COMPLEX            A( LDA, * ), D( * ), DL( * ), DR( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      REAL               ONE
      const              ONE = 1.0E0 ;
      COMPLEX            CONE
      const              CONE = ( 1.0E0, 0.0E0 ) ;
      COMPLEX            CZERO
      const              CZERO = ( 0.0E0, 0.0E0 ) ;
      // ..
      // .. Local Scalars ..
      bool               BADPVT, DZERO, FULBND;
      int                I, IDIST, IGRADE, IISUB, IPACK, IPVTNG, IRSIGN, ISUB, ISYM, J, JJSUB, JSUB, K, KLL, KUU, MNMIN, MNSUB, MXSUB, NPVTS;
      REAL               ONORM, TEMP
      COMPLEX            CALPHA, CTEMP
      // ..
      // .. Local Arrays ..
      REAL               TEMPA( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGB, CLANGE, CLANSB, CLANSP, CLANSY
      COMPLEX            CLATM2, CLATM3
      // EXTERNAL LSAME, CLANGB, CLANGE, CLANSB, CLANSP, CLANSY, CLATM2, CLATM3
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLATM1, CSSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, MIN, MOD, REAL
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
      } else if ( LSAME( DIST, 'D' ) ) {
         IDIST = 4
      } else {
         IDIST = -1
      }

      // Decode SYM

      if ( LSAME( SYM, 'H' ) ) {
         ISYM = 0
      } else if ( LSAME( SYM, 'N' ) ) {
         ISYM = 1
      } else if ( LSAME( SYM, 'S' ) ) {
         ISYM = 2
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
      } else if ( LSAME( GRADE, 'H' ) ) {
         IGRADE = 5
      } else if ( LSAME( GRADE, 'S' ) ) {
         IGRADE = 6
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
         for (I = 1; I <= M; I++) { // 10
            IF( DL( I ).EQ.CZERO ) DZERO = .TRUE.
   10    CONTINUE
      }

      // Check values in IPIVOT

      BADPVT = .FALSE.
      if ( IPVTNG.GT.0 ) {
         for (J = 1; J <= NPVTS; J++) { // 20
            IF( IPIVOT( J ).LE.0 .OR. IPIVOT( J ).GT.NPVTS ) BADPVT = .TRUE.
   20    CONTINUE
      }

      // Set INFO if an error

      if ( M.LT.0 ) {
         INFO = -1
      } else if ( M.NE.N .AND. ( ISYM.EQ.0 .OR. ISYM.EQ.2 ) ) {
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
      } else if ( IGRADE.EQ.-1 .OR. ( IGRADE.EQ.4 .AND. M.NE.N ) .OR. ( ( IGRADE.EQ.1 .OR. IGRADE.EQ.2 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ.6 ) .AND. ISYM.EQ.0 ) .OR. ( ( IGRADE.EQ.1 .OR. IGRADE.EQ.2 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ.5 ) .AND. ISYM.EQ.2 ) ) {
         INFO = -11
      } else if ( IGRADE.EQ.4 .AND. DZERO ) {
         INFO = -12
      } else if ( ( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ.5 .OR. IGRADE.EQ.6 ) .AND. ( MODEL.LT.-6 .OR. MODEL.GT.6 ) ) {
         INFO = -13
      } else if ( ( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ.5 .OR. IGRADE.EQ.6 ) .AND. ( MODEL.NE.-6 .AND. MODEL.NE.0 .AND. MODEL.NE.6 ) .AND. CONDL.LT.ONE ) {
         INFO = -14
      } else if ( ( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) .AND. ( MODER.LT.-6 .OR. MODER.GT.6 ) ) {
         INFO = -16
      } else if ( ( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) .AND. ( MODER.NE.-6 .AND. MODER.NE.0 .AND. MODER.NE.6 ) .AND. CONDR.LT.ONE ) {
         INFO = -17
      } else if ( IPVTNG.EQ.-1 .OR. ( IPVTNG.EQ.3 .AND. M.NE.N ) .OR. ( ( IPVTNG.EQ.1 .OR. IPVTNG.EQ.2 ) .AND. ( ISYM.EQ.0 .OR. ISYM.EQ.2 ) ) ) {
         INFO = -18
      } else if ( IPVTNG.NE.0 .AND. BADPVT ) {
         INFO = -19
      } else if ( KL.LT.0 ) {
         INFO = -20
      } else if ( KU.LT.0 .OR. ( ( ISYM.EQ.0 .OR. ISYM.EQ.2 ) .AND. KL.NE. KU ) ) {
         INFO = -21
      } else if ( SPARSE.LT.ZERO .OR. SPARSE.GT.ONE ) {
         INFO = -22
      } else if ( IPACK.EQ.-1 .OR. ( ( IPACK.EQ.1 .OR. IPACK.EQ.2 .OR. IPACK.EQ.5 .OR. IPACK.EQ.6 ) .AND. ISYM.EQ.1 ) .OR. ( IPACK.EQ.3 .AND. ISYM.EQ.1 .AND. ( KL.NE.0 .OR. M.NE. N ) ) .OR. ( IPACK.EQ.4 .AND. ISYM.EQ.1 .AND. ( KU.NE. 0 .OR. M.NE.N ) ) ) {
         INFO = -24
      } else if ( ( ( IPACK.EQ.0 .OR. IPACK.EQ.1 .OR. IPACK.EQ.2 ) .AND. LDA.LT.MAX( 1, M ) ) .OR. ( ( IPACK.EQ.3 .OR. IPACK.EQ. 4 ) .AND. LDA.LT.1 ) .OR. ( ( IPACK.EQ.5 .OR. IPACK.EQ. 6 ) .AND. LDA.LT.KUU+1 ) .OR. ( IPACK.EQ.7 .AND. LDA.LT.KLL+KUU+1 ) ) {
         INFO = -26
      }

      if ( INFO.NE.0 ) {
         xerbla('CLATMR', -INFO );
         RETURN
      }

      // Decide if we can pivot consistently

      FULBND = .FALSE.
      IF( KUU.EQ.N-1 .AND. KLL.EQ.M-1 ) FULBND = .TRUE.

      // Initialize random number generator

      for (I = 1; I <= 4; I++) { // 30
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
   30 CONTINUE

      ISEED( 4 ) = 2*( ISEED( 4 ) / 2 ) + 1

      // 2)      Set up D, DL, and DR, if indicated.

              // Compute D according to COND and MODE

      clatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, INFO );
      if ( INFO.NE.0 ) {
         INFO = 1
         RETURN
      }
      if ( MODE.NE.0 .AND. MODE.NE.-6 .AND. MODE.NE.6 ) {

         // Scale by DMAX

         TEMP = ABS( D( 1 ) )
         for (I = 2; I <= MNMIN; I++) { // 40
            TEMP = MAX( TEMP, ABS( D( I ) ) )
   40    CONTINUE
         if ( TEMP.EQ.ZERO .AND. DMAX.NE.CZERO ) {
            INFO = 2
            RETURN
         }
         if ( TEMP.NE.ZERO ) {
            CALPHA = DMAX / TEMP
         } else {
            CALPHA = CONE
         }
         for (I = 1; I <= MNMIN; I++) { // 50
            D( I ) = CALPHA*D( I )
   50    CONTINUE

      }

      // If matrix Hermitian, make D real

      if ( ISYM.EQ.0 ) {
         for (I = 1; I <= MNMIN; I++) { // 60
            D( I ) = REAL( D( I ) )
   60    CONTINUE
      }

      // Compute DL if grading set

      if ( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ. 5 .OR. IGRADE.EQ.6 ) {
         clatm1(MODEL, CONDL, 0, IDIST, ISEED, DL, M, INFO );
         if ( INFO.NE.0 ) {
            INFO = 3
            RETURN
         }
      }

      // Compute DR if grading set

      if ( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) {
         clatm1(MODER, CONDR, 0, IDIST, ISEED, DR, N, INFO );
         if ( INFO.NE.0 ) {
            INFO = 4
            RETURN
         }
      }

      // 3)     Generate IWORK if pivoting

      if ( IPVTNG.GT.0 ) {
         for (I = 1; I <= NPVTS; I++) { // 70
            IWORK( I ) = I
   70    CONTINUE
         if ( FULBND ) {
            for (I = 1; I <= NPVTS; I++) { // 80
               K = IPIVOT( I )
               J = IWORK( I )
               IWORK( I ) = IWORK( K )
               IWORK( K ) = J
   80       CONTINUE
         } else {
            DO 90 I = NPVTS, 1, -1
               K = IPIVOT( I )
               J = IWORK( I )
               IWORK( I ) = IWORK( K )
               IWORK( K ) = J
   90       CONTINUE
         }
      }

      // 4)      Generate matrices for each kind of PACKing
              // Always sweep matrix columnwise (if symmetric, upper
              // half only) so that matrix generated does not depend
              // on PACK

      if ( FULBND ) {

         // Use CLATM3 so matrices generated with differing PIVOTing only
         // differ only in the order of their rows and/or columns.

         if ( IPACK.EQ.0 ) {
            if ( ISYM.EQ.0 ) {
               for (J = 1; J <= N; J++) { // 110
                  for (I = 1; I <= J; I++) { // 100
                     CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = CTEMP
                     A( JSUB, ISUB ) = CONJG( CTEMP )
  100             CONTINUE
  110          CONTINUE
            } else if ( ISYM.EQ.1 ) {
               for (J = 1; J <= N; J++) { // 130
                  for (I = 1; I <= M; I++) { // 120
                     CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = CTEMP
  120             CONTINUE
  130          CONTINUE
            } else if ( ISYM.EQ.2 ) {
               for (J = 1; J <= N; J++) { // 150
                  for (I = 1; I <= J; I++) { // 140
                     CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = CTEMP
                     A( JSUB, ISUB ) = CTEMP
  140             CONTINUE
  150          CONTINUE
            }

         } else if ( IPACK.EQ.1 ) {

            for (J = 1; J <= N; J++) { // 170
               for (I = 1; I <= J; I++) { // 160
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  if ( MXSUB.EQ.ISUB .AND. ISYM.EQ.0 ) {
                     A( MNSUB, MXSUB ) = CONJG( CTEMP )
                  } else {
                     A( MNSUB, MXSUB ) = CTEMP
                  }
                  IF( MNSUB.NE.MXSUB ) A( MXSUB, MNSUB ) = CZERO
  160          CONTINUE
  170       CONTINUE

         } else if ( IPACK.EQ.2 ) {

            for (J = 1; J <= N; J++) { // 190
               for (I = 1; I <= J; I++) { // 180
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  if ( MXSUB.EQ.JSUB .AND. ISYM.EQ.0 ) {
                     A( MXSUB, MNSUB ) = CONJG( CTEMP )
                  } else {
                     A( MXSUB, MNSUB ) = CTEMP
                  }
                  IF( MNSUB.NE.MXSUB ) A( MNSUB, MXSUB ) = CZERO
  180          CONTINUE
  190       CONTINUE

         } else if ( IPACK.EQ.3 ) {

            for (J = 1; J <= N; J++) { // 210
               for (I = 1; I <= J; I++) { // 200
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )

                  // Compute K = location of (ISUB,JSUB) entry in packed
                  // array

                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  K = MXSUB*( MXSUB-1 ) / 2 + MNSUB

                  // Convert K to (IISUB,JJSUB) location

                  JJSUB = ( K-1 ) / LDA + 1
                  IISUB = K - LDA*( JJSUB-1 )

                  if ( MXSUB.EQ.ISUB .AND. ISYM.EQ.0 ) {
                     A( IISUB, JJSUB ) = CONJG( CTEMP )
                  } else {
                     A( IISUB, JJSUB ) = CTEMP
                  }
  200          CONTINUE
  210       CONTINUE

         } else if ( IPACK.EQ.4 ) {

            for (J = 1; J <= N; J++) { // 230
               for (I = 1; I <= J; I++) { // 220
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )

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

                  if ( MXSUB.EQ.JSUB .AND. ISYM.EQ.0 ) {
                     A( IISUB, JJSUB ) = CONJG( CTEMP )
                  } else {
                     A( IISUB, JJSUB ) = CTEMP
                  }
  220          CONTINUE
  230       CONTINUE

         } else if ( IPACK.EQ.5 ) {

            for (J = 1; J <= N; J++) { // 250
               DO 240 I = J - KUU, J
                  if ( I.LT.1 ) {
                     A( J-I+1, I+N ) = CZERO
                  } else {
                     CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     MNSUB = MIN( ISUB, JSUB )
                     MXSUB = MAX( ISUB, JSUB )
                     if ( MXSUB.EQ.JSUB .AND. ISYM.EQ.0 ) {
                        A( MXSUB-MNSUB+1, MNSUB ) = CONJG( CTEMP )
                     } else {
                        A( MXSUB-MNSUB+1, MNSUB ) = CTEMP
                     }
                  }
  240          CONTINUE
  250       CONTINUE

         } else if ( IPACK.EQ.6 ) {

            for (J = 1; J <= N; J++) { // 270
               DO 260 I = J - KUU, J
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  if ( MXSUB.EQ.ISUB .AND. ISYM.EQ.0 ) {
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = CONJG( CTEMP )
                  } else {
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = CTEMP
                  }
  260          CONTINUE
  270       CONTINUE

         } else if ( IPACK.EQ.7 ) {

            if ( ISYM.NE.1 ) {
               for (J = 1; J <= N; J++) { // 290
                  DO 280 I = J - KUU, J
                     CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     MNSUB = MIN( ISUB, JSUB )
                     MXSUB = MAX( ISUB, JSUB )
                     IF( I.LT.1 ) A( J-I+1+KUU, I+N ) = CZERO
                     if ( MXSUB.EQ.ISUB .AND. ISYM.EQ.0 ) {
                        A( MNSUB-MXSUB+KUU+1, MXSUB ) = CONJG( CTEMP )
                     } else {
                        A( MNSUB-MXSUB+KUU+1, MXSUB ) = CTEMP
                     }
                     if ( I.GE.1 .AND. MNSUB.NE.MXSUB ) {
                        if ( MNSUB.EQ.ISUB .AND. ISYM.EQ.0 ) {
                           A( MXSUB-MNSUB+1+KUU, MNSUB ) = CONJG( CTEMP )
                        } else {
                           A( MXSUB-MNSUB+1+KUU, MNSUB ) = CTEMP
                        }
                     }
  280             CONTINUE
  290          CONTINUE
            } else if ( ISYM.EQ.1 ) {
               for (J = 1; J <= N; J++) { // 310
                  DO 300 I = J - KUU, J + KLL
                     CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB-JSUB+KUU+1, JSUB ) = CTEMP
  300             CONTINUE
  310          CONTINUE
            }

         }

      } else {

         // Use CLATM2

         if ( IPACK.EQ.0 ) {
            if ( ISYM.EQ.0 ) {
               for (J = 1; J <= N; J++) { // 330
                  for (I = 1; I <= J; I++) { // 320
                     A( I, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( J, I ) = CONJG( A( I, J ) )
  320             CONTINUE
  330          CONTINUE
            } else if ( ISYM.EQ.1 ) {
               for (J = 1; J <= N; J++) { // 350
                  for (I = 1; I <= M; I++) { // 340
                     A( I, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  340             CONTINUE
  350          CONTINUE
            } else if ( ISYM.EQ.2 ) {
               for (J = 1; J <= N; J++) { // 370
                  for (I = 1; I <= J; I++) { // 360
                     A( I, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( J, I ) = A( I, J )
  360             CONTINUE
  370          CONTINUE
            }

         } else if ( IPACK.EQ.1 ) {

            for (J = 1; J <= N; J++) { // 390
               for (I = 1; I <= J; I++) { // 380
                  A( I, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )                   IF( I.NE.J ) A( J, I ) = CZERO
  380          CONTINUE
  390       CONTINUE

         } else if ( IPACK.EQ.2 ) {

            for (J = 1; J <= N; J++) { // 410
               for (I = 1; I <= J; I++) { // 400
                  if ( ISYM.EQ.0 ) {
                     A( J, I ) = CONJG( CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE ) )
                  } else {
                     A( J, I ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  }
                  IF( I.NE.J ) A( I, J ) = CZERO
  400          CONTINUE
  410       CONTINUE

         } else if ( IPACK.EQ.3 ) {

            ISUB = 0
            JSUB = 1
            for (J = 1; J <= N; J++) { // 430
               for (I = 1; I <= J; I++) { // 420
                  ISUB = ISUB + 1
                  if ( ISUB.GT.LDA ) {
                     ISUB = 1
                     JSUB = JSUB + 1
                  }
                  A( ISUB, JSUB ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  420          CONTINUE
  430       CONTINUE

         } else if ( IPACK.EQ.4 ) {

            if ( ISYM.EQ.0 .OR. ISYM.EQ.2 ) {
               for (J = 1; J <= N; J++) { // 450
                  for (I = 1; I <= J; I++) { // 440

                     // Compute K = location of (I,J) entry in packed array

                     if ( I.EQ.1 ) {
                        K = J
                     } else {
                        K = N*( N+1 ) / 2 - ( N-I+1 )*( N-I+2 ) / 2 + J - I + 1
                     }

                     // Convert K to (ISUB,JSUB) location

                     JSUB = ( K-1 ) / LDA + 1
                     ISUB = K - LDA*( JSUB-1 )

                     A( ISUB, JSUB ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     IF( ISYM.EQ.0 ) A( ISUB, JSUB ) = CONJG( A( ISUB, JSUB ) )
  440             CONTINUE
  450          CONTINUE
            } else {
               ISUB = 0
               JSUB = 1
               for (J = 1; J <= N; J++) { // 470
                  for (I = J; I <= M; I++) { // 460
                     ISUB = ISUB + 1
                     if ( ISUB.GT.LDA ) {
                        ISUB = 1
                        JSUB = JSUB + 1
                     }
                     A( ISUB, JSUB ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  460             CONTINUE
  470          CONTINUE
            }

         } else if ( IPACK.EQ.5 ) {

            for (J = 1; J <= N; J++) { // 490
               DO 480 I = J - KUU, J
                  if ( I.LT.1 ) {
                     A( J-I+1, I+N ) = CZERO
                  } else {
                     if ( ISYM.EQ.0 ) {
                        A( J-I+1, I ) = CONJG( CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE ) )
                     } else {
                        A( J-I+1, I ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     }
                  }
  480          CONTINUE
  490       CONTINUE

         } else if ( IPACK.EQ.6 ) {

            for (J = 1; J <= N; J++) { // 510
               DO 500 I = J - KUU, J
                  A( I-J+KUU+1, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  500          CONTINUE
  510       CONTINUE

         } else if ( IPACK.EQ.7 ) {

            if ( ISYM.NE.1 ) {
               for (J = 1; J <= N; J++) { // 530
                  DO 520 I = J - KUU, J
                     A( I-J+KUU+1, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     IF( I.LT.1 ) A( J-I+1+KUU, I+N ) = CZERO
                     if ( I.GE.1 .AND. I.NE.J ) {
                        if ( ISYM.EQ.0 ) {
                           A( J-I+1+KUU, I ) = CONJG( A( I-J+KUU+1, J ) )
                        } else {
                           A( J-I+1+KUU, I ) = A( I-J+KUU+1, J )
                        }
                     }
  520             CONTINUE
  530          CONTINUE
            } else if ( ISYM.EQ.1 ) {
               for (J = 1; J <= N; J++) { // 550
                  DO 540 I = J - KUU, J + KLL
                     A( I-J+KUU+1, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  540             CONTINUE
  550          CONTINUE
            }

         }

      }

      // 5)      Scaling the norm

      if ( IPACK.EQ.0 ) {
         ONORM = CLANGE( 'M', M, N, A, LDA, TEMPA )
      } else if ( IPACK.EQ.1 ) {
         ONORM = CLANSY( 'M', 'U', N, A, LDA, TEMPA )
      } else if ( IPACK.EQ.2 ) {
         ONORM = CLANSY( 'M', 'L', N, A, LDA, TEMPA )
      } else if ( IPACK.EQ.3 ) {
         ONORM = CLANSP( 'M', 'U', N, A, TEMPA )
      } else if ( IPACK.EQ.4 ) {
         ONORM = CLANSP( 'M', 'L', N, A, TEMPA )
      } else if ( IPACK.EQ.5 ) {
         ONORM = CLANSB( 'M', 'L', N, KLL, A, LDA, TEMPA )
      } else if ( IPACK.EQ.6 ) {
         ONORM = CLANSB( 'M', 'U', N, KUU, A, LDA, TEMPA )
      } else if ( IPACK.EQ.7 ) {
         ONORM = CLANGB( 'M', N, KLL, KUU, A, LDA, TEMPA )
      }

      if ( ANORM.GE.ZERO ) {

         if ( ANORM.GT.ZERO .AND. ONORM.EQ.ZERO ) {

            // Desired scaling impossible

            INFO = 5
            RETURN

         } else if ( ( ANORM.GT.ONE .AND. ONORM.LT.ONE ) .OR. ( ANORM.LT.ONE .AND. ONORM.GT.ONE ) ) {

            // Scale carefully to avoid over / underflow

            if ( IPACK.LE.2 ) {
               for (J = 1; J <= N; J++) { // 560
                  csscal(M, ONE / ONORM, A( 1, J ), 1 );
                  csscal(M, ANORM, A( 1, J ), 1 );
  560          CONTINUE

            } else if ( IPACK.EQ.3 .OR. IPACK.EQ.4 ) {

               csscal(N*( N+1 ) / 2, ONE / ONORM, A, 1 );
               csscal(N*( N+1 ) / 2, ANORM, A, 1 );

            } else if ( IPACK.GE.5 ) {

               for (J = 1; J <= N; J++) { // 570
                  csscal(KLL+KUU+1, ONE / ONORM, A( 1, J ), 1 );
                  csscal(KLL+KUU+1, ANORM, A( 1, J ), 1 );
  570          CONTINUE

            }

         } else {

            // Scale straightforwardly

            if ( IPACK.LE.2 ) {
               for (J = 1; J <= N; J++) { // 580
                  csscal(M, ANORM / ONORM, A( 1, J ), 1 );
  580          CONTINUE

            } else if ( IPACK.EQ.3 .OR. IPACK.EQ.4 ) {

               csscal(N*( N+1 ) / 2, ANORM / ONORM, A, 1 );

            } else if ( IPACK.GE.5 ) {

               for (J = 1; J <= N; J++) { // 590
                  csscal(KLL+KUU+1, ANORM / ONORM, A( 1, J ), 1 );
  590          CONTINUE
            }

         }

      }

      // End of CLATMR

      }
