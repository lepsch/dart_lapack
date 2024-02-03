      SUBROUTINE ZLATMR( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RSIGN, GRADE, DL, MODEL, CONDL, DR, MODER, CONDR, PIVTNG, IPIVOT, KL, KU, SPARSE, ANORM, PACK, A, LDA, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIST, GRADE, PACK, PIVTNG, RSIGN, SYM;
      int                INFO, KL, KU, LDA, M, MODE, MODEL, MODER, N;
      double             ANORM, COND, CONDL, CONDR, SPARSE;
      COMPLEX*16         DMAX
      // ..
      // .. Array Arguments ..
      int                IPIVOT( * ), ISEED( 4 ), IWORK( * );
      COMPLEX*16         A( LDA, * ), D( * ), DL( * ), DR( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      double             ONE;
      const              ONE = 1.0D0 ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D0, 0.0D0 ) ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D0, 0.0D0 ) ;
      // ..
      // .. Local Scalars ..
      bool               BADPVT, DZERO, FULBND;
      int                I, IDIST, IGRADE, IISUB, IPACK, IPVTNG, IRSIGN, ISUB, ISYM, J, JJSUB, JSUB, K, KLL, KUU, MNMIN, MNSUB, MXSUB, NPVTS;
      double             ONORM, TEMP;
      COMPLEX*16         CALPHA, CTEMP
      // ..
      // .. Local Arrays ..
      double             TEMPA( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             ZLANGB, ZLANGE, ZLANSB, ZLANSP, ZLANSY;
      COMPLEX*16         ZLATM2, ZLATM3
      // EXTERNAL LSAME, ZLANGB, ZLANGE, ZLANSB, ZLANSP, ZLANSY, ZLATM2, ZLATM3
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZLATM1
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, MAX, MIN, MOD
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
         DO 10 I = 1, M
            IF( DL( I ).EQ.CZERO ) DZERO = .TRUE.
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
         xerbla('ZLATMR', -INFO );
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

      zlatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, INFO );
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
         if ( TEMP.EQ.ZERO .AND. DMAX.NE.CZERO ) {
            INFO = 2
            RETURN
         }
         if ( TEMP.NE.ZERO ) {
            CALPHA = DMAX / TEMP
         } else {
            CALPHA = CONE
         }
         DO 50 I = 1, MNMIN
            D( I ) = CALPHA*D( I )
   50    CONTINUE

      }

      // If matrix Hermitian, make D real

      if ( ISYM.EQ.0 ) {
         DO 60 I = 1, MNMIN
            D( I ) = DBLE( D( I ) )
   60    CONTINUE
      }

      // Compute DL if grading set

      if ( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ. 5 .OR. IGRADE.EQ.6 ) {
         zlatm1(MODEL, CONDL, 0, IDIST, ISEED, DL, M, INFO );
         if ( INFO.NE.0 ) {
            INFO = 3
            RETURN
         }
      }

      // Compute DR if grading set

      if ( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) {
         zlatm1(MODER, CONDR, 0, IDIST, ISEED, DR, N, INFO );
         if ( INFO.NE.0 ) {
            INFO = 4
            RETURN
         }
      }

      // 3)     Generate IWORK if pivoting

      if ( IPVTNG.GT.0 ) {
         DO 70 I = 1, NPVTS
            IWORK( I ) = I
   70    CONTINUE
         if ( FULBND ) {
            DO 80 I = 1, NPVTS
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

         // Use ZLATM3 so matrices generated with differing PIVOTing only
         // differ only in the order of their rows and/or columns.

         if ( IPACK.EQ.0 ) {
            if ( ISYM.EQ.0 ) {
               DO 110 J = 1, N
                  DO 100 I = 1, J
                     CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = CTEMP
                     A( JSUB, ISUB ) = DCONJG( CTEMP )
  100             CONTINUE
  110          CONTINUE
            } else if ( ISYM.EQ.1 ) {
               DO 130 J = 1, N
                  DO 120 I = 1, M
                     CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = CTEMP
  120             CONTINUE
  130          CONTINUE
            } else if ( ISYM.EQ.2 ) {
               DO 150 J = 1, N
                  DO 140 I = 1, J
                     CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = CTEMP
                     A( JSUB, ISUB ) = CTEMP
  140             CONTINUE
  150          CONTINUE
            }

         } else if ( IPACK.EQ.1 ) {

            DO 170 J = 1, N
               DO 160 I = 1, J
                  CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  if ( MXSUB.EQ.ISUB .AND. ISYM.EQ.0 ) {
                     A( MNSUB, MXSUB ) = DCONJG( CTEMP )
                  } else {
                     A( MNSUB, MXSUB ) = CTEMP
                  }
                  IF( MNSUB.NE.MXSUB ) A( MXSUB, MNSUB ) = CZERO
  160          CONTINUE
  170       CONTINUE

         } else if ( IPACK.EQ.2 ) {

            DO 190 J = 1, N
               DO 180 I = 1, J
                  CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  if ( MXSUB.EQ.JSUB .AND. ISYM.EQ.0 ) {
                     A( MXSUB, MNSUB ) = DCONJG( CTEMP )
                  } else {
                     A( MXSUB, MNSUB ) = CTEMP
                  }
                  IF( MNSUB.NE.MXSUB ) A( MNSUB, MXSUB ) = CZERO
  180          CONTINUE
  190       CONTINUE

         } else if ( IPACK.EQ.3 ) {

            DO 210 J = 1, N
               DO 200 I = 1, J
                  CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )

                  // Compute K = location of (ISUB,JSUB) entry in packed
                  // array

                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  K = MXSUB*( MXSUB-1 ) / 2 + MNSUB

                  // Convert K to (IISUB,JJSUB) location

                  JJSUB = ( K-1 ) / LDA + 1
                  IISUB = K - LDA*( JJSUB-1 )

                  if ( MXSUB.EQ.ISUB .AND. ISYM.EQ.0 ) {
                     A( IISUB, JJSUB ) = DCONJG( CTEMP )
                  } else {
                     A( IISUB, JJSUB ) = CTEMP
                  }
  200          CONTINUE
  210       CONTINUE

         } else if ( IPACK.EQ.4 ) {

            DO 230 J = 1, N
               DO 220 I = 1, J
                  CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )

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
                     A( IISUB, JJSUB ) = DCONJG( CTEMP )
                  } else {
                     A( IISUB, JJSUB ) = CTEMP
                  }
  220          CONTINUE
  230       CONTINUE

         } else if ( IPACK.EQ.5 ) {

            DO 250 J = 1, N
               DO 240 I = J - KUU, J
                  if ( I.LT.1 ) {
                     A( J-I+1, I+N ) = CZERO
                  } else {
                     CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     MNSUB = MIN( ISUB, JSUB )
                     MXSUB = MAX( ISUB, JSUB )
                     if ( MXSUB.EQ.JSUB .AND. ISYM.EQ.0 ) {
                        A( MXSUB-MNSUB+1, MNSUB ) = DCONJG( CTEMP )
                     } else {
                        A( MXSUB-MNSUB+1, MNSUB ) = CTEMP
                     }
                  }
  240          CONTINUE
  250       CONTINUE

         } else if ( IPACK.EQ.6 ) {

            DO 270 J = 1, N
               DO 260 I = J - KUU, J
                  CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  if ( MXSUB.EQ.ISUB .AND. ISYM.EQ.0 ) {
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = DCONJG( CTEMP )
                  } else {
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = CTEMP
                  }
  260          CONTINUE
  270       CONTINUE

         } else if ( IPACK.EQ.7 ) {

            if ( ISYM.NE.1 ) {
               DO 290 J = 1, N
                  DO 280 I = J - KUU, J
                     CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     MNSUB = MIN( ISUB, JSUB )
                     MXSUB = MAX( ISUB, JSUB )
                     IF( I.LT.1 ) A( J-I+1+KUU, I+N ) = CZERO
                     if ( MXSUB.EQ.ISUB .AND. ISYM.EQ.0 ) {
                        A( MNSUB-MXSUB+KUU+1, MXSUB ) = DCONJG( CTEMP )
                     } else {
                        A( MNSUB-MXSUB+KUU+1, MXSUB ) = CTEMP
                     }
                     if ( I.GE.1 .AND. MNSUB.NE.MXSUB ) {
                        if ( MNSUB.EQ.ISUB .AND. ISYM.EQ.0 ) {
                           A( MXSUB-MNSUB+1+KUU, MNSUB ) = DCONJG( CTEMP )
                        } else {
                           A( MXSUB-MNSUB+1+KUU, MNSUB ) = CTEMP
                        }
                     }
  280             CONTINUE
  290          CONTINUE
            } else if ( ISYM.EQ.1 ) {
               DO 310 J = 1, N
                  DO 300 I = J - KUU, J + KLL
                     CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB-JSUB+KUU+1, JSUB ) = CTEMP
  300             CONTINUE
  310          CONTINUE
            }

         }

      } else {

         // Use ZLATM2

         if ( IPACK.EQ.0 ) {
            if ( ISYM.EQ.0 ) {
               DO 330 J = 1, N
                  DO 320 I = 1, J
                     A( I, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( J, I ) = DCONJG( A( I, J ) )
  320             CONTINUE
  330          CONTINUE
            } else if ( ISYM.EQ.1 ) {
               DO 350 J = 1, N
                  DO 340 I = 1, M
                     A( I, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  340             CONTINUE
  350          CONTINUE
            } else if ( ISYM.EQ.2 ) {
               DO 370 J = 1, N
                  DO 360 I = 1, J
                     A( I, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( J, I ) = A( I, J )
  360             CONTINUE
  370          CONTINUE
            }

         } else if ( IPACK.EQ.1 ) {

            DO 390 J = 1, N
               DO 380 I = 1, J
                  A( I, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )                   IF( I.NE.J ) A( J, I ) = CZERO
  380          CONTINUE
  390       CONTINUE

         } else if ( IPACK.EQ.2 ) {

            DO 410 J = 1, N
               DO 400 I = 1, J
                  if ( ISYM.EQ.0 ) {
                     A( J, I ) = DCONJG( ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE ) )
                  } else {
                     A( J, I ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  }
                  IF( I.NE.J ) A( I, J ) = CZERO
  400          CONTINUE
  410       CONTINUE

         } else if ( IPACK.EQ.3 ) {

            ISUB = 0
            JSUB = 1
            DO 430 J = 1, N
               DO 420 I = 1, J
                  ISUB = ISUB + 1
                  if ( ISUB.GT.LDA ) {
                     ISUB = 1
                     JSUB = JSUB + 1
                  }
                  A( ISUB, JSUB ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  420          CONTINUE
  430       CONTINUE

         } else if ( IPACK.EQ.4 ) {

            if ( ISYM.EQ.0 .OR. ISYM.EQ.2 ) {
               DO 450 J = 1, N
                  DO 440 I = 1, J

                     // Compute K = location of (I,J) entry in packed array

                     if ( I.EQ.1 ) {
                        K = J
                     } else {
                        K = N*( N+1 ) / 2 - ( N-I+1 )*( N-I+2 ) / 2 + J - I + 1
                     }

                     // Convert K to (ISUB,JSUB) location

                     JSUB = ( K-1 ) / LDA + 1
                     ISUB = K - LDA*( JSUB-1 )

                     A( ISUB, JSUB ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     IF( ISYM.EQ.0 ) A( ISUB, JSUB ) = DCONJG( A( ISUB, JSUB ) )
  440             CONTINUE
  450          CONTINUE
            } else {
               ISUB = 0
               JSUB = 1
               DO 470 J = 1, N
                  DO 460 I = J, M
                     ISUB = ISUB + 1
                     if ( ISUB.GT.LDA ) {
                        ISUB = 1
                        JSUB = JSUB + 1
                     }
                     A( ISUB, JSUB ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  460             CONTINUE
  470          CONTINUE
            }

         } else if ( IPACK.EQ.5 ) {

            DO 490 J = 1, N
               DO 480 I = J - KUU, J
                  if ( I.LT.1 ) {
                     A( J-I+1, I+N ) = CZERO
                  } else {
                     if ( ISYM.EQ.0 ) {
                        A( J-I+1, I ) = DCONJG( ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE ) )
                     } else {
                        A( J-I+1, I ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     }
                  }
  480          CONTINUE
  490       CONTINUE

         } else if ( IPACK.EQ.6 ) {

            DO 510 J = 1, N
               DO 500 I = J - KUU, J
                  A( I-J+KUU+1, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  500          CONTINUE
  510       CONTINUE

         } else if ( IPACK.EQ.7 ) {

            if ( ISYM.NE.1 ) {
               DO 530 J = 1, N
                  DO 520 I = J - KUU, J
                     A( I-J+KUU+1, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     IF( I.LT.1 ) A( J-I+1+KUU, I+N ) = CZERO
                     if ( I.GE.1 .AND. I.NE.J ) {
                        if ( ISYM.EQ.0 ) {
                           A( J-I+1+KUU, I ) = DCONJG( A( I-J+KUU+1, J ) )
                        } else {
                           A( J-I+1+KUU, I ) = A( I-J+KUU+1, J )
                        }
                     }
  520             CONTINUE
  530          CONTINUE
            } else if ( ISYM.EQ.1 ) {
               DO 550 J = 1, N
                  DO 540 I = J - KUU, J + KLL
                     A( I-J+KUU+1, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  540             CONTINUE
  550          CONTINUE
            }

         }

      }

      // 5)      Scaling the norm

      if ( IPACK.EQ.0 ) {
         ONORM = ZLANGE( 'M', M, N, A, LDA, TEMPA )
      } else if ( IPACK.EQ.1 ) {
         ONORM = ZLANSY( 'M', 'U', N, A, LDA, TEMPA )
      } else if ( IPACK.EQ.2 ) {
         ONORM = ZLANSY( 'M', 'L', N, A, LDA, TEMPA )
      } else if ( IPACK.EQ.3 ) {
         ONORM = ZLANSP( 'M', 'U', N, A, TEMPA )
      } else if ( IPACK.EQ.4 ) {
         ONORM = ZLANSP( 'M', 'L', N, A, TEMPA )
      } else if ( IPACK.EQ.5 ) {
         ONORM = ZLANSB( 'M', 'L', N, KLL, A, LDA, TEMPA )
      } else if ( IPACK.EQ.6 ) {
         ONORM = ZLANSB( 'M', 'U', N, KUU, A, LDA, TEMPA )
      } else if ( IPACK.EQ.7 ) {
         ONORM = ZLANGB( 'M', N, KLL, KUU, A, LDA, TEMPA )
      }

      if ( ANORM.GE.ZERO ) {

         if ( ANORM.GT.ZERO .AND. ONORM.EQ.ZERO ) {

            // Desired scaling impossible

            INFO = 5
            RETURN

         } else if ( ( ANORM.GT.ONE .AND. ONORM.LT.ONE ) .OR. ( ANORM.LT.ONE .AND. ONORM.GT.ONE ) ) {

            // Scale carefully to avoid over / underflow

            if ( IPACK.LE.2 ) {
               DO 560 J = 1, N
                  zdscal(M, ONE / ONORM, A( 1, J ), 1 );
                  zdscal(M, ANORM, A( 1, J ), 1 );
  560          CONTINUE

            } else if ( IPACK.EQ.3 .OR. IPACK.EQ.4 ) {

               zdscal(N*( N+1 ) / 2, ONE / ONORM, A, 1 );
               zdscal(N*( N+1 ) / 2, ANORM, A, 1 );

            } else if ( IPACK.GE.5 ) {

               DO 570 J = 1, N
                  zdscal(KLL+KUU+1, ONE / ONORM, A( 1, J ), 1 );
                  zdscal(KLL+KUU+1, ANORM, A( 1, J ), 1 );
  570          CONTINUE

            }

         } else {

            // Scale straightforwardly

            if ( IPACK.LE.2 ) {
               DO 580 J = 1, N
                  zdscal(M, ANORM / ONORM, A( 1, J ), 1 );
  580          CONTINUE

            } else if ( IPACK.EQ.3 .OR. IPACK.EQ.4 ) {

               zdscal(N*( N+1 ) / 2, ANORM / ONORM, A, 1 );

            } else if ( IPACK.GE.5 ) {

               DO 590 J = 1, N
                  zdscal(KLL+KUU+1, ANORM / ONORM, A( 1, J ), 1 );
  590          CONTINUE
            }

         }

      }

      // End of ZLATMR

      }
