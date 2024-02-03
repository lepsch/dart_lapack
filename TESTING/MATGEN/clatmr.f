      void clatmr(M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RSIGN, GRADE, DL, MODEL, CONDL, DR, MODER, CONDR, PIVTNG, IPIVOT, KL, KU, SPARSE, ANORM, PACK, A, LDA, IWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIST, GRADE, PACK, PIVTNG, RSIGN, SYM;
      int                INFO, KL, KU, LDA, M, MODE, MODEL, MODER, N;
      REAL               ANORM, COND, CONDL, CONDR, SPARSE;
      COMPLEX            DMAX;
      // ..
      // .. Array Arguments ..
      int                IPIVOT( * ), ISEED( 4 ), IWORK( * );
      COMPLEX            A( LDA, * ), D( * ), DL( * ), DR( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      REAL               ONE;
      const              ONE = 1.0 ;
      COMPLEX            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      COMPLEX            CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               BADPVT, DZERO, FULBND;
      int                I, IDIST, IGRADE, IISUB, IPACK, IPVTNG, IRSIGN, ISUB, ISYM, J, JJSUB, JSUB, K, KLL, KUU, MNMIN, MNSUB, MXSUB, NPVTS;
      REAL               ONORM, TEMP;
      COMPLEX            CALPHA, CTEMP;
      // ..
      // .. Local Arrays ..
      REAL               TEMPA( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGB, CLANGE, CLANSB, CLANSP, CLANSY;
      COMPLEX            CLATM2, CLATM3;
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

      INFO = 0;

      // Quick return if possible

      if (M == 0 || N == 0) return;

      // Decode DIST

      if ( LSAME( DIST, 'U' ) ) {
         IDIST = 1;
      } else if ( LSAME( DIST, 'S' ) ) {
         IDIST = 2;
      } else if ( LSAME( DIST, 'N' ) ) {
         IDIST = 3;
      } else if ( LSAME( DIST, 'D' ) ) {
         IDIST = 4;
      } else {
         IDIST = -1;
      }

      // Decode SYM

      if ( LSAME( SYM, 'H' ) ) {
         ISYM = 0;
      } else if ( LSAME( SYM, 'N' ) ) {
         ISYM = 1;
      } else if ( LSAME( SYM, 'S' ) ) {
         ISYM = 2;
      } else {
         ISYM = -1;
      }

      // Decode RSIGN

      if ( LSAME( RSIGN, 'F' ) ) {
         IRSIGN = 0;
      } else if ( LSAME( RSIGN, 'T' ) ) {
         IRSIGN = 1;
      } else {
         IRSIGN = -1;
      }

      // Decode PIVTNG

      if ( LSAME( PIVTNG, 'N' ) ) {
         IPVTNG = 0;
      } else if ( LSAME( PIVTNG, ' ' ) ) {
         IPVTNG = 0;
      } else if ( LSAME( PIVTNG, 'L' ) ) {
         IPVTNG = 1;
         NPVTS = M;
      } else if ( LSAME( PIVTNG, 'R' ) ) {
         IPVTNG = 2;
         NPVTS = N;
      } else if ( LSAME( PIVTNG, 'B' ) ) {
         IPVTNG = 3;
         NPVTS = min( N, M );
      } else if ( LSAME( PIVTNG, 'F' ) ) {
         IPVTNG = 3;
         NPVTS = min( N, M );
      } else {
         IPVTNG = -1;
      }

      // Decode GRADE

      if ( LSAME( GRADE, 'N' ) ) {
         IGRADE = 0;
      } else if ( LSAME( GRADE, 'L' ) ) {
         IGRADE = 1;
      } else if ( LSAME( GRADE, 'R' ) ) {
         IGRADE = 2;
      } else if ( LSAME( GRADE, 'B' ) ) {
         IGRADE = 3;
      } else if ( LSAME( GRADE, 'E' ) ) {
         IGRADE = 4;
      } else if ( LSAME( GRADE, 'H' ) ) {
         IGRADE = 5;
      } else if ( LSAME( GRADE, 'S' ) ) {
         IGRADE = 6;
      } else {
         IGRADE = -1;
      }

      // Decode PACK

      if ( LSAME( PACK, 'N' ) ) {
         IPACK = 0;
      } else if ( LSAME( PACK, 'U' ) ) {
         IPACK = 1;
      } else if ( LSAME( PACK, 'L' ) ) {
         IPACK = 2;
      } else if ( LSAME( PACK, 'C' ) ) {
         IPACK = 3;
      } else if ( LSAME( PACK, 'R' ) ) {
         IPACK = 4;
      } else if ( LSAME( PACK, 'B' ) ) {
         IPACK = 5;
      } else if ( LSAME( PACK, 'Q' ) ) {
         IPACK = 6;
      } else if ( LSAME( PACK, 'Z' ) ) {
         IPACK = 7;
      } else {
         IPACK = -1;
      }

      // Set certain internal parameters

      MNMIN = min( M, N );
      KLL = min( KL, M-1 );
      KUU = min( KU, N-1 );

      // If inv(DL) is used, check to see if DL has a zero entry.

      DZERO = false;
      if ( IGRADE == 4 && MODEL == 0 ) {
         for (I = 1; I <= M; I++) { // 10
            if( DL( I ) == CZERO ) DZERO = true;
         } // 10
      }

      // Check values in IPIVOT

      BADPVT = false;
      if ( IPVTNG > 0 ) {
         for (J = 1; J <= NPVTS; J++) { // 20
            if( IPIVOT( J ) <= 0 || IPIVOT( J ) > NPVTS ) BADPVT = true;
         } // 20
      }

      // Set INFO if an error

      if ( M < 0 ) {
         INFO = -1;
      } else if ( M != N && ( ISYM == 0 || ISYM == 2 ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( IDIST == -1 ) {
         INFO = -3;
      } else if ( ISYM == -1 ) {
         INFO = -5;
      } else if ( MODE < -6 || MODE > 6 ) {
         INFO = -7;
      } else if ( ( MODE != -6 && MODE != 0 && MODE != 6 ) && COND < ONE ) {
         INFO = -8;
      } else if ( ( MODE != -6 && MODE != 0 && MODE != 6 ) && IRSIGN == -1 ) {
         INFO = -10;
      } else if ( IGRADE == -1 || ( IGRADE == 4 && M != N ) || ( ( IGRADE == 1 || IGRADE == 2 || IGRADE == 3 || IGRADE == 4 || IGRADE == 6 ) && ISYM == 0 ) || ( ( IGRADE == 1 || IGRADE == 2 || IGRADE == 3 || IGRADE == 4 || IGRADE == 5 ) && ISYM == 2 ) ) {
         INFO = -11;
      } else if ( IGRADE == 4 && DZERO ) {
         INFO = -12;
      } else if ( ( IGRADE == 1 || IGRADE == 3 || IGRADE == 4 || IGRADE == 5 || IGRADE == 6 ) && ( MODEL < -6 || MODEL > 6 ) ) {
         INFO = -13;
      } else if ( ( IGRADE == 1 || IGRADE == 3 || IGRADE == 4 || IGRADE == 5 || IGRADE == 6 ) && ( MODEL != -6 && MODEL != 0 && MODEL != 6 ) && CONDL < ONE ) {
         INFO = -14;
      } else if ( ( IGRADE == 2 || IGRADE == 3 ) && ( MODER < -6 || MODER > 6 ) ) {
         INFO = -16;
      } else if ( ( IGRADE == 2 || IGRADE == 3 ) && ( MODER != -6 && MODER != 0 && MODER != 6 ) && CONDR < ONE ) {
         INFO = -17;
      } else if ( IPVTNG == -1 || ( IPVTNG == 3 && M != N ) || ( ( IPVTNG == 1 || IPVTNG == 2 ) && ( ISYM == 0 || ISYM == 2 ) ) ) {
         INFO = -18;
      } else if ( IPVTNG != 0 && BADPVT ) {
         INFO = -19;
      } else if ( KL < 0 ) {
         INFO = -20;
      } else if ( KU < 0 || ( ( ISYM == 0 || ISYM == 2 ) && KL != KU ) ) {
         INFO = -21;
      } else if ( SPARSE < ZERO || SPARSE > ONE ) {
         INFO = -22;
      } else if ( IPACK == -1 || ( ( IPACK == 1 || IPACK == 2 || IPACK == 5 || IPACK == 6 ) && ISYM == 1 ) || ( IPACK == 3 && ISYM == 1 && ( KL != 0 || M != N ) ) || ( IPACK == 4 && ISYM == 1 && ( KU != 0 || M != N ) ) ) {
         INFO = -24;
      } else if ( ( ( IPACK == 0 || IPACK == 1 || IPACK == 2 ) && LDA < max( 1, M ) ) || ( ( IPACK == 3 || IPACK == 4 ) && LDA < 1 ) || ( ( IPACK == 5 || IPACK == 6 ) && LDA < KUU+1 ) || ( IPACK == 7 && LDA < KLL+KUU+1 ) ) {
         INFO = -26;
      }

      if ( INFO != 0 ) {
         xerbla('CLATMR', -INFO );
         return;
      }

      // Decide if we can pivot consistently

      FULBND = false;
      if (KUU == N-1 && KLL == M-1) FULBND = true ;

      // Initialize random number generator

      for (I = 1; I <= 4; I++) { // 30
         ISEED( I ) = MOD( ( ISEED( I ) ).abs(), 4096 );
      } // 30

      ISEED( 4 ) = 2*( ISEED( 4 ) / 2 ) + 1;

      // 2)      Set up D, DL, and DR, if indicated.

              // Compute D according to COND and MODE

      clatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, INFO );
      if ( INFO != 0 ) {
         INFO = 1;
         return;
      }
      if ( MODE != 0 && MODE != -6 && MODE != 6 ) {

         // Scale by DMAX

         TEMP = ( D( 1 ) ).abs();
         for (I = 2; I <= MNMIN; I++) { // 40
            TEMP = max( TEMP, ( D( I ) ) ).abs();
         } // 40
         if ( TEMP == ZERO && DMAX != CZERO ) {
            INFO = 2;
            return;
         }
         if ( TEMP != ZERO ) {
            CALPHA = DMAX / TEMP;
         } else {
            CALPHA = CONE;
         }
         for (I = 1; I <= MNMIN; I++) { // 50
            D( I ) = CALPHA*D( I );
         } // 50

      }

      // If matrix Hermitian, make D real

      if ( ISYM == 0 ) {
         for (I = 1; I <= MNMIN; I++) { // 60
            D( I ) = REAL( D( I ) );
         } // 60
      }

      // Compute DL if grading set

      if ( IGRADE == 1 || IGRADE == 3 || IGRADE == 4 || IGRADE == 5 || IGRADE == 6 ) {
         clatm1(MODEL, CONDL, 0, IDIST, ISEED, DL, M, INFO );
         if ( INFO != 0 ) {
            INFO = 3;
            return;
         }
      }

      // Compute DR if grading set

      if ( IGRADE == 2 || IGRADE == 3 ) {
         clatm1(MODER, CONDR, 0, IDIST, ISEED, DR, N, INFO );
         if ( INFO != 0 ) {
            INFO = 4;
            return;
         }
      }

      // 3)     Generate IWORK if pivoting

      if ( IPVTNG > 0 ) {
         for (I = 1; I <= NPVTS; I++) { // 70
            IWORK( I ) = I;
         } // 70
         if ( FULBND ) {
            for (I = 1; I <= NPVTS; I++) { // 80
               K = IPIVOT( I );
               J = IWORK( I );
               IWORK( I ) = IWORK( K );
               IWORK( K ) = J;
            } // 80
         } else {
            DO 90 I = NPVTS, 1, -1;
               K = IPIVOT( I );
               J = IWORK( I );
               IWORK( I ) = IWORK( K );
               IWORK( K ) = J;
            } // 90
         }
      }

      // 4)      Generate matrices for each kind of PACKing
              // Always sweep matrix columnwise (if symmetric, upper
              // half only) so that matrix generated does not depend
              // on PACK

      if ( FULBND ) {

         // Use CLATM3 so matrices generated with differing PIVOTing only
         // differ only in the order of their rows and/or columns.

         if ( IPACK == 0 ) {
            if ( ISYM == 0 ) {
               for (J = 1; J <= N; J++) { // 110
                  for (I = 1; I <= J; I++) { // 100
                     CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                     A( ISUB, JSUB ) = CTEMP;
                     A( JSUB, ISUB ) = CONJG( CTEMP );
                  } // 100
               } // 110
            } else if ( ISYM == 1 ) {
               for (J = 1; J <= N; J++) { // 130
                  for (I = 1; I <= M; I++) { // 120
                     CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                     A( ISUB, JSUB ) = CTEMP;
                  } // 120
               } // 130
            } else if ( ISYM == 2 ) {
               for (J = 1; J <= N; J++) { // 150
                  for (I = 1; I <= J; I++) { // 140
                     CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                     A( ISUB, JSUB ) = CTEMP;
                     A( JSUB, ISUB ) = CTEMP;
                  } // 140
               } // 150
            }

         } else if ( IPACK == 1 ) {

            for (J = 1; J <= N; J++) { // 170
               for (I = 1; I <= J; I++) { // 160
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                  MNSUB = min( ISUB, JSUB );
                  MXSUB = max( ISUB, JSUB );
                  if ( MXSUB == ISUB && ISYM == 0 ) {
                     A( MNSUB, MXSUB ) = CONJG( CTEMP );
                  } else {
                     A( MNSUB, MXSUB ) = CTEMP;
                  }
                  if (MNSUB != MXSUB) A( MXSUB, MNSUB ) = CZERO;
               } // 160
            } // 170

         } else if ( IPACK == 2 ) {

            for (J = 1; J <= N; J++) { // 190
               for (I = 1; I <= J; I++) { // 180
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                  MNSUB = min( ISUB, JSUB );
                  MXSUB = max( ISUB, JSUB );
                  if ( MXSUB == JSUB && ISYM == 0 ) {
                     A( MXSUB, MNSUB ) = CONJG( CTEMP );
                  } else {
                     A( MXSUB, MNSUB ) = CTEMP;
                  }
                  if (MNSUB != MXSUB) A( MNSUB, MXSUB ) = CZERO;
               } // 180
            } // 190

         } else if ( IPACK == 3 ) {

            for (J = 1; J <= N; J++) { // 210
               for (I = 1; I <= J; I++) { // 200
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );

                  // Compute K = location of (ISUB,JSUB) entry in packed
                  // array

                  MNSUB = min( ISUB, JSUB );
                  MXSUB = max( ISUB, JSUB );
                  K = MXSUB*( MXSUB-1 ) / 2 + MNSUB;

                  // Convert K to (IISUB,JJSUB) location

                  JJSUB = ( K-1 ) / LDA + 1;
                  IISUB = K - LDA*( JJSUB-1 );

                  if ( MXSUB == ISUB && ISYM == 0 ) {
                     A( IISUB, JJSUB ) = CONJG( CTEMP );
                  } else {
                     A( IISUB, JJSUB ) = CTEMP;
                  }
               } // 200
            } // 210

         } else if ( IPACK == 4 ) {

            for (J = 1; J <= N; J++) { // 230
               for (I = 1; I <= J; I++) { // 220
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );

                  // Compute K = location of (I,J) entry in packed array

                  MNSUB = min( ISUB, JSUB );
                  MXSUB = max( ISUB, JSUB );
                  if ( MNSUB == 1 ) {
                     K = MXSUB;
                  } else {
                     K = N*( N+1 ) / 2 - ( N-MNSUB+1 )*( N-MNSUB+2 ) / 2 + MXSUB - MNSUB + 1;
                  }

                  // Convert K to (IISUB,JJSUB) location

                  JJSUB = ( K-1 ) / LDA + 1;
                  IISUB = K - LDA*( JJSUB-1 );

                  if ( MXSUB == JSUB && ISYM == 0 ) {
                     A( IISUB, JJSUB ) = CONJG( CTEMP );
                  } else {
                     A( IISUB, JJSUB ) = CTEMP;
                  }
               } // 220
            } // 230

         } else if ( IPACK == 5 ) {

            for (J = 1; J <= N; J++) { // 250
               for (I = J - KUU; I <= J; I++) { // 240
                  if ( I < 1 ) {
                     A( J-I+1, I+N ) = CZERO;
                  } else {
                     CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                     MNSUB = min( ISUB, JSUB );
                     MXSUB = max( ISUB, JSUB );
                     if ( MXSUB == JSUB && ISYM == 0 ) {
                        A( MXSUB-MNSUB+1, MNSUB ) = CONJG( CTEMP );
                     } else {
                        A( MXSUB-MNSUB+1, MNSUB ) = CTEMP;
                     }
                  }
               } // 240
            } // 250

         } else if ( IPACK == 6 ) {

            for (J = 1; J <= N; J++) { // 270
               for (I = J - KUU; I <= J; I++) { // 260
                  CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                  MNSUB = min( ISUB, JSUB );
                  MXSUB = max( ISUB, JSUB );
                  if ( MXSUB == ISUB && ISYM == 0 ) {
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = CONJG( CTEMP );
                  } else {
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = CTEMP;
                  }
               } // 260
            } // 270

         } else if ( IPACK == 7 ) {

            if ( ISYM != 1 ) {
               for (J = 1; J <= N; J++) { // 290
                  for (I = J - KUU; I <= J; I++) { // 280
                     CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                     MNSUB = min( ISUB, JSUB );
                     MXSUB = max( ISUB, JSUB );
                     if (I < 1) A( J-I+1+KUU, I+N ) = CZERO;
                     if ( MXSUB == ISUB && ISYM == 0 ) {
                        A( MNSUB-MXSUB+KUU+1, MXSUB ) = CONJG( CTEMP );
                     } else {
                        A( MNSUB-MXSUB+KUU+1, MXSUB ) = CTEMP;
                     }
                     if ( I >= 1 && MNSUB != MXSUB ) {
                        if ( MNSUB == ISUB && ISYM == 0 ) {
                           A( MXSUB-MNSUB+1+KUU, MNSUB ) = CONJG( CTEMP );
                        } else {
                           A( MXSUB-MNSUB+1+KUU, MNSUB ) = CTEMP;
                        }
                     }
                  } // 280
               } // 290
            } else if ( ISYM == 1 ) {
               for (J = 1; J <= N; J++) { // 310
                  for (I = J - KUU; I <= J + KLL; I++) { // 300
                     CTEMP = CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                     A( ISUB-JSUB+KUU+1, JSUB ) = CTEMP;
                  } // 300
               } // 310
            }

         }

      } else {

         // Use CLATM2

         if ( IPACK == 0 ) {
            if ( ISYM == 0 ) {
               for (J = 1; J <= N; J++) { // 330
                  for (I = 1; I <= J; I++) { // 320
                     A( I, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                     A( J, I ) = CONJG( A( I, J ) );
                  } // 320
               } // 330
            } else if ( ISYM == 1 ) {
               for (J = 1; J <= N; J++) { // 350
                  for (I = 1; I <= M; I++) { // 340
                     A( I, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                  } // 340
               } // 350
            } else if ( ISYM == 2 ) {
               for (J = 1; J <= N; J++) { // 370
                  for (I = 1; I <= J; I++) { // 360
                     A( I, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                     A( J, I ) = A( I, J );
                  } // 360
               } // 370
            }

         } else if ( IPACK == 1 ) {

            for (J = 1; J <= N; J++) { // 390
               for (I = 1; I <= J; I++) { // 380
                  A( I, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )                   IF( I != J ) A( J, I ) = CZERO;
               } // 380
            } // 390

         } else if ( IPACK == 2 ) {

            for (J = 1; J <= N; J++) { // 410
               for (I = 1; I <= J; I++) { // 400
                  if ( ISYM == 0 ) {
                     A( J, I ) = CONJG( CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE ) );
                  } else {
                     A( J, I ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                  }
                  if (I != J) A( I, J ) = CZERO;
               } // 400
            } // 410

         } else if ( IPACK == 3 ) {

            ISUB = 0;
            JSUB = 1;
            for (J = 1; J <= N; J++) { // 430
               for (I = 1; I <= J; I++) { // 420
                  ISUB = ISUB + 1;
                  if ( ISUB > LDA ) {
                     ISUB = 1;
                     JSUB = JSUB + 1;
                  }
                  A( ISUB, JSUB ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
               } // 420
            } // 430

         } else if ( IPACK == 4 ) {

            if ( ISYM == 0 || ISYM == 2 ) {
               for (J = 1; J <= N; J++) { // 450
                  for (I = 1; I <= J; I++) { // 440

                     // Compute K = location of (I,J) entry in packed array

                     if ( I == 1 ) {
                        K = J;
                     } else {
                        K = N*( N+1 ) / 2 - ( N-I+1 )*( N-I+2 ) / 2 + J - I + 1;
                     }

                     // Convert K to (ISUB,JSUB) location

                     JSUB = ( K-1 ) / LDA + 1;
                     ISUB = K - LDA*( JSUB-1 );

                     A( ISUB, JSUB ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                     if (ISYM == 0) A( ISUB, JSUB ) = CONJG( A( ISUB, JSUB ) );
                  } // 440
               } // 450
            } else {
               ISUB = 0;
               JSUB = 1;
               for (J = 1; J <= N; J++) { // 470
                  for (I = J; I <= M; I++) { // 460
                     ISUB = ISUB + 1;
                     if ( ISUB > LDA ) {
                        ISUB = 1;
                        JSUB = JSUB + 1;
                     }
                     A( ISUB, JSUB ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                  } // 460
               } // 470
            }

         } else if ( IPACK == 5 ) {

            for (J = 1; J <= N; J++) { // 490
               for (I = J - KUU; I <= J; I++) { // 480
                  if ( I < 1 ) {
                     A( J-I+1, I+N ) = CZERO;
                  } else {
                     if ( ISYM == 0 ) {
                        A( J-I+1, I ) = CONJG( CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE ) );
                     } else {
                        A( J-I+1, I ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                     }
                  }
               } // 480
            } // 490

         } else if ( IPACK == 6 ) {

            for (J = 1; J <= N; J++) { // 510
               for (I = J - KUU; I <= J; I++) { // 500
                  A( I-J+KUU+1, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
               } // 500
            } // 510

         } else if ( IPACK == 7 ) {

            if ( ISYM != 1 ) {
               for (J = 1; J <= N; J++) { // 530
                  for (I = J - KUU; I <= J; I++) { // 520
                     A( I-J+KUU+1, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                     if (I < 1) A( J-I+1+KUU, I+N ) = CZERO;
                     if ( I >= 1 && I != J ) {
                        if ( ISYM == 0 ) {
                           A( J-I+1+KUU, I ) = CONJG( A( I-J+KUU+1, J ) );
                        } else {
                           A( J-I+1+KUU, I ) = A( I-J+KUU+1, J );
                        }
                     }
                  } // 520
               } // 530
            } else if ( ISYM == 1 ) {
               for (J = 1; J <= N; J++) { // 550
                  for (I = J - KUU; I <= J + KLL; I++) { // 540
                     A( I-J+KUU+1, J ) = CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );
                  } // 540
               } // 550
            }

         }

      }

      // 5)      Scaling the norm

      if ( IPACK == 0 ) {
         ONORM = CLANGE( 'M', M, N, A, LDA, TEMPA );
      } else if ( IPACK == 1 ) {
         ONORM = CLANSY( 'M', 'U', N, A, LDA, TEMPA );
      } else if ( IPACK == 2 ) {
         ONORM = CLANSY( 'M', 'L', N, A, LDA, TEMPA );
      } else if ( IPACK == 3 ) {
         ONORM = CLANSP( 'M', 'U', N, A, TEMPA );
      } else if ( IPACK == 4 ) {
         ONORM = CLANSP( 'M', 'L', N, A, TEMPA );
      } else if ( IPACK == 5 ) {
         ONORM = CLANSB( 'M', 'L', N, KLL, A, LDA, TEMPA );
      } else if ( IPACK == 6 ) {
         ONORM = CLANSB( 'M', 'U', N, KUU, A, LDA, TEMPA );
      } else if ( IPACK == 7 ) {
         ONORM = CLANGB( 'M', N, KLL, KUU, A, LDA, TEMPA );
      }

      if ( ANORM >= ZERO ) {

         if ( ANORM > ZERO && ONORM == ZERO ) {

            // Desired scaling impossible

            INFO = 5;
            return;

         } else if ( ( ANORM > ONE && ONORM < ONE ) || ( ANORM < ONE && ONORM > ONE ) ) {

            // Scale carefully to avoid over / underflow

            if ( IPACK <= 2 ) {
               for (J = 1; J <= N; J++) { // 560
                  csscal(M, ONE / ONORM, A( 1, J ), 1 );
                  csscal(M, ANORM, A( 1, J ), 1 );
               } // 560

            } else if ( IPACK == 3 || IPACK == 4 ) {

               csscal(N*( N+1 ) / 2, ONE / ONORM, A, 1 );
               csscal(N*( N+1 ) / 2, ANORM, A, 1 );

            } else if ( IPACK >= 5 ) {

               for (J = 1; J <= N; J++) { // 570
                  csscal(KLL+KUU+1, ONE / ONORM, A( 1, J ), 1 );
                  csscal(KLL+KUU+1, ANORM, A( 1, J ), 1 );
               } // 570

            }

         } else {

            // Scale straightforwardly

            if ( IPACK <= 2 ) {
               for (J = 1; J <= N; J++) { // 580
                  csscal(M, ANORM / ONORM, A( 1, J ), 1 );
               } // 580

            } else if ( IPACK == 3 || IPACK == 4 ) {

               csscal(N*( N+1 ) / 2, ANORM / ONORM, A, 1 );

            } else if ( IPACK >= 5 ) {

               for (J = 1; J <= N; J++) { // 590
                  csscal(KLL+KUU+1, ANORM / ONORM, A( 1, J ), 1 );
               } // 590
            }

         }

      }
      }
