      void slatmt(M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RANK, KL, KU, PACK, A, LDA, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double               COND, DMAX;
      int                INFO, KL, KU, LDA, M, MODE, N, RANK;
      String             DIST, PACK, SYM;
      // ..
      // .. Array Arguments ..
      double               A( LDA, * ), D( * ), WORK( * );
      int                ISEED( 4 );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      double               ONE;
      const              ONE = 1.0 ;
      double               TWOPI;
      const      TWOPI = 6.28318530717958647692528676655900576839 ;
      // ..
      // .. Local Scalars ..
      double               ALPHA, ANGLE, C, DUMMY, EXTRA, S, TEMP;
      int                I, IC, ICOL, IDIST, IENDCH, IINFO, IL, ILDA, IOFFG, IOFFST, IPACK, IPACKG, IR, IR1, IR2, IROW, IRSIGN, ISKEW, ISYM, ISYMPK, J, JC, JCH, JKL, JKU, JR, K, LLB, MINLDA, MNMIN, MR, NC, UUB;
      bool               GIVENS, ILEXTR, ILTEMP, TOPDWN;
      // ..
      // .. External Functions ..
      //- REAL               SLARND;
      //- bool               lsame;
      // EXTERNAL SLARND, lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLATM7, SCOPY, SLAGGE, SLAGSY, SLAROT, SLARTG, SLASET, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, COS, MAX, MIN, MOD, REAL, SIN
      // ..
      // .. Executable Statements ..

      // 1)      Decode and Test the input parameters.
              // Initialize flags & seed.

      INFO = 0;

      // Quick return if possible

      if (M == 0 || N == 0) return;

      // Decode DIST

      if ( lsame( DIST, 'U' ) ) {
         IDIST = 1;
      } else if ( lsame( DIST, 'S' ) ) {
         IDIST = 2;
      } else if ( lsame( DIST, 'N' ) ) {
         IDIST = 3;
      } else {
         IDIST = -1;
      }

      // Decode SYM

      if ( lsame( SYM, 'N' ) ) {
         ISYM = 1;
         IRSIGN = 0;
      } else if ( lsame( SYM, 'P' ) ) {
         ISYM = 2;
         IRSIGN = 0;
      } else if ( lsame( SYM, 'S' ) ) {
         ISYM = 2;
         IRSIGN = 1;
      } else if ( lsame( SYM, 'H' ) ) {
         ISYM = 2;
         IRSIGN = 1;
      } else {
         ISYM = -1;
      }

      // Decode PACK

      ISYMPK = 0;
      if ( lsame( PACK, 'N' ) ) {
         IPACK = 0;
      } else if ( lsame( PACK, 'U' ) ) {
         IPACK = 1;
         ISYMPK = 1;
      } else if ( lsame( PACK, 'L' ) ) {
         IPACK = 2;
         ISYMPK = 1;
      } else if ( lsame( PACK, 'C' ) ) {
         IPACK = 3;
         ISYMPK = 2;
      } else if ( lsame( PACK, 'R' ) ) {
         IPACK = 4;
         ISYMPK = 3;
      } else if ( lsame( PACK, 'B' ) ) {
         IPACK = 5;
         ISYMPK = 3;
      } else if ( lsame( PACK, 'Q' ) ) {
         IPACK = 6;
         ISYMPK = 2;
      } else if ( lsame( PACK, 'Z' ) ) {
         IPACK = 7;
      } else {
         IPACK = -1;
      }

      // Set certain internal parameters

      MNMIN = min( M, N );
      LLB = min( KL, M-1 );
      UUB = min( KU, N-1 );
      MR = min( M, N+LLB );
      NC = min( N, M+UUB );

      if ( IPACK == 5 || IPACK == 6 ) {
         MINLDA = UUB + 1;
      } else if ( IPACK == 7 ) {
         MINLDA = LLB + UUB + 1;
      } else {
         MINLDA = M;
      }

      // Use Givens rotation method if bandwidth small enough,
      // or if LDA is too small to store the matrix unpacked.

      GIVENS = false;
      if ( ISYM == 1 ) {
         if( REAL( LLB+UUB ) < 0.3*double( max( 1, MR+NC ) ) ) GIVENS = true;
      } else {
         if (2*LLB < M) GIVENS = true ;
      }
      if (LDA < M && LDA >= MINLDA) GIVENS = true ;

      // Set INFO if an error

      if ( M < 0 ) {
         INFO = -1;
      } else if ( M != N && ISYM != 1 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( IDIST == -1 ) {
         INFO = -3;
      } else if ( ISYM == -1 ) {
         INFO = -5;
      } else if ( ( MODE ).abs() > 6 ) {
         INFO = -7;
      } else if ( ( MODE != 0 && ( MODE ).abs() != 6 ) && COND < ONE ) {
         INFO = -8;
      } else if ( KL < 0 ) {
         INFO = -10;
      } else if ( KU < 0 || ( ISYM != 1 && KL != KU ) ) {
         INFO = -11;
      } else if ( IPACK == -1 || ( ISYMPK == 1 && ISYM == 1 ) || ( ISYMPK == 2 && ISYM == 1 && KL > 0 ) || ( ISYMPK == 3 && ISYM == 1 && KU > 0 ) || ( ISYMPK != 0 && M != N ) ) {
         INFO = -12;
      } else if ( LDA < max( 1, MINLDA ) ) {
         INFO = -14;
      }

      if ( INFO != 0 ) {
         xerbla('SLATMT', -INFO );
         return;
      }

      // Initialize random number generator

      for (I = 1; I <= 4; I++) { // 100
         ISEED[I] = (( ISEED( I ) ).abs() % 4096);
      } // 100

      if( (ISEED( 4 ) % 2) != 1 ) ISEED( 4 ) = ISEED( 4 ) + 1;

      // 2)      Set up D  if indicated.

              // Compute D according to COND and MODE

      slatm7(MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, RANK, IINFO );
      if ( IINFO != 0 ) {
         INFO = 1;
         return;
      }

      // Choose Top-Down if D is (apparently) increasing,
      // Bottom-Up if D is (apparently) decreasing.

      if ( ( D( 1 ) ).abs() <= ( D( RANK ) ).abs() ) {
         TOPDWN = true;
      } else {
         TOPDWN = false;
      }

      if ( MODE != 0 && ( MODE ).abs() != 6 ) {

         // Scale by DMAX

         TEMP = ( D( 1 ) ).abs();
         for (I = 2; I <= RANK; I++) { // 110
            TEMP = max( TEMP, ( D( I ) ) ).abs();
         } // 110

         if ( TEMP > ZERO ) {
            ALPHA = DMAX / TEMP;
         } else {
            INFO = 2;
            return;
         }

         sscal(RANK, ALPHA, D, 1 );

      }

      // 3)      Generate Banded Matrix using Givens rotations.
              // Also the special case of UUB=LLB=0

                // Compute Addressing constants to cover all
                // storage formats.  Whether GE, SY, GB, or SB,
                // upper or lower triangle or both,
                // the (i,j)-th element is in
                // A( i - ISKEW*j + IOFFST, j )

      if ( IPACK > 4 ) {
         ILDA = LDA - 1;
         ISKEW = 1;
         if ( IPACK > 5 ) {
            IOFFST = UUB + 1;
         } else {
            IOFFST = 1;
         }
      } else {
         ILDA = LDA;
         ISKEW = 0;
         IOFFST = 0;
      }

      // IPACKG is the format that the matrix is generated in. If this is
      // different from IPACK, then the matrix must be repacked at the
      // end.  It also signals how to compute the norm, for scaling.

      IPACKG = 0;
      slaset('Full', LDA, N, ZERO, ZERO, A, LDA );

      // Diagonal Matrix -- We are done, unless it
      // is to be stored SP/PP/TP (PACK='R' or 'C')

      if ( LLB == 0 && UUB == 0 ) {
         scopy(MNMIN, D, 1, A( 1-ISKEW+IOFFST, 1 ), ILDA+1 );
         if (IPACK <= 2 || IPACK >= 5) IPACKG = IPACK;

      } else if ( GIVENS ) {

         // Check whether to use Givens rotations,
         // Householder transformations, or nothing.

         if ( ISYM == 1 ) {

            // Non-symmetric -- A = U D V

            if ( IPACK > 4 ) {
               IPACKG = IPACK;
            } else {
               IPACKG = 0;
            }

            scopy(MNMIN, D, 1, A( 1-ISKEW+IOFFST, 1 ), ILDA+1 );

            if ( TOPDWN ) {
               JKL = 0;
               for (JKU = 1; JKU <= UUB; JKU++) { // 140

                  // Transform from bandwidth JKL, JKU-1 to JKL, JKU

                  // Last row actually rotated is M
                  // Last column actually rotated is min( M+JKU, N )

                  for (JR = 1; JR <= min( M+JKU, N ) + JKL - 1; JR++) { // 130
                     EXTRA = ZERO;
                     ANGLE = TWOPI*SLARND( 1, ISEED );
                     C = COS( ANGLE );
                     S = SIN( ANGLE );
                     ICOL = max( 1, JR-JKL );
                     if ( JR < M ) {
                        IL = min( N, JR+JKU ) + 1 - ICOL;
                        slarot( true , JR > JKL, false , IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, EXTRA, DUMMY );
                     }

                     // Chase "EXTRA" back up

                     IR = JR;
                     IC = ICOL;
                     for (JCH = JR - JKL; -JKL - JKU < 0 ? JCH >= 1 : JCH <= 1; JCH += -JKL - JKU) { // 120
                        if ( IR < M ) {
                           slartg(A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, C, S, DUMMY );
                        }
                        IROW = max( 1, JCH-JKU );
                        IL = IR + 2 - IROW;
                        TEMP = ZERO;
                        ILTEMP = JCH > JKU;
                        slarot( false , ILTEMP, true , IL, C, -S, A( IROW-ISKEW*IC+IOFFST, IC ), ILDA, TEMP, EXTRA );
                        if ( ILTEMP ) {
                           slartg(A( IROW+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), TEMP, C, S, DUMMY );
                           ICOL = max( 1, JCH-JKU-JKL );
                           IL = IC + 2 - ICOL;
                           EXTRA = ZERO;
                           slarot( true , JCH > JKU+JKL, true , IL, C, -S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, TEMP );
                           IC = ICOL;
                           IR = IROW;
                        }
                     } // 120
                  } // 130
               } // 140

               JKU = UUB;
               for (JKL = 1; JKL <= LLB; JKL++) { // 170

                  // Transform from bandwidth JKL-1, JKU to JKL, JKU

                  for (JC = 1; JC <= min( N+JKL, M ) + JKU - 1; JC++) { // 160
                     EXTRA = ZERO;
                     ANGLE = TWOPI*SLARND( 1, ISEED );
                     C = COS( ANGLE );
                     S = SIN( ANGLE );
                     IROW = max( 1, JC-JKU );
                     if ( JC < N ) {
                        IL = min( M, JC+JKL ) + 1 - IROW;
                        slarot( false , JC > JKU, false , IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, EXTRA, DUMMY );
                     }

                     // Chase "EXTRA" back up

                     IC = JC;
                     IR = IROW;
                     for (JCH = JC - JKU; -JKL - JKU < 0 ? JCH >= 1 : JCH <= 1; JCH += -JKL - JKU) { // 150
                        if ( IC < N ) {
                           slartg(A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, C, S, DUMMY );
                        }
                        ICOL = max( 1, JCH-JKL );
                        IL = IC + 2 - ICOL;
                        TEMP = ZERO;
                        ILTEMP = JCH > JKL;
                        slarot( true , ILTEMP, true , IL, C, -S, A( IR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, TEMP, EXTRA );
                        if ( ILTEMP ) {
                           slartg(A( IR+1-ISKEW*( ICOL+1 )+IOFFST, ICOL+1 ), TEMP, C, S, DUMMY );
                           IROW = max( 1, JCH-JKL-JKU );
                           IL = IR + 2 - IROW;
                           EXTRA = ZERO;
                           slarot( false , JCH > JKL+JKU, true , IL, C, -S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, TEMP );
                           IC = ICOL;
                           IR = IROW;
                        }
                     } // 150
                  } // 160
               } // 170

            } else {

               // Bottom-Up -- Start at the bottom right.

               JKL = 0;
               for (JKU = 1; JKU <= UUB; JKU++) { // 200

                  // Transform from bandwidth JKL, JKU-1 to JKL, JKU

                  // First row actually rotated is M
                  // First column actually rotated is min( M+JKU, N )

                  IENDCH = min( M, N+JKL ) - 1;
                  for (JC = min( M+JKU, N ) - 1; JC >= 1 - JKL; JC--) { // 190
                     EXTRA = ZERO;
                     ANGLE = TWOPI*SLARND( 1, ISEED );
                     C = COS( ANGLE );
                     S = SIN( ANGLE );
                     IROW = max( 1, JC-JKU+1 );
                     if ( JC > 0 ) {
                        IL = min( M, JC+JKL+1 ) + 1 - IROW;
                        slarot( false , false , JC+JKL < M, IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, DUMMY, EXTRA );
                     }

                     // Chase "EXTRA" back down

                     IC = JC;
                     for (JCH = JC + JKL; JKL + JKU < 0 ? JCH >= IENDCH : JCH <= IENDCH; JCH += JKL + JKU) { // 180
                        ILEXTR = IC > 0;
                        if ( ILEXTR ) {
                           slartg(A( JCH-ISKEW*IC+IOFFST, IC ), EXTRA, C, S, DUMMY );
                        }
                        IC = max( 1, IC );
                        ICOL = min( N-1, JCH+JKU );
                        ILTEMP = JCH + JKU < N;
                        TEMP = ZERO;
                        slarot( true , ILEXTR, ILTEMP, ICOL+2-IC, C, S, A( JCH-ISKEW*IC+IOFFST, IC ), ILDA, EXTRA, TEMP );
                        if ( ILTEMP ) {
                           slartg(A( JCH-ISKEW*ICOL+IOFFST, ICOL ), TEMP, C, S, DUMMY );
                           IL = min( IENDCH, JCH+JKL+JKU ) + 2 - JCH;
                           EXTRA = ZERO;
                           slarot( false , true , JCH+JKL+JKU <= IENDCH, IL, C, S, A( JCH-ISKEW*ICOL+IOFFST, ICOL ), ILDA, TEMP, EXTRA );
                           IC = ICOL;
                        }
                     } // 180
                  } // 190
               } // 200

               JKU = UUB;
               for (JKL = 1; JKL <= LLB; JKL++) { // 230

                  // Transform from bandwidth JKL-1, JKU to JKL, JKU

                  // First row actually rotated is min( N+JKL, M )
                  // First column actually rotated is N

                  IENDCH = min( N, M+JKU ) - 1;
                  for (JR = min( N+JKL, M ) - 1; JR >= 1 - JKU; JR--) { // 220
                     EXTRA = ZERO;
                     ANGLE = TWOPI*SLARND( 1, ISEED );
                     C = COS( ANGLE );
                     S = SIN( ANGLE );
                     ICOL = max( 1, JR-JKL+1 );
                     if ( JR > 0 ) {
                        IL = min( N, JR+JKU+1 ) + 1 - ICOL;
                        slarot( true , false , JR+JKU < N, IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, DUMMY, EXTRA );
                     }

                     // Chase "EXTRA" back down

                     IR = JR;
                     for (JCH = JR + JKU; JKL + JKU < 0 ? JCH >= IENDCH : JCH <= IENDCH; JCH += JKL + JKU) { // 210
                        ILEXTR = IR > 0;
                        if ( ILEXTR ) {
                           slartg(A( IR-ISKEW*JCH+IOFFST, JCH ), EXTRA, C, S, DUMMY );
                        }
                        IR = max( 1, IR );
                        IROW = min( M-1, JCH+JKL );
                        ILTEMP = JCH + JKL < M;
                        TEMP = ZERO;
                        slarot( false , ILEXTR, ILTEMP, IROW+2-IR, C, S, A( IR-ISKEW*JCH+IOFFST, JCH ), ILDA, EXTRA, TEMP );
                        if ( ILTEMP ) {
                           slartg(A( IROW-ISKEW*JCH+IOFFST, JCH ), TEMP, C, S, DUMMY );
                           IL = min( IENDCH, JCH+JKL+JKU ) + 2 - JCH;
                           EXTRA = ZERO;
                           slarot( true , true , JCH+JKL+JKU <= IENDCH, IL, C, S, A( IROW-ISKEW*JCH+IOFFST, JCH ), ILDA, TEMP, EXTRA );
                           IR = IROW;
                        }
                     } // 210
                  } // 220
               } // 230
            }

         } else {

            // Symmetric -- A = U D U'

            IPACKG = IPACK;
            IOFFG = IOFFST;

            if ( TOPDWN ) {

               // Top-Down -- Generate Upper triangle only

               if ( IPACK >= 5 ) {
                  IPACKG = 6;
                  IOFFG = UUB + 1;
               } else {
                  IPACKG = 1;
               }
               scopy(MNMIN, D, 1, A( 1-ISKEW+IOFFG, 1 ), ILDA+1 );

               for (K = 1; K <= UUB; K++) { // 260
                  for (JC = 1; JC <= N - 1; JC++) { // 250
                     IROW = max( 1, JC-K );
                     IL = min( JC+1, K+2 );
                     EXTRA = ZERO;
                     TEMP = A( JC-ISKEW*( JC+1 )+IOFFG, JC+1 );
                     ANGLE = TWOPI*SLARND( 1, ISEED );
                     C = COS( ANGLE );
                     S = SIN( ANGLE );
                     slarot( false , JC > K, true , IL, C, S, A( IROW-ISKEW*JC+IOFFG, JC ), ILDA, EXTRA, TEMP );
                     slarot( true , true , false , min( K, N-JC )+1, C, S, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, TEMP, DUMMY );

                     // Chase EXTRA back up the matrix

                     ICOL = JC;
                     for (JCH = JC - K; -K < 0 ? JCH >= 1 : JCH <= 1; JCH += -K) { // 240
                        slartg(A( JCH+1-ISKEW*( ICOL+1 )+IOFFG, ICOL+1 ), EXTRA, C, S, DUMMY );
                        TEMP = A( JCH-ISKEW*( JCH+1 )+IOFFG, JCH+1 );
                        slarot( true , true , true , K+2, C, -S, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, TEMP, EXTRA );
                        IROW = max( 1, JCH-K );
                        IL = min( JCH+1, K+2 );
                        EXTRA = ZERO;
                        slarot( false , JCH > K, true , IL, C, -S, A( IROW-ISKEW*JCH+IOFFG, JCH ), ILDA, EXTRA, TEMP );
                        ICOL = JCH;
                     } // 240
                  } // 250
               } // 260

               // If we need lower triangle, copy from upper. Note that
               // the order of copying is chosen to work for 'q' -> 'b'

               if ( IPACK != IPACKG && IPACK != 3 ) {
                  for (JC = 1; JC <= N; JC++) { // 280
                     IROW = IOFFST - ISKEW*JC;
                     for (JR = JC; JR <= min( N, JC+UUB ); JR++) { // 270
                        A[JR+IROW, JC] = A( JC-ISKEW*JR+IOFFG, JR );
                     } // 270
                  } // 280
                  if ( IPACK == 5 ) {
                     for (JC = N - UUB + 1; JC <= N; JC++) { // 300
                        for (JR = N + 2 - JC; JR <= UUB + 1; JR++) { // 290
                           A[JR][JC] = ZERO;
                        } // 290
                     } // 300
                  }
                  if ( IPACKG == 6 ) {
                     IPACKG = IPACK;
                  } else {
                     IPACKG = 0;
                  }
               }
            } else {

               // Bottom-Up -- Generate Lower triangle only

               if ( IPACK >= 5 ) {
                  IPACKG = 5;
                  if (IPACK == 6) IOFFG = 1;
               } else {
                  IPACKG = 2;
               }
               scopy(MNMIN, D, 1, A( 1-ISKEW+IOFFG, 1 ), ILDA+1 );

               for (K = 1; K <= UUB; K++) { // 330
                  for (JC = N - 1; JC >= 1; JC--) { // 320
                     IL = min( N+1-JC, K+2 );
                     EXTRA = ZERO;
                     TEMP = A( 1+( 1-ISKEW )*JC+IOFFG, JC );
                     ANGLE = TWOPI*SLARND( 1, ISEED );
                     C = COS( ANGLE );
                     S = -SIN( ANGLE );
                     slarot( false , true , N-JC > K, IL, C, S, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, TEMP, EXTRA );
                     ICOL = max( 1, JC-K+1 );
                     slarot( true , false , true , JC+2-ICOL, C, S, A( JC-ISKEW*ICOL+IOFFG, ICOL ), ILDA, DUMMY, TEMP );

                     // Chase EXTRA back down the matrix

                     ICOL = JC;
                     for (JCH = JC + K; K < 0 ? JCH >= N - 1 : JCH <= N - 1; JCH += K) { // 310
                        slartg(A( JCH-ISKEW*ICOL+IOFFG, ICOL ), EXTRA, C, S, DUMMY );
                        TEMP = A( 1+( 1-ISKEW )*JCH+IOFFG, JCH );
                        slarot( true , true , true , K+2, C, S, A( JCH-ISKEW*ICOL+IOFFG, ICOL ), ILDA, EXTRA, TEMP );
                        IL = min( N+1-JCH, K+2 );
                        EXTRA = ZERO;
                        slarot( false , true , N-JCH > K, IL, C, S, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, TEMP, EXTRA );
                        ICOL = JCH;
                     } // 310
                  } // 320
               } // 330

               // If we need upper triangle, copy from lower. Note that
               // the order of copying is chosen to work for 'b' -> 'q'

               if ( IPACK != IPACKG && IPACK != 4 ) {
                  for (JC = N; JC >= 1; JC--) { // 350
                     IROW = IOFFST - ISKEW*JC;
                     for (JR = JC; JR >= max( 1, JC-UUB ); JR--) { // 340
                        A[JR+IROW, JC] = A( JC-ISKEW*JR+IOFFG, JR );
                     } // 340
                  } // 350
                  if ( IPACK == 6 ) {
                     for (JC = 1; JC <= UUB; JC++) { // 370
                        for (JR = 1; JR <= UUB + 1 - JC; JR++) { // 360
                           A[JR][JC] = ZERO;
                        } // 360
                     } // 370
                  }
                  if ( IPACKG == 5 ) {
                     IPACKG = IPACK;
                  } else {
                     IPACKG = 0;
                  }
               }
            }
         }

      } else {

         // 4)      Generate Banded Matrix by first
                 // Rotating by random Unitary matrices,
                 // then reducing the bandwidth using Householder
                 // transformations.

                 // Note: we should get here only if LDA >= N

         if ( ISYM == 1 ) {

            // Non-symmetric -- A = U D V

            slagge(MR, NC, LLB, UUB, D, A, LDA, ISEED, WORK, IINFO );
         } else {

            // Symmetric -- A = U D U'

            slagsy(M, LLB, D, A, LDA, ISEED, WORK, IINFO );

         }
         if ( IINFO != 0 ) {
            INFO = 3;
            return;
         }
      }

      // 5)      Pack the matrix

      if ( IPACK != IPACKG ) {
         if ( IPACK == 1 ) {

            // 'U' -- Upper triangular, not packed

            for (J = 1; J <= M; J++) { // 390
               for (I = J + 1; I <= M; I++) { // 380
                  A[I][J] = ZERO;
               } // 380
            } // 390

         } else if ( IPACK == 2 ) {

            // 'L' -- Lower triangular, not packed

            for (J = 2; J <= M; J++) { // 410
               for (I = 1; I <= J - 1; I++) { // 400
                  A[I][J] = ZERO;
               } // 400
            } // 410

         } else if ( IPACK == 3 ) {

            // 'C' -- Upper triangle packed Columnwise.

            ICOL = 1;
            IROW = 0;
            for (J = 1; J <= M; J++) { // 430
               for (I = 1; I <= J; I++) { // 420
                  IROW = IROW + 1;
                  if ( IROW > LDA ) {
                     IROW = 1;
                     ICOL = ICOL + 1;
                  }
                  A[IROW][ICOL] = A( I, J );
               } // 420
            } // 430

         } else if ( IPACK == 4 ) {

            // 'R' -- Lower triangle packed Columnwise.

            ICOL = 1;
            IROW = 0;
            for (J = 1; J <= M; J++) { // 450
               for (I = J; I <= M; I++) { // 440
                  IROW = IROW + 1;
                  if ( IROW > LDA ) {
                     IROW = 1;
                     ICOL = ICOL + 1;
                  }
                  A[IROW][ICOL] = A( I, J );
               } // 440
            } // 450

         } else if ( IPACK >= 5 ) {

            // 'B' -- The lower triangle is packed as a band matrix.
            // 'Q' -- The upper triangle is packed as a band matrix.
            // 'Z' -- The whole matrix is packed as a band matrix.

            if (IPACK == 5) UUB = 0;
            IF( IPACK == 6 ) LLB = 0;

            for (J = 1; J <= UUB; J++) { // 470
               for (I = min( J+LLB, M ); I >= 1; I--) { // 460
                  A[I-J+UUB+1, J] = A( I, J );
               } // 460
            } // 470

            for (J = UUB + 2; J <= N; J++) { // 490
               for (I = J - UUB; I <= min( J+LLB, M ); I++) { // 480
                  A[I-J+UUB+1, J] = A( I, J );
               } // 480
            } // 490
         }

         // If packed, zero out extraneous elements.

         // Symmetric/Triangular Packed --
         // zero out everything after A(IROW,ICOL)

         if ( IPACK == 3 || IPACK == 4 ) {
            for (JC = ICOL; JC <= M; JC++) { // 510
               for (JR = IROW + 1; JR <= LDA; JR++) { // 500
                  A[JR][JC] = ZERO;
               } // 500
               IROW = 0;
            } // 510

         } else if ( IPACK >= 5 ) {

            // Packed Band --
               // 1st row is now in A( UUB+2-j, j), zero above it
               // m-th row is now in A( M+UUB-j,j), zero below it
               // last non-zero diagonal is now in A( UUB+LLB+1,j ),
                  // zero below it, too.

            IR1 = UUB + LLB + 2;
            IR2 = UUB + M + 2;
            for (JC = 1; JC <= N; JC++) { // 540
               for (JR = 1; JR <= UUB + 1 - JC; JR++) { // 520
                  A[JR][JC] = ZERO;
               } // 520
               for (JR = max( 1, min( IR1; LDA < 0 ? JR >= IR2-JC ) ) : JR <= IR2-JC ) ); JR += LDA) { // 530
                  A[JR][JC] = ZERO;
               } // 530
            } // 540
         }
      }

      return;
      }
