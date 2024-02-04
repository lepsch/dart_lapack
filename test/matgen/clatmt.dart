      void clatmt(M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RANK, KL, KU, PACK, A, LDA, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double               COND, DMAX;
      int                INFO, KL, KU, LDA, M, MODE, N, RANK;
      String             DIST, PACK, SYM;
      // ..
      // .. Array Arguments ..
      Complex            A( LDA, * ), WORK( * );
      double               D( * );
      int                ISEED( 4 );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      double               ONE;
      const              ONE = 1.0 ;
      Complex            CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      double               TWOPI;
      const      TWOPI = 6.28318530717958647692528676655900576839 ;
      // ..
      // .. Local Scalars ..
      Complex            C, CT, CTEMP, DUMMY, EXTRA, S, ST;
      double               ALPHA, ANGLE, REALC, TEMP;
      int                I, IC, ICOL, IDIST, IENDCH, IINFO, IL, ILDA, IOFFG, IOFFST, IPACK, IPACKG, IR, IR1, IR2, IROW, IRSIGN, ISKEW, ISYM, ISYMPK, J, JC, JCH, JKL, JKU, JR, K, LLB, MINLDA, MNMIN, MR, NC, UUB;
      bool               CSYM, GIVENS, ILEXTR, ILTEMP, TOPDWN;
      // ..
      // .. External Functions ..
      //- COMPLEX            CLARND;
      //- REAL               SLARND;
      //- bool               lsame;
      // EXTERNAL CLARND, SLARND, lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLAGGE, CLAGHE, CLAGSY, CLAROT, CLARTG, CLASET, SLATM7, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, CONJG, COS, MAX, MIN, MOD, REAL, SIN
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
         CSYM = false;
      } else if ( lsame( SYM, 'P' ) ) {
         ISYM = 2;
         IRSIGN = 0;
         CSYM = false;
      } else if ( lsame( SYM, 'S' ) ) {
         ISYM = 2;
         IRSIGN = 0;
         CSYM = true;
      } else if ( lsame( SYM, 'H' ) ) {
         ISYM = 2;
         IRSIGN = 1;
         CSYM = false;
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
         xerbla('CLATMT', -INFO );
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

      if ( ( D( 1 ) ).abs() <= ( D( RANK ) ) ).abs() {
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

      claset('Full', LDA, N, CZERO, CZERO, A, LDA );

      // 3)      Generate Banded Matrix using Givens rotations.
              // Also the special case of UUB=LLB=0

                // Compute Addressing constants to cover all
                // storage formats.  Whether GE, HE, SY, GB, HB, or SB,
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

      // Diagonal Matrix -- We are done, unless it
      // is to be stored HP/SP/PP/TP (PACK='R' or 'C')

      if ( LLB == 0 && UUB == 0 ) {
         for (J = 1; J <= MNMIN; J++) { // 120
            A[( 1-ISKEW )*J+IOFFST, J] = CMPLX( D( J ) );
         } // 120

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

            for (J = 1; J <= MNMIN; J++) { // 130
               A[( 1-ISKEW )*J+IOFFST, J] = CMPLX( D( J ) );
            } // 130

            if ( TOPDWN ) {
               JKL = 0;
               for (JKU = 1; JKU <= UUB; JKU++) { // 160

                  // Transform from bandwidth JKL, JKU-1 to JKL, JKU

                  // Last row actually rotated is M
                  // Last column actually rotated is min( M+JKU, N )

                  for (JR = 1; JR <= min( M+JKU, N ) + JKL - 1; JR++) { // 150
                     EXTRA = CZERO;
                     ANGLE = TWOPI*SLARND( 1, ISEED );
                     C = COS( ANGLE )*CLARND( 5, ISEED );
                     S = SIN( ANGLE )*CLARND( 5, ISEED );
                     ICOL = max( 1, JR-JKL );
                     if ( JR < M ) {
                        IL = min( N, JR+JKU ) + 1 - ICOL;
                        clarot( true , JR > JKL, false , IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, EXTRA, DUMMY );
                     }

                     // Chase "EXTRA" back up

                     IR = JR;
                     IC = ICOL;
                     for (JCH = JR - JKL; -JKL - JKU < 0 ? JCH >= 1 : JCH <= 1; JCH += -JKL - JKU) { // 140
                        if ( IR < M ) {
                           clartg(A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED );
                           C = CONJG( REALC*DUMMY );
                           S = CONJG( -S*DUMMY );
                        }
                        IROW = max( 1, JCH-JKU );
                        IL = IR + 2 - IROW;
                        CTEMP = CZERO;
                        ILTEMP = JCH > JKU;
                        clarot( false , ILTEMP, true , IL, C, S, A( IROW-ISKEW*IC+IOFFST, IC ), ILDA, CTEMP, EXTRA );
                        if ( ILTEMP ) {
                           clartg(A( IROW+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), CTEMP, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED );
                           C = CONJG( REALC*DUMMY );
                           S = CONJG( -S*DUMMY );

                           ICOL = max( 1, JCH-JKU-JKL );
                           IL = IC + 2 - ICOL;
                           EXTRA = CZERO;
                           clarot( true , JCH > JKU+JKL, true , IL, C, S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, CTEMP );
                           IC = ICOL;
                           IR = IROW;
                        }
                     } // 140
                  } // 150
               } // 160

               JKU = UUB;
               for (JKL = 1; JKL <= LLB; JKL++) { // 190

                  // Transform from bandwidth JKL-1, JKU to JKL, JKU

                  for (JC = 1; JC <= min( N+JKL, M ) + JKU - 1; JC++) { // 180
                     EXTRA = CZERO;
                     ANGLE = TWOPI*SLARND( 1, ISEED );
                     C = COS( ANGLE )*CLARND( 5, ISEED );
                     S = SIN( ANGLE )*CLARND( 5, ISEED );
                     IROW = max( 1, JC-JKU );
                     if ( JC < N ) {
                        IL = min( M, JC+JKL ) + 1 - IROW;
                        clarot( false , JC > JKU, false , IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, EXTRA, DUMMY );
                     }

                     // Chase "EXTRA" back up

                     IC = JC;
                     IR = IROW;
                     for (JCH = JC - JKU; -JKL - JKU < 0 ? JCH >= 1 : JCH <= 1; JCH += -JKL - JKU) { // 170
                        if ( IC < N ) {
                           clartg(A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED );
                           C = CONJG( REALC*DUMMY );
                           S = CONJG( -S*DUMMY );
                        }
                        ICOL = max( 1, JCH-JKL );
                        IL = IC + 2 - ICOL;
                        CTEMP = CZERO;
                        ILTEMP = JCH > JKL;
                        clarot( true , ILTEMP, true , IL, C, S, A( IR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, CTEMP, EXTRA );
                        if ( ILTEMP ) {
                           clartg(A( IR+1-ISKEW*( ICOL+1 )+IOFFST, ICOL+1 ), CTEMP, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED );
                           C = CONJG( REALC*DUMMY );
                           S = CONJG( -S*DUMMY );
                           IROW = max( 1, JCH-JKL-JKU );
                           IL = IR + 2 - IROW;
                           EXTRA = CZERO;
                           clarot( false , JCH > JKL+JKU, true , IL, C, S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, CTEMP );
                           IC = ICOL;
                           IR = IROW;
                        }
                     } // 170
                  } // 180
               } // 190

            } else {

               // Bottom-Up -- Start at the bottom right.

               JKL = 0;
               for (JKU = 1; JKU <= UUB; JKU++) { // 220

                  // Transform from bandwidth JKL, JKU-1 to JKL, JKU

                  // First row actually rotated is M
                  // First column actually rotated is min( M+JKU, N )

                  IENDCH = min( M, N+JKL ) - 1;
                  for (JC = min( M+JKU, N ) - 1; JC >= 1 - JKL; JC--) { // 210
                     EXTRA = CZERO;
                     ANGLE = TWOPI*SLARND( 1, ISEED );
                     C = COS( ANGLE )*CLARND( 5, ISEED );
                     S = SIN( ANGLE )*CLARND( 5, ISEED );
                     IROW = max( 1, JC-JKU+1 );
                     if ( JC > 0 ) {
                        IL = min( M, JC+JKL+1 ) + 1 - IROW;
                        clarot( false , false , JC+JKL < M, IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, DUMMY, EXTRA );
                     }

                     // Chase "EXTRA" back down

                     IC = JC;
                     for (JCH = JC + JKL; JKL + JKU < 0 ? JCH >= IENDCH : JCH <= IENDCH; JCH += JKL + JKU) { // 200
                        ILEXTR = IC > 0;
                        if ( ILEXTR ) {
                           clartg(A( JCH-ISKEW*IC+IOFFST, IC ), EXTRA, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED );
                           C = REALC*DUMMY;
                           S = S*DUMMY;
                        }
                        IC = max( 1, IC );
                        ICOL = min( N-1, JCH+JKU );
                        ILTEMP = JCH + JKU < N;
                        CTEMP = CZERO;
                        clarot( true , ILEXTR, ILTEMP, ICOL+2-IC, C, S, A( JCH-ISKEW*IC+IOFFST, IC ), ILDA, EXTRA, CTEMP );
                        if ( ILTEMP ) {
                           clartg(A( JCH-ISKEW*ICOL+IOFFST, ICOL ), CTEMP, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED );
                           C = REALC*DUMMY;
                           S = S*DUMMY;
                           IL = min( IENDCH, JCH+JKL+JKU ) + 2 - JCH;
                           EXTRA = CZERO;
                           clarot( false , true , JCH+JKL+JKU <= IENDCH, IL, C, S, A( JCH-ISKEW*ICOL+IOFFST, ICOL ), ILDA, CTEMP, EXTRA );
                           IC = ICOL;
                        }
                     } // 200
                  } // 210
               } // 220

               JKU = UUB;
               for (JKL = 1; JKL <= LLB; JKL++) { // 250

                  // Transform from bandwidth JKL-1, JKU to JKL, JKU

                  // First row actually rotated is min( N+JKL, M )
                  // First column actually rotated is N

                  IENDCH = min( N, M+JKU ) - 1;
                  for (JR = min( N+JKL, M ) - 1; JR >= 1 - JKU; JR--) { // 240
                     EXTRA = CZERO;
                     ANGLE = TWOPI*SLARND( 1, ISEED );
                     C = COS( ANGLE )*CLARND( 5, ISEED );
                     S = SIN( ANGLE )*CLARND( 5, ISEED );
                     ICOL = max( 1, JR-JKL+1 );
                     if ( JR > 0 ) {
                        IL = min( N, JR+JKU+1 ) + 1 - ICOL;
                        clarot( true , false , JR+JKU < N, IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, DUMMY, EXTRA );
                     }

                     // Chase "EXTRA" back down

                     IR = JR;
                     for (JCH = JR + JKU; JKL + JKU < 0 ? JCH >= IENDCH : JCH <= IENDCH; JCH += JKL + JKU) { // 230
                        ILEXTR = IR > 0;
                        if ( ILEXTR ) {
                           clartg(A( IR-ISKEW*JCH+IOFFST, JCH ), EXTRA, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED );
                           C = REALC*DUMMY;
                           S = S*DUMMY;
                        }
                        IR = max( 1, IR );
                        IROW = min( M-1, JCH+JKL );
                        ILTEMP = JCH + JKL < M;
                        CTEMP = CZERO;
                        clarot( false , ILEXTR, ILTEMP, IROW+2-IR, C, S, A( IR-ISKEW*JCH+IOFFST, JCH ), ILDA, EXTRA, CTEMP );
                        if ( ILTEMP ) {
                           clartg(A( IROW-ISKEW*JCH+IOFFST, JCH ), CTEMP, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED );
                           C = REALC*DUMMY;
                           S = S*DUMMY;
                           IL = min( IENDCH, JCH+JKL+JKU ) + 2 - JCH;
                           EXTRA = CZERO;
                           clarot( true , true , JCH+JKL+JKU <= IENDCH, IL, C, S, A( IROW-ISKEW*JCH+IOFFST, JCH ), ILDA, CTEMP, EXTRA );
                           IR = IROW;
                        }
                     } // 230
                  } // 240
               } // 250

            }

         } else {

            // Symmetric -- A = U D U'
            // Hermitian -- A = U D U*

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

               for (J = 1; J <= MNMIN; J++) { // 260
                  A[( 1-ISKEW )*J+IOFFG, J] = CMPLX( D( J ) );
               } // 260

               for (K = 1; K <= UUB; K++) { // 290
                  for (JC = 1; JC <= N - 1; JC++) { // 280
                     IROW = max( 1, JC-K );
                     IL = min( JC+1, K+2 );
                     EXTRA = CZERO;
                     CTEMP = A( JC-ISKEW*( JC+1 )+IOFFG, JC+1 );
                     ANGLE = TWOPI*SLARND( 1, ISEED );
                     C = COS( ANGLE )*CLARND( 5, ISEED );
                     S = SIN( ANGLE )*CLARND( 5, ISEED );
                     if ( CSYM ) {
                        CT = C;
                        ST = S;
                     } else {
                        CTEMP = CONJG( CTEMP );
                        CT = CONJG( C );
                        ST = CONJG( S );
                     }
                     clarot( false , JC > K, true , IL, C, S, A( IROW-ISKEW*JC+IOFFG, JC ), ILDA, EXTRA, CTEMP );
                     clarot( true , true , false , min( K, N-JC )+1, CT, ST, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, CTEMP, DUMMY );

                     // Chase EXTRA back up the matrix

                     ICOL = JC;
                     for (JCH = JC - K; -K < 0 ? JCH >= 1 : JCH <= 1; JCH += -K) { // 270
                        clartg(A( JCH+1-ISKEW*( ICOL+1 )+IOFFG, ICOL+1 ), EXTRA, REALC, S, DUMMY );
                        DUMMY = CLARND( 5, ISEED );
                        C = CONJG( REALC*DUMMY );
                        S = CONJG( -S*DUMMY );
                        CTEMP = A( JCH-ISKEW*( JCH+1 )+IOFFG, JCH+1 );
                        if ( CSYM ) {
                           CT = C;
                           ST = S;
                        } else {
                           CTEMP = CONJG( CTEMP );
                           CT = CONJG( C );
                           ST = CONJG( S );
                        }
                        clarot( true , true , true , K+2, C, S, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, CTEMP, EXTRA );
                        IROW = max( 1, JCH-K );
                        IL = min( JCH+1, K+2 );
                        EXTRA = CZERO;
                        clarot( false , JCH > K, true , IL, CT, ST, A( IROW-ISKEW*JCH+IOFFG, JCH ), ILDA, EXTRA, CTEMP );
                        ICOL = JCH;
                     } // 270
                  } // 280
               } // 290

               // If we need lower triangle, copy from upper. Note that
               // the order of copying is chosen to work for 'q' -> 'b'

               if ( IPACK != IPACKG && IPACK != 3 ) {
                  for (JC = 1; JC <= N; JC++) { // 320
                     IROW = IOFFST - ISKEW*JC;
                     if ( CSYM ) {
                        for (JR = JC; JR <= min( N, JC+UUB ); JR++) { // 300
                           A[JR+IROW, JC] = A( JC-ISKEW*JR+IOFFG, JR );
                        } // 300
                     } else {
                        for (JR = JC; JR <= min( N, JC+UUB ); JR++) { // 310
                           A[JR+IROW, JC] = CONJG( A( JC-ISKEW*JR+ IOFFG, JR ) );
                        } // 310
                     }
                  } // 320
                  if ( IPACK == 5 ) {
                     for (JC = N - UUB + 1; JC <= N; JC++) { // 340
                        for (JR = N + 2 - JC; JR <= UUB + 1; JR++) { // 330
                           A[JR, JC] = CZERO;
                        } // 330
                     } // 340
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

               for (J = 1; J <= MNMIN; J++) { // 350
                  A[( 1-ISKEW )*J+IOFFG, J] = CMPLX( D( J ) );
               } // 350

               for (K = 1; K <= UUB; K++) { // 380
                  for (JC = N - 1; JC >= 1; JC--) { // 370
                     IL = min( N+1-JC, K+2 );
                     EXTRA = CZERO;
                     CTEMP = A( 1+( 1-ISKEW )*JC+IOFFG, JC );
                     ANGLE = TWOPI*SLARND( 1, ISEED );
                     C = COS( ANGLE )*CLARND( 5, ISEED );
                     S = SIN( ANGLE )*CLARND( 5, ISEED );
                     if ( CSYM ) {
                        CT = C;
                        ST = S;
                     } else {
                        CTEMP = CONJG( CTEMP );
                        CT = CONJG( C );
                        ST = CONJG( S );
                     }
                     clarot( false , true , N-JC > K, IL, C, S, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, CTEMP, EXTRA );
                     ICOL = max( 1, JC-K+1 );
                     clarot( true , false , true , JC+2-ICOL, CT, ST, A( JC-ISKEW*ICOL+IOFFG, ICOL ), ILDA, DUMMY, CTEMP );

                     // Chase EXTRA back down the matrix

                     ICOL = JC;
                     for (JCH = JC + K; K < 0 ? JCH >= N - 1 : JCH <= N - 1; JCH += K) { // 360
                        clartg(A( JCH-ISKEW*ICOL+IOFFG, ICOL ), EXTRA, REALC, S, DUMMY );
                        DUMMY = CLARND( 5, ISEED );
                        C = REALC*DUMMY;
                        S = S*DUMMY;
                        CTEMP = A( 1+( 1-ISKEW )*JCH+IOFFG, JCH );
                        if ( CSYM ) {
                           CT = C;
                           ST = S;
                        } else {
                           CTEMP = CONJG( CTEMP );
                           CT = CONJG( C );
                           ST = CONJG( S );
                        }
                        clarot( true , true , true , K+2, C, S, A( JCH-ISKEW*ICOL+IOFFG, ICOL ), ILDA, EXTRA, CTEMP );
                        IL = min( N+1-JCH, K+2 );
                        EXTRA = CZERO;
                        clarot( false , true , N-JCH > K, IL, CT, ST, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, CTEMP, EXTRA );
                        ICOL = JCH;
                     } // 360
                  } // 370
               } // 380

               // If we need upper triangle, copy from lower. Note that
               // the order of copying is chosen to work for 'b' -> 'q'

               if ( IPACK != IPACKG && IPACK != 4 ) {
                  for (JC = N; JC >= 1; JC--) { // 410
                     IROW = IOFFST - ISKEW*JC;
                     if ( CSYM ) {
                        for (JR = JC; JR >= max( 1, JC-UUB ); JR--) { // 390
                           A[JR+IROW, JC] = A( JC-ISKEW*JR+IOFFG, JR );
                        } // 390
                     } else {
                        for (JR = JC; JR >= max( 1, JC-UUB ); JR--) { // 400
                           A[JR+IROW, JC] = CONJG( A( JC-ISKEW*JR+ IOFFG, JR ) );
                        } // 400
                     }
                  } // 410
                  if ( IPACK == 6 ) {
                     for (JC = 1; JC <= UUB; JC++) { // 430
                        for (JR = 1; JR <= UUB + 1 - JC; JR++) { // 420
                           A[JR, JC] = CZERO;
                        } // 420
                     } // 430
                  }
                  if ( IPACKG == 5 ) {
                     IPACKG = IPACK;
                  } else {
                     IPACKG = 0;
                  }
               }
            }

            // Ensure that the diagonal is real if Hermitian

            if ( !CSYM ) {
               for (JC = 1; JC <= N; JC++) { // 440
                  IROW = IOFFST + ( 1-ISKEW )*JC;
                  A[IROW, JC] = CMPLX( double( A( IROW, JC ) ) );
               } // 440
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

            clagge(MR, NC, LLB, UUB, D, A, LDA, ISEED, WORK, IINFO );
         } else {

            // Symmetric -- A = U D U' or
            // Hermitian -- A = U D U*

            if ( CSYM ) {
               clagsy(M, LLB, D, A, LDA, ISEED, WORK, IINFO );
            } else {
               claghe(M, LLB, D, A, LDA, ISEED, WORK, IINFO );
            }
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

            for (J = 1; J <= M; J++) { // 460
               for (I = J + 1; I <= M; I++) { // 450
                  A[I, J] = CZERO;
               } // 450
            } // 460

         } else if ( IPACK == 2 ) {

            // 'L' -- Lower triangular, not packed

            for (J = 2; J <= M; J++) { // 480
               for (I = 1; I <= J - 1; I++) { // 470
                  A[I, J] = CZERO;
               } // 470
            } // 480

         } else if ( IPACK == 3 ) {

            // 'C' -- Upper triangle packed Columnwise.

            ICOL = 1;
            IROW = 0;
            for (J = 1; J <= M; J++) { // 500
               for (I = 1; I <= J; I++) { // 490
                  IROW = IROW + 1;
                  if ( IROW > LDA ) {
                     IROW = 1;
                     ICOL = ICOL + 1;
                  }
                  A[IROW, ICOL] = A( I, J );
               } // 490
            } // 500

         } else if ( IPACK == 4 ) {

            // 'R' -- Lower triangle packed Columnwise.

            ICOL = 1;
            IROW = 0;
            for (J = 1; J <= M; J++) { // 520
               for (I = J; I <= M; I++) { // 510
                  IROW = IROW + 1;
                  if ( IROW > LDA ) {
                     IROW = 1;
                     ICOL = ICOL + 1;
                  }
                  A[IROW, ICOL] = A( I, J );
               } // 510
            } // 520

         } else if ( IPACK >= 5 ) {

            // 'B' -- The lower triangle is packed as a band matrix.
            // 'Q' -- The upper triangle is packed as a band matrix.
            // 'Z' -- The whole matrix is packed as a band matrix.

            if (IPACK == 5) UUB = 0;
            IF( IPACK == 6 ) LLB = 0;

            for (J = 1; J <= UUB; J++) { // 540
               for (I = min( J+LLB, M ); I >= 1; I--) { // 530
                  A[I-J+UUB+1, J] = A( I, J );
               } // 530
            } // 540

            for (J = UUB + 2; J <= N; J++) { // 560
               for (I = J - UUB; I <= min( J+LLB, M ); I++) { // 550
                  A[I-J+UUB+1, J] = A( I, J );
               } // 550
            } // 560
         }

         // If packed, zero out extraneous elements.

         // Symmetric/Triangular Packed --
         // zero out everything after A(IROW,ICOL)

         if ( IPACK == 3 || IPACK == 4 ) {
            for (JC = ICOL; JC <= M; JC++) { // 580
               for (JR = IROW + 1; JR <= LDA; JR++) { // 570
                  A[JR, JC] = CZERO;
               } // 570
               IROW = 0;
            } // 580

         } else if ( IPACK >= 5 ) {

            // Packed Band --
               // 1st row is now in A( UUB+2-j, j), zero above it
               // m-th row is now in A( M+UUB-j,j), zero below it
               // last non-zero diagonal is now in A( UUB+LLB+1,j ),
                  // zero below it, too.

            IR1 = UUB + LLB + 2;
            IR2 = UUB + M + 2;
            for (JC = 1; JC <= N; JC++) { // 610
               for (JR = 1; JR <= UUB + 1 - JC; JR++) { // 590
                  A[JR, JC] = CZERO;
               } // 590
               for (JR = max( 1, min( IR1; LDA < 0 ? JR >= IR2-JC ) ) : JR <= IR2-JC ) ); JR += LDA) { // 600
                  A[JR, JC] = CZERO;
               } // 600
            } // 610
         }
      }

      return;
      }