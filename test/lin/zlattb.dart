      void zlattb(final int IMAT, final int UPLO, final int TRANS, final int DIAG, final Array<int> ISEED_, final int N, final int KD, final Matrix<double> AB_, final int LDAB, final int B, final Array<double> WORK_, final Array<double> RWORK_, final Box<int> INFO,) {
  final ISEED = ISEED_.having();
  final AB = AB_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                IMAT, INFO, KD, LDAB, N;
      int                ISEED( 4 );
      double             RWORK( * );
      Complex         AB( LDAB, * ), B( * ), WORK( * );
      // ..

      double             ONE, TWO, ZERO;
      const              ONE = 1.0, TWO = 2.0, ZERO = 0.0 ;
      bool               UPPER;
      String             DIST, PACKIT, TYPE;
      String             PATH;
      int                I, IOFF, IY, J, JCOUNT, KL, KU, LENJ, MODE;
      double             ANORM, BIGNUM, BNORM, BSCAL, CNDNUM, REXP, SFAC, SMLNUM, TEXP, TLEFT, TNORM, TSCAL, ULP, UNFL;
      Complex         PLUS1, PLUS2, STAR1;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                IZAMAX;
      //- double             DLAMCH, DLARND;
      //- Complex         ZLARND;
      // EXTERNAL lsame, IZAMAX, DLAMCH, DLARND, ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARNV, ZCOPY, ZDSCAL, ZLARNV, ZLATB4, ZLATMS, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, MAX, MIN, SQRT

      PATH[1: 1] = 'Zomplex precision';
      PATH[2: 3] = 'TB';
      UNFL = dlamch( 'Safe minimum' );
      ULP = dlamch( 'Epsilon' )*dlamch( 'Base' );
      SMLNUM = UNFL;
      BIGNUM = ( ONE-ULP ) / SMLNUM;
      if ( ( IMAT >= 6 && IMAT <= 9 ) || IMAT == 17 ) {
         DIAG = 'U';
      } else {
         DIAG = 'N';
      }
      INFO = 0;

      // Quick return if N <= 0.

      if (N <= 0) return;

      // Call ZLATB4 to set parameters for ZLATMS.

      UPPER = lsame( UPLO, 'U' );
      if ( UPPER ) {
         zlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
         KU = KD;
         IOFF = 1 + max( 0, KD-N+1 );
         KL = 0;
         PACKIT = 'Q';
      } else {
         zlatb4(PATH, -IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
         KL = KD;
         IOFF = 1;
         KU = 0;
         PACKIT = 'B';
      }

      // IMAT <= 5:  Non-unit triangular matrix

      if ( IMAT <= 5 ) {
         zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, AB( IOFF, 1 ), LDAB, WORK, INFO );

      // IMAT > 5:  Unit triangular matrix
      // The diagonal is deliberately set to something other than 1.

      // IMAT = 6:  Matrix is the identity

      } else if ( IMAT == 6 ) {
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 20
               for (I = max( 1, KD+2-J ); I <= KD; I++) { // 10
                  AB[I][J] = ZERO;
               } // 10
               AB[KD+1][J] = J;
            } // 20
         } else {
            for (J = 1; J <= N; J++) { // 40
               AB[1][J] = J;
               for (I = 2; I <= min( KD+1, N-J+1 ); I++) { // 30
                  AB[I][J] = ZERO;
               } // 30
            } // 40
         }

      // IMAT > 6:  Non-trivial unit triangular matrix

      // A unit triangular matrix T with condition CNDNUM is formed.
      // In this version, T only has bandwidth 2, the rest of it is zero.

      } else if ( IMAT <= 9 ) {
         TNORM = sqrt( CNDNUM );

         // Initialize AB to zero.

         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 60
               for (I = max( 1, KD+2-J ); I <= KD; I++) { // 50
                  AB[I][J] = ZERO;
               } // 50
               AB[KD+1][J] = J.toDouble();
            } // 60
         } else {
            for (J = 1; J <= N; J++) { // 80
               for (I = 2; I <= min( KD+1, N-J+1 ); I++) { // 70
                  AB[I][J] = ZERO;
               } // 70
               AB[1][J] = J.toDouble();
            } // 80
         }

         // Special case:  T is tridiagonal.  Set every other offdiagonal
         // so that the matrix has norm TNORM+1.

         if ( KD == 1 ) {
            if ( UPPER ) {
               AB[1][2] = TNORM*ZLARND( 5, ISEED );
               LENJ = ( N-3 ) / 2;
               zlarnv(2, ISEED, LENJ, WORK );
               for (J = 1; J <= LENJ; J++) { // 90
                  AB[1][2*( J+1 )] = TNORM*WORK( J );
               } // 90
            } else {
               AB[2][1] = TNORM*ZLARND( 5, ISEED );
               LENJ = ( N-3 ) / 2;
               zlarnv(2, ISEED, LENJ, WORK );
               for (J = 1; J <= LENJ; J++) { // 100
                  AB[2][2*J+1] = TNORM*WORK( J );
               } // 100
            }
         } else if ( KD > 1 ) {

            // Form a unit triangular matrix T with condition CNDNUM.  T is
            // given by
            //         | 1   +   *                      |
            //         |     1   +                      |
            //     T = |         1   +   *              |
            //         |             1   +              |
            //         |                 1   +   *      |
            //         |                     1   +      |
            //         |                          . . . |
         // Each element marked with a '*' is formed by taking the product
         // of the adjacent elements marked with '+'.  The '*'s can be
         // chosen freely, and the '+'s are chosen so that the inverse of
         // T will have elements of the same magnitude as T.

         // The two offdiagonals of T are stored in WORK.

            STAR1 = TNORM*ZLARND( 5, ISEED );
            SFAC = sqrt( TNORM );
            PLUS1 = SFAC*ZLARND( 5, ISEED );
            for (J = 1; J <= N; J += 2) { // 110
               PLUS2 = STAR1 / PLUS1;
               WORK[J] = PLUS1;
               WORK[N+J] = STAR1;
               if ( J+1 <= N ) {
                  WORK[J+1] = PLUS2;
                  WORK[N+J+1] = ZERO;
                  PLUS1 = STAR1 / PLUS2;

                  // Generate a new *-value with norm between sqrt(TNORM)
                  // and TNORM.

                  REXP = dlarnd( 2, ISEED );
                  if ( REXP < ZERO ) {
                     STAR1 = -SFAC**( ONE-REXP )*ZLARND( 5, ISEED );
                  } else {
                     STAR1 = SFAC**( ONE+REXP )*ZLARND( 5, ISEED );
                  }
               }
            } // 110

            // Copy the tridiagonal T to AB.

            if ( UPPER ) {
               zcopy(N-1, WORK, 1, AB( KD, 2 ), LDAB );
               zcopy(N-2, WORK( N+1 ), 1, AB( KD-1, 3 ), LDAB );
            } else {
               zcopy(N-1, WORK, 1, AB( 2, 1 ), LDAB );
               zcopy(N-2, WORK( N+1 ), 1, AB( 3, 1 ), LDAB );
            }
         }

      // IMAT > 9:  Pathological test cases.  These triangular matrices
      // are badly scaled or badly conditioned, so when used in solving a
      // triangular system they may cause overflow in the solution vector.

      } else if ( IMAT == 10 ) {

         // Type 10:  Generate a triangular matrix with elements between
         // -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
         // Make the right hand side large so that it requires scaling.

         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 120
               LENJ = min( J-1, KD );
               zlarnv(4, ISEED, LENJ, AB( KD+1-LENJ, J ) );
               AB[KD+1][J] = ZLARND( 5, ISEED )*TWO;
            } // 120
         } else {
            for (J = 1; J <= N; J++) { // 130
               LENJ = min( N-J, KD );
               if (LENJ > 0) zlarnv( 4, ISEED, LENJ, AB( 2, J ) );
               AB[1][J] = ZLARND( 5, ISEED )*TWO;
            } // 130
         }

         // Set the right hand side so that the largest value is BIGNUM.

         zlarnv(2, ISEED, N, B );
         IY = IZAMAX( N, B, 1 );
         BNORM = ( B( IY ) ).abs();
         BSCAL = BIGNUM / max( ONE, BNORM );
         zdscal(N, BSCAL, B, 1 );

      } else if ( IMAT == 11 ) {

         // Type 11:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 11, the offdiagonal elements are small (CNORM(j) < 1).

         zlarnv(2, ISEED, N, B );
         TSCAL = ONE / (KD+1).toDouble();
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 140
               LENJ = min( J-1, KD );
               if ( LENJ > 0 ) {
                  zlarnv(4, ISEED, LENJ, AB( KD+2-LENJ, J ) );
                  zdscal(LENJ, TSCAL, AB( KD+2-LENJ, J ), 1 );
               }
               AB[KD+1][J] = ZLARND( 5, ISEED );
            } // 140
            AB[KD+1][N] = SMLNUM*AB( KD+1, N );
         } else {
            for (J = 1; J <= N; J++) { // 150
               LENJ = min( N-J, KD );
               if ( LENJ > 0 ) {
                  zlarnv(4, ISEED, LENJ, AB( 2, J ) );
                  zdscal(LENJ, TSCAL, AB( 2, J ), 1 );
               }
               AB[1][J] = ZLARND( 5, ISEED );
            } // 150
            AB[1][1] = SMLNUM*AB( 1, 1 );
         }

      } else if ( IMAT == 12 ) {

         // Type 12:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 12, the offdiagonal elements are O(1) (CNORM(j) > 1).

         zlarnv(2, ISEED, N, B );
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 160
               LENJ = min( J-1, KD );
               if (LENJ > 0) zlarnv( 4, ISEED, LENJ, AB( KD+2-LENJ, J ) );
               AB[KD+1][J] = ZLARND( 5, ISEED );
            } // 160
            AB[KD+1][N] = SMLNUM*AB( KD+1, N );
         } else {
            for (J = 1; J <= N; J++) { // 170
               LENJ = min( N-J, KD );
               if (LENJ > 0) zlarnv( 4, ISEED, LENJ, AB( 2, J ) );
               AB[1][J] = ZLARND( 5, ISEED );
            } // 170
            AB[1][1] = SMLNUM*AB( 1, 1 );
         }

      } else if ( IMAT == 13 ) {

         // Type 13:  T is diagonal with small numbers on the diagonal to
         // make the growth factor underflow, but a small right hand side
         // chosen so that the solution does not overflow.

         if ( UPPER ) {
            JCOUNT = 1;
            for (J = N; J >= 1; J--) { // 190
               DO 180 I = max( 1, KD+1-( J-1 ) ), KD;
                  AB[I][J] = ZERO;
               } // 180
               if ( JCOUNT <= 2 ) {
                  AB[KD+1][J] = SMLNUM*ZLARND( 5, ISEED );
               } else {
                  AB[KD+1][J] = ZLARND( 5, ISEED );
               }
               JCOUNT = JCOUNT + 1;
               if (JCOUNT > 4) JCOUNT = 1;
            } // 190
         } else {
            JCOUNT = 1;
            for (J = 1; J <= N; J++) { // 210
               for (I = 2; I <= min( N-J+1, KD+1 ); I++) { // 200
                  AB[I][J] = ZERO;
               } // 200
               if ( JCOUNT <= 2 ) {
                  AB[1][J] = SMLNUM*ZLARND( 5, ISEED );
               } else {
                  AB[1][J] = ZLARND( 5, ISEED );
               }
               JCOUNT = JCOUNT + 1;
               if (JCOUNT > 4) JCOUNT = 1;
            } // 210
         }

         // Set the right hand side alternately zero and small.

         if ( UPPER ) {
            B[1] = ZERO;
            for (I = N; I >= 2; I -= 2) { // 220
               B[I] = ZERO;
               B[I-1] = SMLNUM*ZLARND( 5, ISEED );
            } // 220
         } else {
            B[N] = ZERO;
            for (I = 1; 2 < 0 ? I >= N - 1 : I <= N - 1; I += 2) { // 230
               B[I] = ZERO;
               B[I+1] = SMLNUM*ZLARND( 5, ISEED );
            } // 230
         }

      } else if ( IMAT == 14 ) {

         // Type 14:  Make the diagonal elements small to cause gradual
         // overflow when dividing by T(j,j).  To control the amount of
         // scaling needed, the matrix is bidiagonal.

         TEXP = ONE / (KD+1).toDouble();
         TSCAL = SMLNUM**TEXP;
         zlarnv(4, ISEED, N, B );
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 250
               for (I = max( 1, KD+2-J ); I <= KD; I++) { // 240
                  AB[I][J] = ZERO;
               } // 240
               if (J > 1 && KD > 0) AB( KD, J ) = DCMPLX( -ONE, -ONE );
               AB[KD+1][J] = TSCAL*ZLARND( 5, ISEED );
            } // 250
            B[N] = DCMPLX( ONE, ONE );
         } else {
            for (J = 1; J <= N; J++) { // 270
               for (I = 3; I <= min( N-J+1, KD+1 ); I++) { // 260
                  AB[I][J] = ZERO;
               } // 260
               if (J < N && KD > 0) AB( 2, J ) = DCMPLX( -ONE, -ONE );
               AB[1][J] = TSCAL*ZLARND( 5, ISEED );
            } // 270
            B[1] = DCMPLX( ONE, ONE );
         }

      } else if ( IMAT == 15 ) {

         // Type 15:  One zero diagonal element.

         IY = N / 2 + 1;
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 280
               LENJ = min( J, KD+1 );
               zlarnv(4, ISEED, LENJ, AB( KD+2-LENJ, J ) );
               if ( J != IY ) {
                  AB[KD+1][J] = ZLARND( 5, ISEED )*TWO;
               } else {
                  AB[KD+1][J] = ZERO;
               }
            } // 280
         } else {
            for (J = 1; J <= N; J++) { // 290
               LENJ = min( N-J+1, KD+1 );
               zlarnv(4, ISEED, LENJ, AB( 1, J ) );
               if ( J != IY ) {
                  AB[1][J] = ZLARND( 5, ISEED )*TWO;
               } else {
                  AB[1][J] = ZERO;
               }
            } // 290
         }
         zlarnv(2, ISEED, N, B );
         zdscal(N, TWO, B, 1 );

      } else if ( IMAT == 16 ) {

         // Type 16:  Make the offdiagonal elements large to cause overflow
         // when adding a column of T.  In the non-transposed case, the
         // matrix is constructed to cause overflow when adding a column in
         // every other step.

         TSCAL = UNFL / ULP;
         TSCAL = ( ONE-ULP ) / TSCAL;
         for (J = 1; J <= N; J++) { // 310
            for (I = 1; I <= KD + 1; I++) { // 300
               AB[I][J] = ZERO;
            } // 300
         } // 310
         TEXP = ONE;
         if ( KD > 0 ) {
            if ( UPPER ) {
               for (J = N; -KD < 0 ? J >= 1 : J <= 1; J += -KD) { // 330
                  for (I = J; I >= max( 1, J-KD+1 ); I -= 2) { // 320
                     AB[1+( J-I )][I] = -TSCAL / (KD+2).toDouble();
                     AB[KD+1][I] = ONE;
                     B[I] = TEXP*( ONE-ULP );
                     if ( I > max( 1, J-KD+1 ) ) {
                        AB[2+( J-I )][I-1] = -( TSCAL / DBLE( KD+2 ) ) / (KD+3).toDouble();
                        AB[KD+1][I-1] = ONE;
                        B[I-1] = TEXP*DBLE( ( KD+1 )*( KD+1 )+KD );
                     }
                     TEXP = TEXP*TWO;
                  } // 320
                  B[max( 1][J-KD+1 )] = ( DBLE( KD+2 ) / (KD+3).toDouble() )*TSCAL;
               } // 330
            } else {
               for (J = 1; KD < 0 ? J >= N : J <= N; J += KD) { // 350
                  TEXP = ONE;
                  LENJ = min( KD+1, N-J+1 );
                  for (I = J; 2 < 0 ? I >= min( N, J+KD-1 ) : I <= min( N, J+KD-1 ); I += 2) { // 340
                     AB[LENJ-( I-J )][J] = -TSCAL / (KD+2).toDouble();
                     AB[1][J] = ONE;
                     B[J] = TEXP*( ONE-ULP );
                     if ( I < min( N, J+KD-1 ) ) {
                        AB[LENJ-( I-J+1 )][I+1] = -( TSCAL / DBLE( KD+2 ) ) / (KD+3).toDouble();
                        AB[1][I+1] = ONE;
                        B[I+1] = TEXP*DBLE( ( KD+1 )*( KD+1 )+KD );
                     }
                     TEXP = TEXP*TWO;
                  } // 340
                  B[min( N][J+KD-1 )] = ( DBLE( KD+2 ) / (KD+3).toDouble() )*TSCAL;
               } // 350
            }
         }

      } else if ( IMAT == 17 ) {

         // Type 17:  Generate a unit triangular matrix with elements
         // between -1 and 1, and make the right hand side large so that it
         // requires scaling.

         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 360
               LENJ = min( J-1, KD );
               zlarnv(4, ISEED, LENJ, AB( KD+1-LENJ, J ) );
               AB[KD+1][J] = J.toDouble();
            } // 360
         } else {
            for (J = 1; J <= N; J++) { // 370
               LENJ = min( N-J, KD );
               if (LENJ > 0) zlarnv( 4, ISEED, LENJ, AB( 2, J ) );
               AB[1][J] = J.toDouble();
            } // 370
         }

         // Set the right hand side so that the largest value is BIGNUM.

         zlarnv(2, ISEED, N, B );
         IY = IZAMAX( N, B, 1 );
         BNORM = ( B( IY ) ).abs();
         BSCAL = BIGNUM / max( ONE, BNORM );
         zdscal(N, BSCAL, B, 1 );

      } else if ( IMAT == 18 ) {

         // Type 18:  Generate a triangular matrix with elements between
         // BIGNUM/(KD+1) and BIGNUM so that at least one of the column
         // norms will exceed BIGNUM.
         // 1/3/91:  ZLATBS no longer can handle this case

         TLEFT = BIGNUM / (KD+1).toDouble();
         TSCAL = BIGNUM*( DBLE( KD+1 ) / (KD+2).toDouble() );
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 390
               LENJ = min( J, KD+1 );
               zlarnv(5, ISEED, LENJ, AB( KD+2-LENJ, J ) );
               dlarnv(1, ISEED, LENJ, RWORK( KD+2-LENJ ) );
               for (I = KD + 2 - LENJ; I <= KD + 1; I++) { // 380
                  AB[I][J] = AB( I, J )*( TLEFT+RWORK( I )*TSCAL );
               } // 380
            } // 390
         } else {
            for (J = 1; J <= N; J++) { // 410
               LENJ = min( N-J+1, KD+1 );
               zlarnv(5, ISEED, LENJ, AB( 1, J ) );
               dlarnv(1, ISEED, LENJ, RWORK );
               for (I = 1; I <= LENJ; I++) { // 400
                  AB[I][J] = AB( I, J )*( TLEFT+RWORK( I )*TSCAL );
               } // 400
            } // 410
         }
         zlarnv(2, ISEED, N, B );
         zdscal(N, TWO, B, 1 );
      }

      // Flip the matrix if the transpose will be used.

      if ( !lsame( TRANS, 'N' ) ) {
         if ( UPPER ) {
            for (J = 1; J <= N / 2; J++) { // 420
               LENJ = min( N-2*J+1, KD+1 );
               zswap(LENJ, AB( KD+1, J ), LDAB-1, AB( KD+2-LENJ, N-J+1 ), -1 );
            } // 420
         } else {
            for (J = 1; J <= N / 2; J++) { // 430
               LENJ = min( N-2*J+1, KD+1 );
               zswap(LENJ, AB( 1, J ), 1, AB( LENJ, N-J+2-LENJ ), -LDAB+1 );
            } // 430
         }
      }

      }
