// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void slattr(final int IMAT, final int UPLO, final int TRANS, final int DIAG, final Array<int> ISEED_, final int N, final Matrix<double> A_, final int LDA, final int B, final Array<double> _WORK_, final Box<int> INFO,) {
  final ISEED = ISEED_.dim();
  final A = A_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                IMAT, INFO, LDA, N;
      int                ISEED( 4 );
      double               A( LDA, * ), B( * ), WORK( * );
      // ..

      double               ONE, TWO, ZERO;
      const              ONE = 1.0, TWO = 2.0, ZERO = 0.0 ;
      bool               UPPER;
      String             DIST, TYPE;
      String             PATH;
      int                I, IY, J, JCOUNT, KL, KU, MODE;
      double               ANORM, BIGNUM, BNORM, BSCAL, C, CNDNUM, PLUS1, PLUS2, RA, RB, REXP, S, SFAC, SMLNUM, STAR1, TEXP, TLEFT, TSCAL, ULP, UNFL, X, Y, Z;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ISAMAX;
      //- REAL               SLAMCH, SLARND;
      // EXTERNAL lsame, ISAMAX, SLAMCH, SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLARNV, SLATB4, SLATMS, SROT, SROTG, SSCAL, SSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL, SIGN, SQRT

      PATH[1: 1] = 'Single precision';
      PATH[2: 3] = 'TR';
      UNFL = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );
      SMLNUM = UNFL;
      BIGNUM = ( ONE-ULP ) / SMLNUM;
      if ( ( IMAT >= 7 && IMAT <= 10 ) || IMAT == 18 ) {
         DIAG = 'U';
      } else {
         DIAG = 'N';
      }
      INFO = 0;

      // Quick return if N <= 0.

      if (N <= 0) return;

      // Call SLATB4 to set parameters for SLATMS.

      UPPER = lsame( UPLO, 'U' );
      if ( UPPER ) {
         slatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
      } else {
         slatb4(PATH, -IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
      }

      // IMAT <= 6:  Non-unit triangular matrix

      if ( IMAT <= 6 ) {
         slatms(N, N, DIST, ISEED, TYPE, B, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO );

      // IMAT > 6:  Unit triangular matrix
      // The diagonal is deliberately set to something other than 1.

      // IMAT = 7:  Matrix is the identity

      } else if ( IMAT == 7 ) {
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 20
               for (I = 1; I <= J - 1; I++) { // 10
                  A[I][J] = ZERO;
               } // 10
               A[J][J] = J;
            } // 20
         } else {
            for (J = 1; J <= N; J++) { // 40
               A[J][J] = J;
               for (I = J + 1; I <= N; I++) { // 30
                  A[I][J] = ZERO;
               } // 30
            } // 40
         }

      // IMAT > 7:  Non-trivial unit triangular matrix

      // Generate a unit triangular matrix T with condition CNDNUM by
      // forming a triangular matrix with known singular values and
      // filling in the zero entries with Givens rotations.

      } else if ( IMAT <= 10 ) {
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 60
               for (I = 1; I <= J - 1; I++) { // 50
                  A[I][J] = ZERO;
               } // 50
               A[J][J] = J;
            } // 60
         } else {
            for (J = 1; J <= N; J++) { // 80
               A[J][J] = J;
               for (I = J + 1; I <= N; I++) { // 70
                  A[I][J] = ZERO;
               } // 70
            } // 80
         }

         // Since the trace of a unit triangular matrix is 1, the product
         // of its singular values must be 1.  Let s = sqrt(CNDNUM),
         // x = sqrt(s) - 1/sqrt(s), y = sqrt(2/(n-2))*x, and z = x**2.
         // The following triangular matrix has singular values s, 1, 1,
         // ..., 1, 1/s:

         // 1  y  y  y  ...  y  y  z
         //    1  0  0  ...  0  0  y
         //       1  0  ...  0  0  y
         //          .  ...  .  .  .
         //              .   .  .  .
         //                  1  0  y
         //                     1  y
         //                        1

         // To fill in the zeros, we first multiply by a matrix with small
         // condition number of the form

         // 1  0  0  0  0  ...
         //    1  +  *  0  0  ...
         //       1  +  0  0  0
         //          1  +  *  0  0
         //             1  +  0  0
         //                ...
         //                   1  +  0
         //                      1  0
         //                         1

         // Each element marked with a '*' is formed by taking the product
         // of the adjacent elements marked with '+'.  The '*'s can be
         // chosen freely, and the '+'s are chosen so that the inverse of
         // T will have elements of the same magnitude as T.  If the *'s in
         // both T and inv(T) have small magnitude, T is well conditioned.
         // The two offdiagonals of T are stored in WORK.

         // The product of these two matrices has the form

         // 1  y  y  y  y  y  .  y  y  z
         //    1  +  *  0  0  .  0  0  y
         //       1  +  0  0  .  0  0  y
         //          1  +  *  .  .  .  .
         //             1  +  .  .  .  .
         //                .  .  .  .  .
         //                   .  .  .  .
         //                      1  +  y
         //                         1  y
         //                            1

         // Now we multiply by Givens rotations, using the fact that

               // [  c   s ] [  1   w ] [ -c  -s ] =  [  1  -w ]
               // [ -s   c ] [  0   1 ] [  s  -c ]    [  0   1 ]
         // and
         //       [ -c  -s ] [  1   0 ] [  c   s ] =  [  1   0 ]
         //       [  s  -c ] [  w   1 ] [ -s   c ]    [ -w   1 ]

         // where c = w / sqrt(w**2+4) and s = 2 / sqrt(w**2+4).

         STAR1 = 0.25;
         SFAC = 0.5;
         PLUS1 = SFAC;
         for (J = 1; J <= N; J += 2) { // 90
            PLUS2 = STAR1 / PLUS1;
            WORK[J] = PLUS1;
            WORK[N+J] = STAR1;
            if ( J+1 <= N ) {
               WORK[J+1] = PLUS2;
               WORK[N+J+1] = ZERO;
               PLUS1 = STAR1 / PLUS2;
               REXP = SLARND( 2, ISEED );
               STAR1 = STAR1*( SFAC**REXP );
               if ( REXP < ZERO ) {
                  STAR1 = -SFAC**( ONE-REXP );
               } else {
                  STAR1 = SFAC**( ONE+REXP );
               }
            }
         } // 90

         X = sqrt( CNDNUM ) - 1 / sqrt( CNDNUM );
         if ( N > 2 ) {
            Y = sqrt( 2. / ( N-2 ) )*X;
         } else {
            Y = ZERO;
         }
         Z = X*X;

         if ( UPPER ) {
            if ( N > 3 ) {
               scopy(N-3, WORK, 1, A( 2, 3 ), LDA+1 );
               if (N > 4) scopy( N-4, WORK( N+1 ), 1, A( 2, 4 ), LDA+1 );
            }
            for (J = 2; J <= N - 1; J++) { // 100
               A[1][J] = Y;
               A[J][N] = Y;
            } // 100
            A[1][N] = Z;
         } else {
            if ( N > 3 ) {
               scopy(N-3, WORK, 1, A( 3, 2 ), LDA+1 );
               if (N > 4) scopy( N-4, WORK( N+1 ), 1, A( 4, 2 ), LDA+1 );
            }
            for (J = 2; J <= N - 1; J++) { // 110
               A[J][1] = Y;
               A[N][J] = Y;
            } // 110
            A[N][1] = Z;
         }

         // Fill in the zeros using Givens rotations.

         if ( UPPER ) {
            for (J = 1; J <= N - 1; J++) { // 120
               RA = A( J, J+1 );
               RB = 2.0;
               srotg(RA, RB, C, S );

               // Multiply by [ c  s; -s  c] on the left.

               if (N > J+1) srot( N-J-1, A( J, J+2 ), LDA, A( J+1, J+2 ), LDA, C, S );

               // Multiply by [-c -s;  s -c] on the right.

               if (J > 1) srot( J-1, A( 1, J+1 ), 1, A( 1, J ), 1, -C, -S );

               // Negate A(J,J+1).

               A[J][J+1] = -A( J, J+1 );
            } // 120
         } else {
            for (J = 1; J <= N - 1; J++) { // 130
               RA = A( J+1, J );
               RB = 2.0;
               srotg(RA, RB, C, S );

               // Multiply by [ c -s;  s  c] on the right.

               if (N > J+1) srot( N-J-1, A( J+2, J+1 ), 1, A( J+2, J ), 1, C, -S );

               // Multiply by [-c  s; -s -c] on the left.

               if (J > 1) srot( J-1, A( J, 1 ), LDA, A( J+1, 1 ), LDA, -C, S );

               // Negate A(J+1,J).

               A[J+1][J] = -A( J+1, J );
            } // 130
         }

      // IMAT > 10:  Pathological test cases.  These triangular matrices
      // are badly scaled or badly conditioned, so when used in solving a
      // triangular system they may cause overflow in the solution vector.

      } else if ( IMAT == 11 ) {

         // Type 11:  Generate a triangular matrix with elements between
         // -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
         // Make the right hand side large so that it requires scaling.

         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 140
               slarnv(2, ISEED, J, A( 1, J ) );
               A[J][J] = sign( TWO, A( J, J ) );
            } // 140
         } else {
            for (J = 1; J <= N; J++) { // 150
               slarnv(2, ISEED, N-J+1, A( J, J ) );
               A[J][J] = sign( TWO, A( J, J ) );
            } // 150
         }

         // Set the right hand side so that the largest value is BIGNUM.

         slarnv(2, ISEED, N, B );
         IY = ISAMAX( N, B, 1 );
         BNORM = ( B( IY ) ).abs();
         BSCAL = BIGNUM / max( ONE, BNORM );
         sscal(N, BSCAL, B, 1 );

      } else if ( IMAT == 12 ) {

         // Type 12:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 12, the offdiagonal elements are small (CNORM(j) < 1).

         slarnv(2, ISEED, N, B );
         TSCAL = ONE / max( ONE, REAL( N-1 ) );
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 160
               slarnv(2, ISEED, J, A( 1, J ) );
               sscal(J-1, TSCAL, A( 1, J ), 1 );
               A[J][J] = sign( ONE, A( J, J ) );
            } // 160
            A[N][N] = SMLNUM*A( N, N );
         } else {
            for (J = 1; J <= N; J++) { // 170
               slarnv(2, ISEED, N-J+1, A( J, J ) );
               if (N > J) sscal( N-J, TSCAL, A( J+1, J ), 1 );
               A[J][J] = sign( ONE, A( J, J ) );
            } // 170
            A[1][1] = SMLNUM*A( 1, 1 );
         }

      } else if ( IMAT == 13 ) {

         // Type 13:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 13, the offdiagonal elements are O(1) (CNORM(j) > 1).

         slarnv(2, ISEED, N, B );
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 180
               slarnv(2, ISEED, J, A( 1, J ) );
               A[J][J] = sign( ONE, A( J, J ) );
            } // 180
            A[N][N] = SMLNUM*A( N, N );
         } else {
            for (J = 1; J <= N; J++) { // 190
               slarnv(2, ISEED, N-J+1, A( J, J ) );
               A[J][J] = sign( ONE, A( J, J ) );
            } // 190
            A[1][1] = SMLNUM*A( 1, 1 );
         }

      } else if ( IMAT == 14 ) {

         // Type 14:  T is diagonal with small numbers on the diagonal to
         // make the growth factor underflow, but a small right hand side
         // chosen so that the solution does not overflow.

         if ( UPPER ) {
            JCOUNT = 1;
            for (J = N; J >= 1; J--) { // 210
               for (I = 1; I <= J - 1; I++) { // 200
                  A[I][J] = ZERO;
               } // 200
               if ( JCOUNT <= 2 ) {
                  A[J][J] = SMLNUM;
               } else {
                  A[J][J] = ONE;
               }
               JCOUNT = JCOUNT + 1;
               if (JCOUNT > 4) JCOUNT = 1;
            } // 210
         } else {
            JCOUNT = 1;
            for (J = 1; J <= N; J++) { // 230
               for (I = J + 1; I <= N; I++) { // 220
                  A[I][J] = ZERO;
               } // 220
               if ( JCOUNT <= 2 ) {
                  A[J][J] = SMLNUM;
               } else {
                  A[J][J] = ONE;
               }
               JCOUNT = JCOUNT + 1;
               if (JCOUNT > 4) JCOUNT = 1;
            } // 230
         }

         // Set the right hand side alternately zero and small.

         if ( UPPER ) {
            B[1] = ZERO;
            for (I = N; I >= 2; I -= 2) { // 240
               B[I] = ZERO;
               B[I-1] = SMLNUM;
            } // 240
         } else {
            B[N] = ZERO;
            for (I = 1; 2 < 0 ? I >= N - 1 : I <= N - 1; I += 2) { // 250
               B[I] = ZERO;
               B[I+1] = SMLNUM;
            } // 250
         }

      } else if ( IMAT == 15 ) {

         // Type 15:  Make the diagonal elements small to cause gradual
         // overflow when dividing by T(j,j).  To control the amount of
         // scaling needed, the matrix is bidiagonal.

         TEXP = ONE / max( ONE, REAL( N-1 ) );
         TSCAL = SMLNUM**TEXP;
         slarnv(2, ISEED, N, B );
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 270
               for (I = 1; I <= J - 2; I++) { // 260
                  A[I][J] = 0.;
               } // 260
               if (J > 1) A( J-1, J ) = -ONE;
               A[J][J] = TSCAL;
            } // 270
            B[N] = ONE;
         } else {
            for (J = 1; J <= N; J++) { // 290
               for (I = J + 2; I <= N; I++) { // 280
                  A[I][J] = 0.;
               } // 280
               if (J < N) A( J+1, J ) = -ONE;
               A[J][J] = TSCAL;
            } // 290
            B[1] = ONE;
         }

      } else if ( IMAT == 16 ) {

         // Type 16:  One zero diagonal element.

         IY = N / 2 + 1;
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 300
               slarnv(2, ISEED, J, A( 1, J ) );
               if ( J != IY ) {
                  A[J][J] = sign( TWO, A( J, J ) );
               } else {
                  A[J][J] = ZERO;
               }
            } // 300
         } else {
            for (J = 1; J <= N; J++) { // 310
               slarnv(2, ISEED, N-J+1, A( J, J ) );
               if ( J != IY ) {
                  A[J][J] = sign( TWO, A( J, J ) );
               } else {
                  A[J][J] = ZERO;
               }
            } // 310
         }
         slarnv(2, ISEED, N, B );
         sscal(N, TWO, B, 1 );

      } else if ( IMAT == 17 ) {

         // Type 17:  Make the offdiagonal elements large to cause overflow
         // when adding a column of T.  In the non-transposed case, the
         // matrix is constructed to cause overflow when adding a column in
         // every other step.

         TSCAL = UNFL / ULP;
         TSCAL = ( ONE-ULP ) / TSCAL;
         for (J = 1; J <= N; J++) { // 330
            for (I = 1; I <= N; I++) { // 320
               A[I][J] = 0.;
            } // 320
         } // 330
         TEXP = ONE;
         if ( UPPER ) {
            for (J = N; J >= 2; J -= 2) { // 340
               A[1][J] = -TSCAL / REAL( N+1 );
               A[J][J] = ONE;
               B[J] = TEXP*( ONE-ULP );
               A[1][J-1] = -( TSCAL / REAL( N+1 ) ) / REAL( N+2 );
               A[J-1][J-1] = ONE;
               B[J-1] = TEXP*double( N*N+N-1 );
               TEXP = TEXP*2.;
            } // 340
            B[1] = ( double( N+1 ) / REAL( N+2 ) )*TSCAL;
         } else {
            for (J = 1; 2 < 0 ? J >= N - 1 : J <= N - 1; J += 2) { // 350
               A[N][J] = -TSCAL / REAL( N+1 );
               A[J][J] = ONE;
               B[J] = TEXP*( ONE-ULP );
               A[N][J+1] = -( TSCAL / REAL( N+1 ) ) / REAL( N+2 );
               A[J+1][J+1] = ONE;
               B[J+1] = TEXP*double( N*N+N-1 );
               TEXP = TEXP*2.;
            } // 350
            B[N] = ( double( N+1 ) / REAL( N+2 ) )*TSCAL;
         }

      } else if ( IMAT == 18 ) {

         // Type 18:  Generate a unit triangular matrix with elements
         // between -1 and 1, and make the right hand side large so that it
         // requires scaling.

         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 360
               slarnv(2, ISEED, J-1, A( 1, J ) );
               A[J][J] = ZERO;
            } // 360
         } else {
            for (J = 1; J <= N; J++) { // 370
               if (J < N) slarnv( 2, ISEED, N-J, A( J+1, J ) );
               A[J][J] = ZERO;
            } // 370
         }

         // Set the right hand side so that the largest value is BIGNUM.

         slarnv(2, ISEED, N, B );
         IY = ISAMAX( N, B, 1 );
         BNORM = ( B( IY ) ).abs();
         BSCAL = BIGNUM / max( ONE, BNORM );
         sscal(N, BSCAL, B, 1 );

      } else if ( IMAT == 19 ) {

         // Type 19:  Generate a triangular matrix with elements between
         // BIGNUM/(n-1) and BIGNUM so that at least one of the column
         // norms will exceed BIGNUM.
         // 1/3/91:  SLATRS no longer can handle this case

         TLEFT = BIGNUM / max( ONE, REAL( N-1 ) );
         TSCAL = BIGNUM*( double( N-1 ) / max( ONE, REAL( N ) ) );
         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 390
               slarnv(2, ISEED, J, A( 1, J ) );
               for (I = 1; I <= J; I++) { // 380
                  A[I][J] = sign( TLEFT, A( I, J ) ) + TSCAL*A( I, J );
               } // 380
            } // 390
         } else {
            for (J = 1; J <= N; J++) { // 410
               slarnv(2, ISEED, N-J+1, A( J, J ) );
               for (I = J; I <= N; I++) { // 400
                  A[I][J] = sign( TLEFT, A( I, J ) ) + TSCAL*A( I, J );
               } // 400
            } // 410
         }
         slarnv(2, ISEED, N, B );
         sscal(N, TWO, B, 1 );
      }

      // Flip the matrix if the transpose will be used.

      if ( !lsame( TRANS, 'N' ) ) {
         if ( UPPER ) {
            for (J = 1; J <= N / 2; J++) { // 420
               sswap(N-2*J+1, A( J, J ), LDA, A( J+1, N-J+1 ), -1 );
            } // 420
         } else {
            for (J = 1; J <= N / 2; J++) { // 430
               sswap(N-2*J+1, A( J, J ), 1, A( N-J+1, J+1 ), -LDA );
            } // 430
         }
      }

      }
