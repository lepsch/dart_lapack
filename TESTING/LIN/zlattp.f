      SUBROUTINE ZLATTP( IMAT, UPLO, TRANS, DIAG, ISEED, N, AP, B, WORK, RWORK, INFO );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                IMAT, INFO, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             RWORK( * );
      COMPLEX*16         AP( * ), B( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, TWO, ZERO;
      const              ONE = 1.0, TWO = 2.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      String             DIST, PACKIT, TYPE;
      String             PATH;
      int                I, IY, J, JC, JCNEXT, JCOUNT, JJ, JL, JR, JX, KL, KU, MODE;
      double             ANORM, BIGNUM, BNORM, BSCAL, C, CNDNUM, REXP, SFAC, SMLNUM, T, TEXP, TLEFT, TSCAL, ULP, UNFL, X, Y, Z;
      COMPLEX*16         CTEMP, PLUS1, PLUS2, RA, RB, S, STAR1;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH;
      COMPLEX*16         ZLARND;
      // EXTERNAL LSAME, IZAMAX, DLAMCH, ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARNV, ZDSCAL, ZLARNV, ZLATB4, ZLATMS, ZROT, ZROTG
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, MAX, SQRT
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'Zomplex precision';
      PATH( 2: 3 ) = 'TP';
      UNFL = DLAMCH( 'Safe minimum' );
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' );
      SMLNUM = UNFL;
      BIGNUM = ( ONE-ULP ) / SMLNUM;
      if ( ( IMAT >= 7 && IMAT <= 10 ) || IMAT == 18 ) {
         DIAG = 'U';
      } else {
         DIAG = 'N';
      }
      INFO = 0;

      // Quick return if N <= 0.

      if (N <= 0) RETURN;

      // Call ZLATB4 to set parameters for CLATMS.

      UPPER = LSAME( UPLO, 'U' );
      if ( UPPER ) {
         zlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
         PACKIT = 'C';
      } else {
         zlatb4(PATH, -IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
         PACKIT = 'R';
      }

      // IMAT <= 6:  Non-unit triangular matrix

      if ( IMAT <= 6 ) {
         zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, AP, N, WORK, INFO );

      // IMAT > 6:  Unit triangular matrix
      // The diagonal is deliberately set to something other than 1.

      // IMAT = 7:  Matrix is the identity

      } else if ( IMAT == 7 ) {
         if ( UPPER ) {
            JC = 1;
            for (J = 1; J <= N; J++) { // 20
               for (I = 1; I <= J - 1; I++) { // 10
                  AP( JC+I-1 ) = ZERO;
               } // 10
               AP( JC+J-1 ) = J;
               JC = JC + J;
            } // 20
         } else {
            JC = 1;
            for (J = 1; J <= N; J++) { // 40
               AP( JC ) = J;
               for (I = J + 1; I <= N; I++) { // 30
                  AP( JC+I-J ) = ZERO;
               } // 30
               JC = JC + N - J + 1;
            } // 40
         }

      // IMAT > 7:  Non-trivial unit triangular matrix

      // Generate a unit triangular matrix T with condition CNDNUM by
      // forming a triangular matrix with known singular values and
      // filling in the zero entries with Givens rotations.

      } else if ( IMAT <= 10 ) {
         if ( UPPER ) {
            JC = 0;
            for (J = 1; J <= N; J++) { // 60
               for (I = 1; I <= J - 1; I++) { // 50
                  AP( JC+I ) = ZERO;
               } // 50
               AP( JC+J ) = J;
               JC = JC + J;
            } // 60
         } else {
            JC = 1;
            for (J = 1; J <= N; J++) { // 80
               AP( JC ) = J;
               for (I = J + 1; I <= N; I++) { // 70
                  AP( JC+I-J ) = ZERO;
               } // 70
               JC = JC + N - J + 1;
            } // 80
         }

         // Since the trace of a unit triangular matrix is 1, the product
         // of its singular values must be 1.  Let s = sqrt(CNDNUM),
         // x = sqrt(s) - 1/sqrt(s), y = sqrt(2/(n-2))*x, and z = x**2.
         // The following triangular matrix has singular values s, 1, 1,
         // ..., 1, 1/s:

         // 1  y  y  y  ...  y  y  z
            // 1  0  0  ...  0  0  y
               // 1  0  ...  0  0  y
                  // .  ...  .  .  .
                      // .   .  .  .
                          // 1  0  y
                             // 1  y
                                // 1

         // To fill in the zeros, we first multiply by a matrix with small
         // condition number of the form

         // 1  0  0  0  0  ...
            // 1  +  *  0  0  ...
               // 1  +  0  0  0
                  // 1  +  *  0  0
                     // 1  +  0  0
                        // ...
                           // 1  +  0
                              // 1  0
                                 // 1

         // Each element marked with a '*' is formed by taking the product
         // of the adjacent elements marked with '+'.  The '*'s can be
         // chosen freely, and the '+'s are chosen so that the inverse of
         // T will have elements of the same magnitude as T.  If the *'s in
         // both T and inv(T) have small magnitude, T is well conditioned.
         // The two offdiagonals of T are stored in WORK.

         // The product of these two matrices has the form

         // 1  y  y  y  y  y  .  y  y  z
            // 1  +  *  0  0  .  0  0  y
               // 1  +  0  0  .  0  0  y
                  // 1  +  *  .  .  .  .
                     // 1  +  .  .  .  .
                        // .  .  .  .  .
                           // .  .  .  .
                              // 1  +  y
                                 // 1  y
                                    // 1

         // Now we multiply by Givens rotations, using the fact that

               // [  c   s ] [  1   w ] [ -c  -s ] =  [  1  -w ]
               // [ -s   c ] [  0   1 ] [  s  -c ]    [  0   1 ]
         // and
               // [ -c  -s ] [  1   0 ] [  c   s ] =  [  1   0 ]
               // [  s  -c ] [  w   1 ] [ -s   c ]    [ -w   1 ]

         // where c = w / sqrt(w**2+4) and s = 2 / sqrt(w**2+4).

         STAR1 = 0.25*ZLARND( 5, ISEED );
         SFAC = 0.5;
         PLUS1 = SFAC*ZLARND( 5, ISEED );
         DO 90 J = 1, N, 2;
            PLUS2 = STAR1 / PLUS1;
            WORK( J ) = PLUS1;
            WORK( N+J ) = STAR1;
            if ( J+1 <= N ) {
               WORK( J+1 ) = PLUS2;
               WORK( N+J+1 ) = ZERO;
               PLUS1 = STAR1 / PLUS2;
               REXP = DBLE( ZLARND( 2, ISEED ) );
               if ( REXP < ZERO ) {
                  STAR1 = -SFAC**( ONE-REXP )*ZLARND( 5, ISEED );
               } else {
                  STAR1 = SFAC**( ONE+REXP )*ZLARND( 5, ISEED );
               }
            }
         } // 90

         X = SQRT( CNDNUM ) - ONE / SQRT( CNDNUM );
         if ( N > 2 ) {
            Y = SQRT( TWO / DBLE( N-2 ) )*X;
         } else {
            Y = ZERO;
         }
         Z = X*X;

         if ( UPPER ) {

            // Set the upper triangle of A with a unit triangular matrix
            // of known condition number.

            JC = 1;
            for (J = 2; J <= N; J++) { // 100
               AP( JC+1 ) = Y;
               if (J > 2) AP( JC+J-1 ) = WORK( J-2 )                IF( J > 3 ) AP( JC+J-2 ) = WORK( N+J-3 );
               JC = JC + J;
            } // 100
            JC = JC - N;
            AP( JC+1 ) = Z;
            for (J = 2; J <= N - 1; J++) { // 110
               AP( JC+J ) = Y;
            } // 110
         } else {

            // Set the lower triangle of A with a unit triangular matrix
            // of known condition number.

            for (I = 2; I <= N - 1; I++) { // 120
               AP( I ) = Y;
            } // 120
            AP( N ) = Z;
            JC = N + 1;
            for (J = 2; J <= N - 1; J++) { // 130
               AP( JC+1 ) = WORK( J-1 );
               if (J < N-1) AP( JC+2 ) = WORK( N+J-1 );
               AP( JC+N-J ) = Y;
               JC = JC + N - J + 1;
            } // 130
         }

         // Fill in the zeros using Givens rotations

         if ( UPPER ) {
            JC = 1;
            for (J = 1; J <= N - 1; J++) { // 150
               JCNEXT = JC + J;
               RA = AP( JCNEXT+J-1 );
               RB = TWO;
               zrotg(RA, RB, C, S );

               // Multiply by [ c  s; -conjg(s)  c] on the left.

               if ( N > J+1 ) {
                  JX = JCNEXT + J;
                  for (I = J + 2; I <= N; I++) { // 140
                     CTEMP = C*AP( JX+J ) + S*AP( JX+J+1 );
                     AP( JX+J+1 ) = -DCONJG( S )*AP( JX+J ) + C*AP( JX+J+1 );
                     AP( JX+J ) = CTEMP;
                     JX = JX + I;
                  } // 140
               }

               // Multiply by [-c -s;  conjg(s) -c] on the right.

               if (J > 1) CALL ZROT( J-1, AP( JCNEXT ), 1, AP( JC ), 1, -C, -S );

               // Negate A(J,J+1).

               AP( JCNEXT+J-1 ) = -AP( JCNEXT+J-1 );
               JC = JCNEXT;
            } // 150
         } else {
            JC = 1;
            for (J = 1; J <= N - 1; J++) { // 170
               JCNEXT = JC + N - J + 1;
               RA = AP( JC+1 );
               RB = TWO;
               zrotg(RA, RB, C, S );
               S = DCONJG( S );

               // Multiply by [ c -s;  conjg(s) c] on the right.

               if (N > J+1) CALL ZROT( N-J-1, AP( JCNEXT+1 ), 1, AP( JC+2 ), 1, C, -S );

               // Multiply by [-c  s; -conjg(s) -c] on the left.

               if ( J > 1 ) {
                  JX = 1;
                  for (I = 1; I <= J - 1; I++) { // 160
                     CTEMP = -C*AP( JX+J-I ) + S*AP( JX+J-I+1 );
                     AP( JX+J-I+1 ) = -DCONJG( S )*AP( JX+J-I ) - C*AP( JX+J-I+1 );
                     AP( JX+J-I ) = CTEMP;
                     JX = JX + N - I + 1;
                  } // 160
               }

               // Negate A(J+1,J).

               AP( JC+1 ) = -AP( JC+1 );
               JC = JCNEXT;
            } // 170
         }

      // IMAT > 10:  Pathological test cases.  These triangular matrices
      // are badly scaled or badly conditioned, so when used in solving a
      // triangular system they may cause overflow in the solution vector.

      } else if ( IMAT == 11 ) {

         // Type 11:  Generate a triangular matrix with elements between
         // -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
         // Make the right hand side large so that it requires scaling.

         if ( UPPER ) {
            JC = 1;
            for (J = 1; J <= N; J++) { // 180
               zlarnv(4, ISEED, J-1, AP( JC ) );
               AP( JC+J-1 ) = ZLARND( 5, ISEED )*TWO;
               JC = JC + J;
            } // 180
         } else {
            JC = 1;
            for (J = 1; J <= N; J++) { // 190
               if (J < N) CALL ZLARNV( 4, ISEED, N-J, AP( JC+1 ) );
               AP( JC ) = ZLARND( 5, ISEED )*TWO;
               JC = JC + N - J + 1;
            } // 190
         }

         // Set the right hand side so that the largest value is BIGNUM.

         zlarnv(2, ISEED, N, B );
         IY = IZAMAX( N, B, 1 );
         BNORM = ABS( B( IY ) );
         BSCAL = BIGNUM / MAX( ONE, BNORM );
         zdscal(N, BSCAL, B, 1 );

      } else if ( IMAT == 12 ) {

         // Type 12:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 12, the offdiagonal elements are small (CNORM(j) < 1).

         zlarnv(2, ISEED, N, B );
         TSCAL = ONE / MAX( ONE, DBLE( N-1 ) );
         if ( UPPER ) {
            JC = 1;
            for (J = 1; J <= N; J++) { // 200
               zlarnv(4, ISEED, J-1, AP( JC ) );
               zdscal(J-1, TSCAL, AP( JC ), 1 );
               AP( JC+J-1 ) = ZLARND( 5, ISEED );
               JC = JC + J;
            } // 200
            AP( N*( N+1 ) / 2 ) = SMLNUM*AP( N*( N+1 ) / 2 );
         } else {
            JC = 1;
            for (J = 1; J <= N; J++) { // 210
               zlarnv(2, ISEED, N-J, AP( JC+1 ) );
               zdscal(N-J, TSCAL, AP( JC+1 ), 1 );
               AP( JC ) = ZLARND( 5, ISEED );
               JC = JC + N - J + 1;
            } // 210
            AP( 1 ) = SMLNUM*AP( 1 );
         }

      } else if ( IMAT == 13 ) {

         // Type 13:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 13, the offdiagonal elements are O(1) (CNORM(j) > 1).

         zlarnv(2, ISEED, N, B );
         if ( UPPER ) {
            JC = 1;
            for (J = 1; J <= N; J++) { // 220
               zlarnv(4, ISEED, J-1, AP( JC ) );
               AP( JC+J-1 ) = ZLARND( 5, ISEED );
               JC = JC + J;
            } // 220
            AP( N*( N+1 ) / 2 ) = SMLNUM*AP( N*( N+1 ) / 2 );
         } else {
            JC = 1;
            for (J = 1; J <= N; J++) { // 230
               zlarnv(4, ISEED, N-J, AP( JC+1 ) );
               AP( JC ) = ZLARND( 5, ISEED );
               JC = JC + N - J + 1;
            } // 230
            AP( 1 ) = SMLNUM*AP( 1 );
         }

      } else if ( IMAT == 14 ) {

         // Type 14:  T is diagonal with small numbers on the diagonal to
         // make the growth factor underflow, but a small right hand side
         // chosen so that the solution does not overflow.

         if ( UPPER ) {
            JCOUNT = 1;
            JC = ( N-1 )*N / 2 + 1;
            DO 250 J = N, 1, -1;
               for (I = 1; I <= J - 1; I++) { // 240
                  AP( JC+I-1 ) = ZERO;
               } // 240
               if ( JCOUNT <= 2 ) {
                  AP( JC+J-1 ) = SMLNUM*ZLARND( 5, ISEED );
               } else {
                  AP( JC+J-1 ) = ZLARND( 5, ISEED );
               }
               JCOUNT = JCOUNT + 1;
               if (JCOUNT > 4) JCOUNT = 1;
               JC = JC - J + 1;
            } // 250
         } else {
            JCOUNT = 1;
            JC = 1;
            for (J = 1; J <= N; J++) { // 270
               for (I = J + 1; I <= N; I++) { // 260
                  AP( JC+I-J ) = ZERO;
               } // 260
               if ( JCOUNT <= 2 ) {
                  AP( JC ) = SMLNUM*ZLARND( 5, ISEED );
               } else {
                  AP( JC ) = ZLARND( 5, ISEED );
               }
               JCOUNT = JCOUNT + 1;
               if (JCOUNT > 4) JCOUNT = 1;
               JC = JC + N - J + 1;
            } // 270
         }

         // Set the right hand side alternately zero and small.

         if ( UPPER ) {
            B( 1 ) = ZERO;
            DO 280 I = N, 2, -2;
               B( I ) = ZERO;
               B( I-1 ) = SMLNUM*ZLARND( 5, ISEED );
            } // 280
         } else {
            B( N ) = ZERO;
            DO 290 I = 1, N - 1, 2;
               B( I ) = ZERO;
               B( I+1 ) = SMLNUM*ZLARND( 5, ISEED );
            } // 290
         }

      } else if ( IMAT == 15 ) {

         // Type 15:  Make the diagonal elements small to cause gradual
         // overflow when dividing by T(j,j).  To control the amount of
         // scaling needed, the matrix is bidiagonal.

         TEXP = ONE / MAX( ONE, DBLE( N-1 ) );
         TSCAL = SMLNUM**TEXP;
         zlarnv(4, ISEED, N, B );
         if ( UPPER ) {
            JC = 1;
            for (J = 1; J <= N; J++) { // 310
               for (I = 1; I <= J - 2; I++) { // 300
                  AP( JC+I-1 ) = ZERO;
               } // 300
               if (J > 1) AP( JC+J-2 ) = DCMPLX( -ONE, -ONE );
               AP( JC+J-1 ) = TSCAL*ZLARND( 5, ISEED );
               JC = JC + J;
            } // 310
            B( N ) = DCMPLX( ONE, ONE );
         } else {
            JC = 1;
            for (J = 1; J <= N; J++) { // 330
               for (I = J + 2; I <= N; I++) { // 320
                  AP( JC+I-J ) = ZERO;
               } // 320
               if (J < N) AP( JC+1 ) = DCMPLX( -ONE, -ONE );
               AP( JC ) = TSCAL*ZLARND( 5, ISEED );
               JC = JC + N - J + 1;
            } // 330
            B( 1 ) = DCMPLX( ONE, ONE );
         }

      } else if ( IMAT == 16 ) {

         // Type 16:  One zero diagonal element.

         IY = N / 2 + 1;
         if ( UPPER ) {
            JC = 1;
            for (J = 1; J <= N; J++) { // 340
               zlarnv(4, ISEED, J, AP( JC ) );
               if ( J != IY ) {
                  AP( JC+J-1 ) = ZLARND( 5, ISEED )*TWO;
               } else {
                  AP( JC+J-1 ) = ZERO;
               }
               JC = JC + J;
            } // 340
         } else {
            JC = 1;
            for (J = 1; J <= N; J++) { // 350
               zlarnv(4, ISEED, N-J+1, AP( JC ) );
               if ( J != IY ) {
                  AP( JC ) = ZLARND( 5, ISEED )*TWO;
               } else {
                  AP( JC ) = ZERO;
               }
               JC = JC + N - J + 1;
            } // 350
         }
         zlarnv(2, ISEED, N, B );
         zdscal(N, TWO, B, 1 );

      } else if ( IMAT == 17 ) {

         // Type 17:  Make the offdiagonal elements large to cause overflow
         // when adding a column of T.  In the non-transposed case, the
         // matrix is constructed to cause overflow when adding a column in
         // every other step.

         TSCAL = UNFL / ULP;
         TSCAL = ( ONE-ULP ) / TSCAL;
         for (J = 1; J <= N*( N+1 ) / 2; J++) { // 360
            AP( J ) = ZERO;
         } // 360
         TEXP = ONE;
         if ( UPPER ) {
            JC = ( N-1 )*N / 2 + 1;
            DO 370 J = N, 2, -2;
               AP( JC ) = -TSCAL / DBLE( N+1 );
               AP( JC+J-1 ) = ONE;
               B( J ) = TEXP*( ONE-ULP );
               JC = JC - J + 1;
               AP( JC ) = -( TSCAL / DBLE( N+1 ) ) / DBLE( N+2 );
               AP( JC+J-2 ) = ONE;
               B( J-1 ) = TEXP*DBLE( N*N+N-1 );
               TEXP = TEXP*TWO;
               JC = JC - J + 2;
            } // 370
            B( 1 ) = ( DBLE( N+1 ) / DBLE( N+2 ) )*TSCAL;
         } else {
            JC = 1;
            DO 380 J = 1, N - 1, 2;
               AP( JC+N-J ) = -TSCAL / DBLE( N+1 );
               AP( JC ) = ONE;
               B( J ) = TEXP*( ONE-ULP );
               JC = JC + N - J + 1;
               AP( JC+N-J-1 ) = -( TSCAL / DBLE( N+1 ) ) / DBLE( N+2 );
               AP( JC ) = ONE;
               B( J+1 ) = TEXP*DBLE( N*N+N-1 );
               TEXP = TEXP*TWO;
               JC = JC + N - J;
            } // 380
            B( N ) = ( DBLE( N+1 ) / DBLE( N+2 ) )*TSCAL;
         }

      } else if ( IMAT == 18 ) {

         // Type 18:  Generate a unit triangular matrix with elements
         // between -1 and 1, and make the right hand side large so that it
         // requires scaling.

         if ( UPPER ) {
            JC = 1;
            for (J = 1; J <= N; J++) { // 390
               zlarnv(4, ISEED, J-1, AP( JC ) );
               AP( JC+J-1 ) = ZERO;
               JC = JC + J;
            } // 390
         } else {
            JC = 1;
            for (J = 1; J <= N; J++) { // 400
               if (J < N) CALL ZLARNV( 4, ISEED, N-J, AP( JC+1 ) );
               AP( JC ) = ZERO;
               JC = JC + N - J + 1;
            } // 400
         }

         // Set the right hand side so that the largest value is BIGNUM.

         zlarnv(2, ISEED, N, B );
         IY = IZAMAX( N, B, 1 );
         BNORM = ABS( B( IY ) );
         BSCAL = BIGNUM / MAX( ONE, BNORM );
         zdscal(N, BSCAL, B, 1 );

      } else if ( IMAT == 19 ) {

         // Type 19:  Generate a triangular matrix with elements between
         // BIGNUM/(n-1) and BIGNUM so that at least one of the column
         // norms will exceed BIGNUM.
         // 1/3/91:  ZLATPS no longer can handle this case

         TLEFT = BIGNUM / MAX( ONE, DBLE( N-1 ) );
         TSCAL = BIGNUM*( DBLE( N-1 ) / MAX( ONE, DBLE( N ) ) );
         if ( UPPER ) {
            JC = 1;
            for (J = 1; J <= N; J++) { // 420
               zlarnv(5, ISEED, J, AP( JC ) );
               dlarnv(1, ISEED, J, RWORK );
               for (I = 1; I <= J; I++) { // 410
                  AP( JC+I-1 ) = AP( JC+I-1 )*( TLEFT+RWORK( I )*TSCAL );
               } // 410
               JC = JC + J;
            } // 420
         } else {
            JC = 1;
            for (J = 1; J <= N; J++) { // 440
               zlarnv(5, ISEED, N-J+1, AP( JC ) );
               dlarnv(1, ISEED, N-J+1, RWORK );
               for (I = J; I <= N; I++) { // 430
                  AP( JC+I-J ) = AP( JC+I-J )* ( TLEFT+RWORK( I-J+1 )*TSCAL );
               } // 430
               JC = JC + N - J + 1;
            } // 440
         }
         zlarnv(2, ISEED, N, B );
         zdscal(N, TWO, B, 1 );
      }

      // Flip the matrix across its counter-diagonal if the transpose will
      // be used.

      if ( !LSAME( TRANS, 'N' ) ) {
         if ( UPPER ) {
            JJ = 1;
            JR = N*( N+1 ) / 2;
            for (J = 1; J <= N / 2; J++) { // 460
               JL = JJ;
               for (I = J; I <= N - J; I++) { // 450
                  T = DBLE( AP( JR-I+J ) );
                  AP( JR-I+J ) = AP( JL );
                  AP( JL ) = T;
                  JL = JL + I;
               } // 450
               JJ = JJ + J + 1;
               JR = JR - ( N-J+1 );
            } // 460
         } else {
            JL = 1;
            JJ = N*( N+1 ) / 2;
            for (J = 1; J <= N / 2; J++) { // 480
               JR = JJ;
               for (I = J; I <= N - J; I++) { // 470
                  T = DBLE( AP( JL+I-J ) );
                  AP( JL+I-J ) = AP( JR );
                  AP( JR ) = T;
                  JR = JR - I;
               } // 470
               JL = JL + N - J + 1;
               JJ = JJ - J - 1;
            } // 480
         }
      }

      RETURN;

      // End of ZLATTP

      }
