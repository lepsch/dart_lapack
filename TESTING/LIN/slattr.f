      SUBROUTINE SLATTR( IMAT, UPLO, TRANS, DIAG, ISEED, N, A, LDA, B, WORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                IMAT, INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      REAL               A( LDA, * ), B( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, TWO, ZERO
      const              ONE = 1.0E+0, TWO = 2.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      String             DIST, TYPE;
      String             PATH;
      int                I, IY, J, JCOUNT, KL, KU, MODE;
      REAL               ANORM, BIGNUM, BNORM, BSCAL, C, CNDNUM, PLUS1, PLUS2, RA, RB, REXP, S, SFAC, SMLNUM, STAR1, TEXP, TLEFT, TSCAL, ULP, UNFL, X, Y, Z
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SLAMCH, SLARND
      // EXTERNAL LSAME, ISAMAX, SLAMCH, SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLARNV, SLATB4, SLATMS, SROT, SROTG, SSCAL, SSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'TR'
      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      SMLNUM = UNFL
      BIGNUM = ( ONE-ULP ) / SMLNUM
      if ( ( IMAT.GE.7 .AND. IMAT.LE.10 ) .OR. IMAT.EQ.18 ) {
         DIAG = 'U'
      } else {
         DIAG = 'N'
      }
      INFO = 0

      // Quick return if N.LE.0.

      IF( N.LE.0 ) RETURN

      // Call SLATB4 to set parameters for SLATMS.

      UPPER = LSAME( UPLO, 'U' )
      if ( UPPER ) {
         slatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
      } else {
         slatb4(PATH, -IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
      }

      // IMAT <= 6:  Non-unit triangular matrix

      if ( IMAT.LE.6 ) {
         slatms(N, N, DIST, ISEED, TYPE, B, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO );

      // IMAT > 6:  Unit triangular matrix
      // The diagonal is deliberately set to something other than 1.

      // IMAT = 7:  Matrix is the identity

      } else if ( IMAT.EQ.7 ) {
         if ( UPPER ) {
            DO 20 J = 1, N
               DO 10 I = 1, J - 1
                  A( I, J ) = ZERO
   10          CONTINUE
               A( J, J ) = J
   20       CONTINUE
         } else {
            DO 40 J = 1, N
               A( J, J ) = J
               DO 30 I = J + 1, N
                  A( I, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
         }

      // IMAT > 7:  Non-trivial unit triangular matrix

      // Generate a unit triangular matrix T with condition CNDNUM by
      // forming a triangular matrix with known singular values and
      // filling in the zero entries with Givens rotations.

      } else if ( IMAT.LE.10 ) {
         if ( UPPER ) {
            DO 60 J = 1, N
               DO 50 I = 1, J - 1
                  A( I, J ) = ZERO
   50          CONTINUE
               A( J, J ) = J
   60       CONTINUE
         } else {
            DO 80 J = 1, N
               A( J, J ) = J
               DO 70 I = J + 1, N
                  A( I, J ) = ZERO
   70          CONTINUE
   80       CONTINUE
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

         STAR1 = 0.25
         SFAC = 0.5
         PLUS1 = SFAC
         DO 90 J = 1, N, 2
            PLUS2 = STAR1 / PLUS1
            WORK( J ) = PLUS1
            WORK( N+J ) = STAR1
            if ( J+1.LE.N ) {
               WORK( J+1 ) = PLUS2
               WORK( N+J+1 ) = ZERO
               PLUS1 = STAR1 / PLUS2
               REXP = SLARND( 2, ISEED )
               STAR1 = STAR1*( SFAC**REXP )
               if ( REXP.LT.ZERO ) {
                  STAR1 = -SFAC**( ONE-REXP )
               } else {
                  STAR1 = SFAC**( ONE+REXP )
               }
            }
   90    CONTINUE

         X = SQRT( CNDNUM ) - 1 / SQRT( CNDNUM )
         if ( N.GT.2 ) {
            Y = SQRT( 2. / ( N-2 ) )*X
         } else {
            Y = ZERO
         }
         Z = X*X

         if ( UPPER ) {
            if ( N.GT.3 ) {
               scopy(N-3, WORK, 1, A( 2, 3 ), LDA+1 );
               IF( N.GT.4 ) CALL SCOPY( N-4, WORK( N+1 ), 1, A( 2, 4 ), LDA+1 )
            }
            DO 100 J = 2, N - 1
               A( 1, J ) = Y
               A( J, N ) = Y
  100       CONTINUE
            A( 1, N ) = Z
         } else {
            if ( N.GT.3 ) {
               scopy(N-3, WORK, 1, A( 3, 2 ), LDA+1 );
               IF( N.GT.4 ) CALL SCOPY( N-4, WORK( N+1 ), 1, A( 4, 2 ), LDA+1 )
            }
            DO 110 J = 2, N - 1
               A( J, 1 ) = Y
               A( N, J ) = Y
  110       CONTINUE
            A( N, 1 ) = Z
         }

         // Fill in the zeros using Givens rotations.

         if ( UPPER ) {
            DO 120 J = 1, N - 1
               RA = A( J, J+1 )
               RB = 2.0
               srotg(RA, RB, C, S );

               // Multiply by [ c  s; -s  c] on the left.

               IF( N.GT.J+1 ) CALL SROT( N-J-1, A( J, J+2 ), LDA, A( J+1, J+2 ), LDA, C, S )

               // Multiply by [-c -s;  s -c] on the right.

               IF( J.GT.1 ) CALL SROT( J-1, A( 1, J+1 ), 1, A( 1, J ), 1, -C, -S )

               // Negate A(J,J+1).

               A( J, J+1 ) = -A( J, J+1 )
  120       CONTINUE
         } else {
            DO 130 J = 1, N - 1
               RA = A( J+1, J )
               RB = 2.0
               srotg(RA, RB, C, S );

               // Multiply by [ c -s;  s  c] on the right.

               IF( N.GT.J+1 ) CALL SROT( N-J-1, A( J+2, J+1 ), 1, A( J+2, J ), 1, C, -S )

               // Multiply by [-c  s; -s -c] on the left.

               IF( J.GT.1 ) CALL SROT( J-1, A( J, 1 ), LDA, A( J+1, 1 ), LDA, -C, S )

               // Negate A(J+1,J).

               A( J+1, J ) = -A( J+1, J )
  130       CONTINUE
         }

      // IMAT > 10:  Pathological test cases.  These triangular matrices
      // are badly scaled or badly conditioned, so when used in solving a
      // triangular system they may cause overflow in the solution vector.

      } else if ( IMAT.EQ.11 ) {

         // Type 11:  Generate a triangular matrix with elements between
         // -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
         // Make the right hand side large so that it requires scaling.

         if ( UPPER ) {
            DO 140 J = 1, N
               slarnv(2, ISEED, J, A( 1, J ) );
               A( J, J ) = SIGN( TWO, A( J, J ) )
  140       CONTINUE
         } else {
            DO 150 J = 1, N
               slarnv(2, ISEED, N-J+1, A( J, J ) );
               A( J, J ) = SIGN( TWO, A( J, J ) )
  150       CONTINUE
         }

         // Set the right hand side so that the largest value is BIGNUM.

         slarnv(2, ISEED, N, B );
         IY = ISAMAX( N, B, 1 )
         BNORM = ABS( B( IY ) )
         BSCAL = BIGNUM / MAX( ONE, BNORM )
         sscal(N, BSCAL, B, 1 );

      } else if ( IMAT.EQ.12 ) {

         // Type 12:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 12, the offdiagonal elements are small (CNORM(j) < 1).

         slarnv(2, ISEED, N, B );
         TSCAL = ONE / MAX( ONE, REAL( N-1 ) )
         if ( UPPER ) {
            DO 160 J = 1, N
               slarnv(2, ISEED, J, A( 1, J ) );
               sscal(J-1, TSCAL, A( 1, J ), 1 );
               A( J, J ) = SIGN( ONE, A( J, J ) )
  160       CONTINUE
            A( N, N ) = SMLNUM*A( N, N )
         } else {
            DO 170 J = 1, N
               slarnv(2, ISEED, N-J+1, A( J, J ) );
               IF( N.GT.J ) CALL SSCAL( N-J, TSCAL, A( J+1, J ), 1 )
               A( J, J ) = SIGN( ONE, A( J, J ) )
  170       CONTINUE
            A( 1, 1 ) = SMLNUM*A( 1, 1 )
         }

      } else if ( IMAT.EQ.13 ) {

         // Type 13:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 13, the offdiagonal elements are O(1) (CNORM(j) > 1).

         slarnv(2, ISEED, N, B );
         if ( UPPER ) {
            DO 180 J = 1, N
               slarnv(2, ISEED, J, A( 1, J ) );
               A( J, J ) = SIGN( ONE, A( J, J ) )
  180       CONTINUE
            A( N, N ) = SMLNUM*A( N, N )
         } else {
            DO 190 J = 1, N
               slarnv(2, ISEED, N-J+1, A( J, J ) );
               A( J, J ) = SIGN( ONE, A( J, J ) )
  190       CONTINUE
            A( 1, 1 ) = SMLNUM*A( 1, 1 )
         }

      } else if ( IMAT.EQ.14 ) {

         // Type 14:  T is diagonal with small numbers on the diagonal to
         // make the growth factor underflow, but a small right hand side
         // chosen so that the solution does not overflow.

         if ( UPPER ) {
            JCOUNT = 1
            DO 210 J = N, 1, -1
               DO 200 I = 1, J - 1
                  A( I, J ) = ZERO
  200          CONTINUE
               if ( JCOUNT.LE.2 ) {
                  A( J, J ) = SMLNUM
               } else {
                  A( J, J ) = ONE
               }
               JCOUNT = JCOUNT + 1
               IF( JCOUNT.GT.4 ) JCOUNT = 1
  210       CONTINUE
         } else {
            JCOUNT = 1
            DO 230 J = 1, N
               DO 220 I = J + 1, N
                  A( I, J ) = ZERO
  220          CONTINUE
               if ( JCOUNT.LE.2 ) {
                  A( J, J ) = SMLNUM
               } else {
                  A( J, J ) = ONE
               }
               JCOUNT = JCOUNT + 1
               IF( JCOUNT.GT.4 ) JCOUNT = 1
  230       CONTINUE
         }

         // Set the right hand side alternately zero and small.

         if ( UPPER ) {
            B( 1 ) = ZERO
            DO 240 I = N, 2, -2
               B( I ) = ZERO
               B( I-1 ) = SMLNUM
  240       CONTINUE
         } else {
            B( N ) = ZERO
            DO 250 I = 1, N - 1, 2
               B( I ) = ZERO
               B( I+1 ) = SMLNUM
  250       CONTINUE
         }

      } else if ( IMAT.EQ.15 ) {

         // Type 15:  Make the diagonal elements small to cause gradual
         // overflow when dividing by T(j,j).  To control the amount of
         // scaling needed, the matrix is bidiagonal.

         TEXP = ONE / MAX( ONE, REAL( N-1 ) )
         TSCAL = SMLNUM**TEXP
         slarnv(2, ISEED, N, B );
         if ( UPPER ) {
            DO 270 J = 1, N
               DO 260 I = 1, J - 2
                  A( I, J ) = 0.
  260          CONTINUE
               IF( J.GT.1 ) A( J-1, J ) = -ONE
               A( J, J ) = TSCAL
  270       CONTINUE
            B( N ) = ONE
         } else {
            DO 290 J = 1, N
               DO 280 I = J + 2, N
                  A( I, J ) = 0.
  280          CONTINUE
               IF( J.LT.N ) A( J+1, J ) = -ONE
               A( J, J ) = TSCAL
  290       CONTINUE
            B( 1 ) = ONE
         }

      } else if ( IMAT.EQ.16 ) {

         // Type 16:  One zero diagonal element.

         IY = N / 2 + 1
         if ( UPPER ) {
            DO 300 J = 1, N
               slarnv(2, ISEED, J, A( 1, J ) );
               if ( J.NE.IY ) {
                  A( J, J ) = SIGN( TWO, A( J, J ) )
               } else {
                  A( J, J ) = ZERO
               }
  300       CONTINUE
         } else {
            DO 310 J = 1, N
               slarnv(2, ISEED, N-J+1, A( J, J ) );
               if ( J.NE.IY ) {
                  A( J, J ) = SIGN( TWO, A( J, J ) )
               } else {
                  A( J, J ) = ZERO
               }
  310       CONTINUE
         }
         slarnv(2, ISEED, N, B );
         sscal(N, TWO, B, 1 );

      } else if ( IMAT.EQ.17 ) {

         // Type 17:  Make the offdiagonal elements large to cause overflow
         // when adding a column of T.  In the non-transposed case, the
         // matrix is constructed to cause overflow when adding a column in
         // every other step.

         TSCAL = UNFL / ULP
         TSCAL = ( ONE-ULP ) / TSCAL
         DO 330 J = 1, N
            DO 320 I = 1, N
               A( I, J ) = 0.
  320       CONTINUE
  330    CONTINUE
         TEXP = ONE
         if ( UPPER ) {
            DO 340 J = N, 2, -2
               A( 1, J ) = -TSCAL / REAL( N+1 )
               A( J, J ) = ONE
               B( J ) = TEXP*( ONE-ULP )
               A( 1, J-1 ) = -( TSCAL / REAL( N+1 ) ) / REAL( N+2 )
               A( J-1, J-1 ) = ONE
               B( J-1 ) = TEXP*REAL( N*N+N-1 )
               TEXP = TEXP*2.
  340       CONTINUE
            B( 1 ) = ( REAL( N+1 ) / REAL( N+2 ) )*TSCAL
         } else {
            DO 350 J = 1, N - 1, 2
               A( N, J ) = -TSCAL / REAL( N+1 )
               A( J, J ) = ONE
               B( J ) = TEXP*( ONE-ULP )
               A( N, J+1 ) = -( TSCAL / REAL( N+1 ) ) / REAL( N+2 )
               A( J+1, J+1 ) = ONE
               B( J+1 ) = TEXP*REAL( N*N+N-1 )
               TEXP = TEXP*2.
  350       CONTINUE
            B( N ) = ( REAL( N+1 ) / REAL( N+2 ) )*TSCAL
         }

      } else if ( IMAT.EQ.18 ) {

         // Type 18:  Generate a unit triangular matrix with elements
         // between -1 and 1, and make the right hand side large so that it
         // requires scaling.

         if ( UPPER ) {
            DO 360 J = 1, N
               slarnv(2, ISEED, J-1, A( 1, J ) );
               A( J, J ) = ZERO
  360       CONTINUE
         } else {
            DO 370 J = 1, N
               IF( J.LT.N ) CALL SLARNV( 2, ISEED, N-J, A( J+1, J ) )
               A( J, J ) = ZERO
  370       CONTINUE
         }

         // Set the right hand side so that the largest value is BIGNUM.

         slarnv(2, ISEED, N, B );
         IY = ISAMAX( N, B, 1 )
         BNORM = ABS( B( IY ) )
         BSCAL = BIGNUM / MAX( ONE, BNORM )
         sscal(N, BSCAL, B, 1 );

      } else if ( IMAT.EQ.19 ) {

         // Type 19:  Generate a triangular matrix with elements between
         // BIGNUM/(n-1) and BIGNUM so that at least one of the column
         // norms will exceed BIGNUM.
         // 1/3/91:  SLATRS no longer can handle this case

         TLEFT = BIGNUM / MAX( ONE, REAL( N-1 ) )
         TSCAL = BIGNUM*( REAL( N-1 ) / MAX( ONE, REAL( N ) ) )
         if ( UPPER ) {
            DO 390 J = 1, N
               slarnv(2, ISEED, J, A( 1, J ) );
               DO 380 I = 1, J
                  A( I, J ) = SIGN( TLEFT, A( I, J ) ) + TSCAL*A( I, J )
  380          CONTINUE
  390       CONTINUE
         } else {
            DO 410 J = 1, N
               slarnv(2, ISEED, N-J+1, A( J, J ) );
               DO 400 I = J, N
                  A( I, J ) = SIGN( TLEFT, A( I, J ) ) + TSCAL*A( I, J )
  400          CONTINUE
  410       CONTINUE
         }
         slarnv(2, ISEED, N, B );
         sscal(N, TWO, B, 1 );
      }

      // Flip the matrix if the transpose will be used.

      if ( .NOT.LSAME( TRANS, 'N' ) ) {
         if ( UPPER ) {
            DO 420 J = 1, N / 2
               sswap(N-2*J+1, A( J, J ), LDA, A( J+1, N-J+1 ), -1 );
  420       CONTINUE
         } else {
            DO 430 J = 1, N / 2
               sswap(N-2*J+1, A( J, J ), 1, A( N-J+1, J+1 ), -LDA );
  430       CONTINUE
         }
      }

      RETURN

      // End of SLATTR

      }
