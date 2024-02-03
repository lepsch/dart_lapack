      SUBROUTINE DLATTP( IMAT, UPLO, TRANS, DIAG, ISEED, N, A, B, WORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                IMAT, INFO, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             A( * ), B( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, TWO, ZERO;
      const              ONE = 1.0D+0, TWO = 2.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      String             DIST, PACKIT, TYPE;
      String             PATH;
      int                I, IY, J, JC, JCNEXT, JCOUNT, JJ, JL, JR, JX, KL, KU, MODE;
      double             ANORM, BIGNUM, BNORM, BSCAL, C, CNDNUM, PLUS1, PLUS2, RA, RB, REXP, S, SFAC, SMLNUM, STAR1, STEMP, T, TEXP, TLEFT, TSCAL, ULP, UNFL, X, Y, Z;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH, DLARND;
      // EXTERNAL LSAME, IDAMAX, DLAMCH, DLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARNV, DLATB4, DLATMS, DROT, DROTG, DSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'TP'
      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
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

      // Call DLATB4 to set parameters for DLATMS.

      UPPER = LSAME( UPLO, 'U' )
      if ( UPPER ) {
         dlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
         PACKIT = 'C'
      } else {
         dlatb4(PATH, -IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );
         PACKIT = 'R'
      }

      // IMAT <= 6:  Non-unit triangular matrix

      if ( IMAT.LE.6 ) {
         dlatms(N, N, DIST, ISEED, TYPE, B, MODE, CNDNUM, ANORM, KL, KU, PACKIT, A, N, WORK, INFO );

      // IMAT > 6:  Unit triangular matrix
      // The diagonal is deliberately set to something other than 1.

      // IMAT = 7:  Matrix is the identity

      } else if ( IMAT.EQ.7 ) {
         if ( UPPER ) {
            JC = 1
            DO 20 J = 1, N
               DO 10 I = 1, J - 1
                  A( JC+I-1 ) = ZERO
   10          CONTINUE
               A( JC+J-1 ) = J
               JC = JC + J
   20       CONTINUE
         } else {
            JC = 1
            DO 40 J = 1, N
               A( JC ) = J
               DO 30 I = J + 1, N
                  A( JC+I-J ) = ZERO
   30          CONTINUE
               JC = JC + N - J + 1
   40       CONTINUE
         }

      // IMAT > 7:  Non-trivial unit triangular matrix

      // Generate a unit triangular matrix T with condition CNDNUM by
      // forming a triangular matrix with known singular values and
      // filling in the zero entries with Givens rotations.

      } else if ( IMAT.LE.10 ) {
         if ( UPPER ) {
            JC = 0
            DO 60 J = 1, N
               DO 50 I = 1, J - 1
                  A( JC+I ) = ZERO
   50          CONTINUE
               A( JC+J ) = J
               JC = JC + J
   60       CONTINUE
         } else {
            JC = 1
            DO 80 J = 1, N
               A( JC ) = J
               DO 70 I = J + 1, N
                  A( JC+I-J ) = ZERO
   70          CONTINUE
               JC = JC + N - J + 1
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

         STAR1 = 0.25D0
         SFAC = 0.5D0
         PLUS1 = SFAC
         DO 90 J = 1, N, 2
            PLUS2 = STAR1 / PLUS1
            WORK( J ) = PLUS1
            WORK( N+J ) = STAR1
            if ( J+1.LE.N ) {
               WORK( J+1 ) = PLUS2
               WORK( N+J+1 ) = ZERO
               PLUS1 = STAR1 / PLUS2
               REXP = DLARND( 2, ISEED )
               STAR1 = STAR1*( SFAC**REXP )
               if ( REXP.LT.ZERO ) {
                  STAR1 = -SFAC**( ONE-REXP )
               } else {
                  STAR1 = SFAC**( ONE+REXP )
               }
            }
   90    CONTINUE

         X = SQRT( CNDNUM ) - ONE / SQRT( CNDNUM )
         if ( N.GT.2 ) {
            Y = SQRT( TWO / DBLE( N-2 ) )*X
         } else {
            Y = ZERO
         }
         Z = X*X

         if ( UPPER ) {

            // Set the upper triangle of A with a unit triangular matrix
            // of known condition number.

            JC = 1
            DO 100 J = 2, N
               A( JC+1 ) = Y
               IF( J.GT.2 ) A( JC+J-1 ) = WORK( J-2 )                IF( J.GT.3 ) A( JC+J-2 ) = WORK( N+J-3 )
               JC = JC + J
  100       CONTINUE
            JC = JC - N
            A( JC+1 ) = Z
            DO 110 J = 2, N - 1
               A( JC+J ) = Y
  110       CONTINUE
         } else {

            // Set the lower triangle of A with a unit triangular matrix
            // of known condition number.

            DO 120 I = 2, N - 1
               A( I ) = Y
  120       CONTINUE
            A( N ) = Z
            JC = N + 1
            DO 130 J = 2, N - 1
               A( JC+1 ) = WORK( J-1 )
               IF( J.LT.N-1 ) A( JC+2 ) = WORK( N+J-1 )
               A( JC+N-J ) = Y
               JC = JC + N - J + 1
  130       CONTINUE
         }

         // Fill in the zeros using Givens rotations

         if ( UPPER ) {
            JC = 1
            DO 150 J = 1, N - 1
               JCNEXT = JC + J
               RA = A( JCNEXT+J-1 )
               RB = TWO
               drotg(RA, RB, C, S );

               // Multiply by [ c  s; -s  c] on the left.

               if ( N.GT.J+1 ) {
                  JX = JCNEXT + J
                  DO 140 I = J + 2, N
                     STEMP = C*A( JX+J ) + S*A( JX+J+1 )
                     A( JX+J+1 ) = -S*A( JX+J ) + C*A( JX+J+1 )
                     A( JX+J ) = STEMP
                     JX = JX + I
  140             CONTINUE
               }

               // Multiply by [-c -s;  s -c] on the right.

               IF( J.GT.1 ) CALL DROT( J-1, A( JCNEXT ), 1, A( JC ), 1, -C, -S )

               // Negate A(J,J+1).

               A( JCNEXT+J-1 ) = -A( JCNEXT+J-1 )
               JC = JCNEXT
  150       CONTINUE
         } else {
            JC = 1
            DO 170 J = 1, N - 1
               JCNEXT = JC + N - J + 1
               RA = A( JC+1 )
               RB = TWO
               drotg(RA, RB, C, S );

               // Multiply by [ c -s;  s  c] on the right.

               IF( N.GT.J+1 ) CALL DROT( N-J-1, A( JCNEXT+1 ), 1, A( JC+2 ), 1, C, -S )

               // Multiply by [-c  s; -s -c] on the left.

               if ( J.GT.1 ) {
                  JX = 1
                  DO 160 I = 1, J - 1
                     STEMP = -C*A( JX+J-I ) + S*A( JX+J-I+1 )
                     A( JX+J-I+1 ) = -S*A( JX+J-I ) - C*A( JX+J-I+1 )
                     A( JX+J-I ) = STEMP
                     JX = JX + N - I + 1
  160             CONTINUE
               }

               // Negate A(J+1,J).

               A( JC+1 ) = -A( JC+1 )
               JC = JCNEXT
  170       CONTINUE
         }

      // IMAT > 10:  Pathological test cases.  These triangular matrices
      // are badly scaled or badly conditioned, so when used in solving a
      // triangular system they may cause overflow in the solution vector.

      } else if ( IMAT.EQ.11 ) {

         // Type 11:  Generate a triangular matrix with elements between
         // -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
         // Make the right hand side large so that it requires scaling.

         if ( UPPER ) {
            JC = 1
            DO 180 J = 1, N
               dlarnv(2, ISEED, J, A( JC ) );
               A( JC+J-1 ) = SIGN( TWO, A( JC+J-1 ) )
               JC = JC + J
  180       CONTINUE
         } else {
            JC = 1
            DO 190 J = 1, N
               dlarnv(2, ISEED, N-J+1, A( JC ) );
               A( JC ) = SIGN( TWO, A( JC ) )
               JC = JC + N - J + 1
  190       CONTINUE
         }

         // Set the right hand side so that the largest value is BIGNUM.

         dlarnv(2, ISEED, N, B );
         IY = IDAMAX( N, B, 1 )
         BNORM = ABS( B( IY ) )
         BSCAL = BIGNUM / MAX( ONE, BNORM )
         dscal(N, BSCAL, B, 1 );

      } else if ( IMAT.EQ.12 ) {

         // Type 12:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 12, the offdiagonal elements are small (CNORM(j) < 1).

         dlarnv(2, ISEED, N, B );
         TSCAL = ONE / MAX( ONE, DBLE( N-1 ) )
         if ( UPPER ) {
            JC = 1
            DO 200 J = 1, N
               dlarnv(2, ISEED, J-1, A( JC ) );
               dscal(J-1, TSCAL, A( JC ), 1 );
               A( JC+J-1 ) = SIGN( ONE, DLARND( 2, ISEED ) )
               JC = JC + J
  200       CONTINUE
            A( N*( N+1 ) / 2 ) = SMLNUM
         } else {
            JC = 1
            DO 210 J = 1, N
               dlarnv(2, ISEED, N-J, A( JC+1 ) );
               dscal(N-J, TSCAL, A( JC+1 ), 1 );
               A( JC ) = SIGN( ONE, DLARND( 2, ISEED ) )
               JC = JC + N - J + 1
  210       CONTINUE
            A( 1 ) = SMLNUM
         }

      } else if ( IMAT.EQ.13 ) {

         // Type 13:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 13, the offdiagonal elements are O(1) (CNORM(j) > 1).

         dlarnv(2, ISEED, N, B );
         if ( UPPER ) {
            JC = 1
            DO 220 J = 1, N
               dlarnv(2, ISEED, J-1, A( JC ) );
               A( JC+J-1 ) = SIGN( ONE, DLARND( 2, ISEED ) )
               JC = JC + J
  220       CONTINUE
            A( N*( N+1 ) / 2 ) = SMLNUM
         } else {
            JC = 1
            DO 230 J = 1, N
               dlarnv(2, ISEED, N-J, A( JC+1 ) );
               A( JC ) = SIGN( ONE, DLARND( 2, ISEED ) )
               JC = JC + N - J + 1
  230       CONTINUE
            A( 1 ) = SMLNUM
         }

      } else if ( IMAT.EQ.14 ) {

         // Type 14:  T is diagonal with small numbers on the diagonal to
         // make the growth factor underflow, but a small right hand side
         // chosen so that the solution does not overflow.

         if ( UPPER ) {
            JCOUNT = 1
            JC = ( N-1 )*N / 2 + 1
            DO 250 J = N, 1, -1
               DO 240 I = 1, J - 1
                  A( JC+I-1 ) = ZERO
  240          CONTINUE
               if ( JCOUNT.LE.2 ) {
                  A( JC+J-1 ) = SMLNUM
               } else {
                  A( JC+J-1 ) = ONE
               }
               JCOUNT = JCOUNT + 1
               IF( JCOUNT.GT.4 ) JCOUNT = 1
               JC = JC - J + 1
  250       CONTINUE
         } else {
            JCOUNT = 1
            JC = 1
            DO 270 J = 1, N
               DO 260 I = J + 1, N
                  A( JC+I-J ) = ZERO
  260          CONTINUE
               if ( JCOUNT.LE.2 ) {
                  A( JC ) = SMLNUM
               } else {
                  A( JC ) = ONE
               }
               JCOUNT = JCOUNT + 1
               IF( JCOUNT.GT.4 ) JCOUNT = 1
               JC = JC + N - J + 1
  270       CONTINUE
         }

         // Set the right hand side alternately zero and small.

         if ( UPPER ) {
            B( 1 ) = ZERO
            DO 280 I = N, 2, -2
               B( I ) = ZERO
               B( I-1 ) = SMLNUM
  280       CONTINUE
         } else {
            B( N ) = ZERO
            DO 290 I = 1, N - 1, 2
               B( I ) = ZERO
               B( I+1 ) = SMLNUM
  290       CONTINUE
         }

      } else if ( IMAT.EQ.15 ) {

         // Type 15:  Make the diagonal elements small to cause gradual
         // overflow when dividing by T(j,j).  To control the amount of
         // scaling needed, the matrix is bidiagonal.

         TEXP = ONE / MAX( ONE, DBLE( N-1 ) )
         TSCAL = SMLNUM**TEXP
         dlarnv(2, ISEED, N, B );
         if ( UPPER ) {
            JC = 1
            DO 310 J = 1, N
               DO 300 I = 1, J - 2
                  A( JC+I-1 ) = ZERO
  300          CONTINUE
               IF( J.GT.1 ) A( JC+J-2 ) = -ONE
               A( JC+J-1 ) = TSCAL
               JC = JC + J
  310       CONTINUE
            B( N ) = ONE
         } else {
            JC = 1
            DO 330 J = 1, N
               DO 320 I = J + 2, N
                  A( JC+I-J ) = ZERO
  320          CONTINUE
               IF( J.LT.N ) A( JC+1 ) = -ONE
               A( JC ) = TSCAL
               JC = JC + N - J + 1
  330       CONTINUE
            B( 1 ) = ONE
         }

      } else if ( IMAT.EQ.16 ) {

         // Type 16:  One zero diagonal element.

         IY = N / 2 + 1
         if ( UPPER ) {
            JC = 1
            DO 340 J = 1, N
               dlarnv(2, ISEED, J, A( JC ) );
               if ( J.NE.IY ) {
                  A( JC+J-1 ) = SIGN( TWO, A( JC+J-1 ) )
               } else {
                  A( JC+J-1 ) = ZERO
               }
               JC = JC + J
  340       CONTINUE
         } else {
            JC = 1
            DO 350 J = 1, N
               dlarnv(2, ISEED, N-J+1, A( JC ) );
               if ( J.NE.IY ) {
                  A( JC ) = SIGN( TWO, A( JC ) )
               } else {
                  A( JC ) = ZERO
               }
               JC = JC + N - J + 1
  350       CONTINUE
         }
         dlarnv(2, ISEED, N, B );
         dscal(N, TWO, B, 1 );

      } else if ( IMAT.EQ.17 ) {

         // Type 17:  Make the offdiagonal elements large to cause overflow
         // when adding a column of T.  In the non-transposed case, the
         // matrix is constructed to cause overflow when adding a column in
         // every other step.

         TSCAL = UNFL / ULP
         TSCAL = ( ONE-ULP ) / TSCAL
         DO 360 J = 1, N*( N+1 ) / 2
            A( J ) = ZERO
  360    CONTINUE
         TEXP = ONE
         if ( UPPER ) {
            JC = ( N-1 )*N / 2 + 1
            DO 370 J = N, 2, -2
               A( JC ) = -TSCAL / DBLE( N+1 )
               A( JC+J-1 ) = ONE
               B( J ) = TEXP*( ONE-ULP )
               JC = JC - J + 1
               A( JC ) = -( TSCAL / DBLE( N+1 ) ) / DBLE( N+2 )
               A( JC+J-2 ) = ONE
               B( J-1 ) = TEXP*DBLE( N*N+N-1 )
               TEXP = TEXP*TWO
               JC = JC - J + 2
  370       CONTINUE
            B( 1 ) = ( DBLE( N+1 ) / DBLE( N+2 ) )*TSCAL
         } else {
            JC = 1
            DO 380 J = 1, N - 1, 2
               A( JC+N-J ) = -TSCAL / DBLE( N+1 )
               A( JC ) = ONE
               B( J ) = TEXP*( ONE-ULP )
               JC = JC + N - J + 1
               A( JC+N-J-1 ) = -( TSCAL / DBLE( N+1 ) ) / DBLE( N+2 )
               A( JC ) = ONE
               B( J+1 ) = TEXP*DBLE( N*N+N-1 )
               TEXP = TEXP*TWO
               JC = JC + N - J
  380       CONTINUE
            B( N ) = ( DBLE( N+1 ) / DBLE( N+2 ) )*TSCAL
         }

      } else if ( IMAT.EQ.18 ) {

         // Type 18:  Generate a unit triangular matrix with elements
         // between -1 and 1, and make the right hand side large so that it
         // requires scaling.

         if ( UPPER ) {
            JC = 1
            DO 390 J = 1, N
               dlarnv(2, ISEED, J-1, A( JC ) );
               A( JC+J-1 ) = ZERO
               JC = JC + J
  390       CONTINUE
         } else {
            JC = 1
            DO 400 J = 1, N
               IF( J.LT.N ) CALL DLARNV( 2, ISEED, N-J, A( JC+1 ) )
               A( JC ) = ZERO
               JC = JC + N - J + 1
  400       CONTINUE
         }

         // Set the right hand side so that the largest value is BIGNUM.

         dlarnv(2, ISEED, N, B );
         IY = IDAMAX( N, B, 1 )
         BNORM = ABS( B( IY ) )
         BSCAL = BIGNUM / MAX( ONE, BNORM )
         dscal(N, BSCAL, B, 1 );

      } else if ( IMAT.EQ.19 ) {

         // Type 19:  Generate a triangular matrix with elements between
         // BIGNUM/(n-1) and BIGNUM so that at least one of the column
         // norms will exceed BIGNUM.

         TLEFT = BIGNUM / MAX( ONE, DBLE( N-1 ) )
         TSCAL = BIGNUM*( DBLE( N-1 ) / MAX( ONE, DBLE( N ) ) )
         if ( UPPER ) {
            JC = 1
            DO 420 J = 1, N
               dlarnv(2, ISEED, J, A( JC ) );
               DO 410 I = 1, J
                  A( JC+I-1 ) = SIGN( TLEFT, A( JC+I-1 ) ) + TSCAL*A( JC+I-1 )
  410          CONTINUE
               JC = JC + J
  420       CONTINUE
         } else {
            JC = 1
            DO 440 J = 1, N
               dlarnv(2, ISEED, N-J+1, A( JC ) );
               DO 430 I = J, N
                  A( JC+I-J ) = SIGN( TLEFT, A( JC+I-J ) ) + TSCAL*A( JC+I-J )
  430          CONTINUE
               JC = JC + N - J + 1
  440       CONTINUE
         }
         dlarnv(2, ISEED, N, B );
         dscal(N, TWO, B, 1 );
      }

      // Flip the matrix across its counter-diagonal if the transpose will
      // be used.

      if ( .NOT.LSAME( TRANS, 'N' ) ) {
         if ( UPPER ) {
            JJ = 1
            JR = N*( N+1 ) / 2
            DO 460 J = 1, N / 2
               JL = JJ
               DO 450 I = J, N - J
                  T = A( JR-I+J )
                  A( JR-I+J ) = A( JL )
                  A( JL ) = T
                  JL = JL + I
  450          CONTINUE
               JJ = JJ + J + 1
               JR = JR - ( N-J+1 )
  460       CONTINUE
         } else {
            JL = 1
            JJ = N*( N+1 ) / 2
            DO 480 J = 1, N / 2
               JR = JJ
               DO 470 I = J, N - J
                  T = A( JL+I-J )
                  A( JL+I-J ) = A( JR )
                  A( JR ) = T
                  JR = JR - I
  470          CONTINUE
               JL = JL + N - J + 1
               JJ = JJ - J - 1
  480       CONTINUE
         }
      }

      RETURN

      // End of DLATTP

      }
