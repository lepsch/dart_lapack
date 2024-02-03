      SUBROUTINE ZLATTB( IMAT, UPLO, TRANS, DIAG, ISEED, N, KD, AB, LDAB, B, WORK, RWORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                IMAT, INFO, KD, LDAB, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             RWORK( * );
      COMPLEX*16         AB( LDAB, * ), B( * ), WORK( * )
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
      int                I, IOFF, IY, J, JCOUNT, KL, KU, LENJ, MODE;
      double             ANORM, BIGNUM, BNORM, BSCAL, CNDNUM, REXP, SFAC, SMLNUM, TEXP, TLEFT, TNORM, TSCAL, ULP, UNFL;
      COMPLEX*16         PLUS1, PLUS2, STAR1
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH, DLARND;
      COMPLEX*16         ZLARND
      // EXTERNAL LSAME, IZAMAX, DLAMCH, DLARND, ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARNV, ZCOPY, ZDSCAL, ZLARNV, ZLATB4, ZLATMS, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'TB'
      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      SMLNUM = UNFL
      BIGNUM = ( ONE-ULP ) / SMLNUM
      if ( ( IMAT.GE.6 .AND. IMAT.LE.9 ) .OR. IMAT.EQ.17 ) {
         DIAG = 'U'
      } else {
         DIAG = 'N'
      }
      INFO = 0

      // Quick return if N.LE.0.

      IF( N.LE.0 ) RETURN

      // Call ZLATB4 to set parameters for ZLATMS.

      UPPER = LSAME( UPLO, 'U' )
      if ( UPPER ) {
         CALL ZLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
         KU = KD
         IOFF = 1 + MAX( 0, KD-N+1 )
         KL = 0
         PACKIT = 'Q'
      } else {
         CALL ZLATB4( PATH, -IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
         KL = KD
         IOFF = 1
         KU = 0
         PACKIT = 'B'
      }

      // IMAT <= 5:  Non-unit triangular matrix

      if ( IMAT.LE.5 ) {
         CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, AB( IOFF, 1 ), LDAB, WORK, INFO )

      // IMAT > 5:  Unit triangular matrix
      // The diagonal is deliberately set to something other than 1.

      // IMAT = 6:  Matrix is the identity

      } else if ( IMAT.EQ.6 ) {
         if ( UPPER ) {
            DO 20 J = 1, N
               DO 10 I = MAX( 1, KD+2-J ), KD
                  AB( I, J ) = ZERO
   10          CONTINUE
               AB( KD+1, J ) = J
   20       CONTINUE
         } else {
            DO 40 J = 1, N
               AB( 1, J ) = J
               DO 30 I = 2, MIN( KD+1, N-J+1 )
                  AB( I, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
         }

      // IMAT > 6:  Non-trivial unit triangular matrix

      // A unit triangular matrix T with condition CNDNUM is formed.
      // In this version, T only has bandwidth 2, the rest of it is zero.

      } else if ( IMAT.LE.9 ) {
         TNORM = SQRT( CNDNUM )

         // Initialize AB to zero.

         if ( UPPER ) {
            DO 60 J = 1, N
               DO 50 I = MAX( 1, KD+2-J ), KD
                  AB( I, J ) = ZERO
   50          CONTINUE
               AB( KD+1, J ) = DBLE( J )
   60       CONTINUE
         } else {
            DO 80 J = 1, N
               DO 70 I = 2, MIN( KD+1, N-J+1 )
                  AB( I, J ) = ZERO
   70          CONTINUE
               AB( 1, J ) = DBLE( J )
   80       CONTINUE
         }

         // Special case:  T is tridiagonal.  Set every other offdiagonal
         // so that the matrix has norm TNORM+1.

         if ( KD.EQ.1 ) {
            if ( UPPER ) {
               AB( 1, 2 ) = TNORM*ZLARND( 5, ISEED )
               LENJ = ( N-3 ) / 2
               CALL ZLARNV( 2, ISEED, LENJ, WORK )
               DO 90 J = 1, LENJ
                  AB( 1, 2*( J+1 ) ) = TNORM*WORK( J )
   90          CONTINUE
            } else {
               AB( 2, 1 ) = TNORM*ZLARND( 5, ISEED )
               LENJ = ( N-3 ) / 2
               CALL ZLARNV( 2, ISEED, LENJ, WORK )
               DO 100 J = 1, LENJ
                  AB( 2, 2*J+1 ) = TNORM*WORK( J )
  100          CONTINUE
            }
         } else if ( KD.GT.1 ) {

            // Form a unit triangular matrix T with condition CNDNUM.  T is
            // given by
                    // | 1   +   *                      |
                    // |     1   +                      |
                // T = |         1   +   *              |
                    // |             1   +              |
                    // |                 1   +   *      |
                    // |                     1   +      |
                    // |                          . . . |
         // Each element marked with a '*' is formed by taking the product
         // of the adjacent elements marked with '+'.  The '*'s can be
         // chosen freely, and the '+'s are chosen so that the inverse of
         // T will have elements of the same magnitude as T.

         // The two offdiagonals of T are stored in WORK.

            STAR1 = TNORM*ZLARND( 5, ISEED )
            SFAC = SQRT( TNORM )
            PLUS1 = SFAC*ZLARND( 5, ISEED )
            DO 110 J = 1, N, 2
               PLUS2 = STAR1 / PLUS1
               WORK( J ) = PLUS1
               WORK( N+J ) = STAR1
               if ( J+1.LE.N ) {
                  WORK( J+1 ) = PLUS2
                  WORK( N+J+1 ) = ZERO
                  PLUS1 = STAR1 / PLUS2

                  // Generate a new *-value with norm between sqrt(TNORM)
                  // and TNORM.

                  REXP = DLARND( 2, ISEED )
                  if ( REXP.LT.ZERO ) {
                     STAR1 = -SFAC**( ONE-REXP )*ZLARND( 5, ISEED )
                  } else {
                     STAR1 = SFAC**( ONE+REXP )*ZLARND( 5, ISEED )
                  }
               }
  110       CONTINUE

            // Copy the tridiagonal T to AB.

            if ( UPPER ) {
               CALL ZCOPY( N-1, WORK, 1, AB( KD, 2 ), LDAB )
               CALL ZCOPY( N-2, WORK( N+1 ), 1, AB( KD-1, 3 ), LDAB )
            } else {
               CALL ZCOPY( N-1, WORK, 1, AB( 2, 1 ), LDAB )
               CALL ZCOPY( N-2, WORK( N+1 ), 1, AB( 3, 1 ), LDAB )
            }
         }

      // IMAT > 9:  Pathological test cases.  These triangular matrices
      // are badly scaled or badly conditioned, so when used in solving a
     t // riangular system they may cause overflow in the solution vector.

      } else if ( IMAT.EQ.10 ) {

         // Type 10:  Generate a triangular matrix with elements between
         // -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
         // Make the right hand side large so that it requires scaling.

         if ( UPPER ) {
            DO 120 J = 1, N
               LENJ = MIN( J-1, KD )
               CALL ZLARNV( 4, ISEED, LENJ, AB( KD+1-LENJ, J ) )
               AB( KD+1, J ) = ZLARND( 5, ISEED )*TWO
  120       CONTINUE
         } else {
            DO 130 J = 1, N
               LENJ = MIN( N-J, KD )
               IF( LENJ.GT.0 ) CALL ZLARNV( 4, ISEED, LENJ, AB( 2, J ) )
               AB( 1, J ) = ZLARND( 5, ISEED )*TWO
  130       CONTINUE
         }

         // Set the right hand side so that the largest value is BIGNUM.

         CALL ZLARNV( 2, ISEED, N, B )
         IY = IZAMAX( N, B, 1 )
         BNORM = ABS( B( IY ) )
         BSCAL = BIGNUM / MAX( ONE, BNORM )
         CALL ZDSCAL( N, BSCAL, B, 1 )

      } else if ( IMAT.EQ.11 ) {

         // Type 11:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 11, the offdiagonal elements are small (CNORM(j) < 1).

         CALL ZLARNV( 2, ISEED, N, B )
         TSCAL = ONE / DBLE( KD+1 )
         if ( UPPER ) {
            DO 140 J = 1, N
               LENJ = MIN( J-1, KD )
               if ( LENJ.GT.0 ) {
                  CALL ZLARNV( 4, ISEED, LENJ, AB( KD+2-LENJ, J ) )
                  CALL ZDSCAL( LENJ, TSCAL, AB( KD+2-LENJ, J ), 1 )
               }
               AB( KD+1, J ) = ZLARND( 5, ISEED )
  140       CONTINUE
            AB( KD+1, N ) = SMLNUM*AB( KD+1, N )
         } else {
            DO 150 J = 1, N
               LENJ = MIN( N-J, KD )
               if ( LENJ.GT.0 ) {
                  CALL ZLARNV( 4, ISEED, LENJ, AB( 2, J ) )
                  CALL ZDSCAL( LENJ, TSCAL, AB( 2, J ), 1 )
               }
               AB( 1, J ) = ZLARND( 5, ISEED )
  150       CONTINUE
            AB( 1, 1 ) = SMLNUM*AB( 1, 1 )
         }

      } else if ( IMAT.EQ.12 ) {

         // Type 12:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 12, the offdiagonal elements are O(1) (CNORM(j) > 1).

         CALL ZLARNV( 2, ISEED, N, B )
         if ( UPPER ) {
            DO 160 J = 1, N
               LENJ = MIN( J-1, KD )
               IF( LENJ.GT.0 ) CALL ZLARNV( 4, ISEED, LENJ, AB( KD+2-LENJ, J ) )
               AB( KD+1, J ) = ZLARND( 5, ISEED )
  160       CONTINUE
            AB( KD+1, N ) = SMLNUM*AB( KD+1, N )
         } else {
            DO 170 J = 1, N
               LENJ = MIN( N-J, KD )
               IF( LENJ.GT.0 ) CALL ZLARNV( 4, ISEED, LENJ, AB( 2, J ) )
               AB( 1, J ) = ZLARND( 5, ISEED )
  170       CONTINUE
            AB( 1, 1 ) = SMLNUM*AB( 1, 1 )
         }

      } else if ( IMAT.EQ.13 ) {

         // Type 13:  T is diagonal with small numbers on the diagonal to
         // make the growth factor underflow, but a small right hand side
         // chosen so that the solution does not overflow.

         if ( UPPER ) {
            JCOUNT = 1
            DO 190 J = N, 1, -1
               DO 180 I = MAX( 1, KD+1-( J-1 ) ), KD
                  AB( I, J ) = ZERO
  180          CONTINUE
               if ( JCOUNT.LE.2 ) {
                  AB( KD+1, J ) = SMLNUM*ZLARND( 5, ISEED )
               } else {
                  AB( KD+1, J ) = ZLARND( 5, ISEED )
               }
               JCOUNT = JCOUNT + 1
               IF( JCOUNT.GT.4 ) JCOUNT = 1
  190       CONTINUE
         } else {
            JCOUNT = 1
            DO 210 J = 1, N
               DO 200 I = 2, MIN( N-J+1, KD+1 )
                  AB( I, J ) = ZERO
  200          CONTINUE
               if ( JCOUNT.LE.2 ) {
                  AB( 1, J ) = SMLNUM*ZLARND( 5, ISEED )
               } else {
                  AB( 1, J ) = ZLARND( 5, ISEED )
               }
               JCOUNT = JCOUNT + 1
               IF( JCOUNT.GT.4 ) JCOUNT = 1
  210       CONTINUE
         }

         // Set the right hand side alternately zero and small.

         if ( UPPER ) {
            B( 1 ) = ZERO
            DO 220 I = N, 2, -2
               B( I ) = ZERO
               B( I-1 ) = SMLNUM*ZLARND( 5, ISEED )
  220       CONTINUE
         } else {
            B( N ) = ZERO
            DO 230 I = 1, N - 1, 2
               B( I ) = ZERO
               B( I+1 ) = SMLNUM*ZLARND( 5, ISEED )
  230       CONTINUE
         }

      } else if ( IMAT.EQ.14 ) {

         // Type 14:  Make the diagonal elements small to cause gradual
         // overflow when dividing by T(j,j).  To control the amount of
         // scaling needed, the matrix is bidiagonal.

         TEXP = ONE / DBLE( KD+1 )
         TSCAL = SMLNUM**TEXP
         CALL ZLARNV( 4, ISEED, N, B )
         if ( UPPER ) {
            DO 250 J = 1, N
               DO 240 I = MAX( 1, KD+2-J ), KD
                  AB( I, J ) = ZERO
  240          CONTINUE
               IF( J.GT.1 .AND. KD.GT.0 ) AB( KD, J ) = DCMPLX( -ONE, -ONE )
               AB( KD+1, J ) = TSCAL*ZLARND( 5, ISEED )
  250       CONTINUE
            B( N ) = DCMPLX( ONE, ONE )
         } else {
            DO 270 J = 1, N
               DO 260 I = 3, MIN( N-J+1, KD+1 )
                  AB( I, J ) = ZERO
  260          CONTINUE
               IF( J.LT.N .AND. KD.GT.0 ) AB( 2, J ) = DCMPLX( -ONE, -ONE )
               AB( 1, J ) = TSCAL*ZLARND( 5, ISEED )
  270       CONTINUE
            B( 1 ) = DCMPLX( ONE, ONE )
         }

      } else if ( IMAT.EQ.15 ) {

         // Type 15:  One zero diagonal element.

         IY = N / 2 + 1
         if ( UPPER ) {
            DO 280 J = 1, N
               LENJ = MIN( J, KD+1 )
               CALL ZLARNV( 4, ISEED, LENJ, AB( KD+2-LENJ, J ) )
               if ( J.NE.IY ) {
                  AB( KD+1, J ) = ZLARND( 5, ISEED )*TWO
               } else {
                  AB( KD+1, J ) = ZERO
               }
  280       CONTINUE
         } else {
            DO 290 J = 1, N
               LENJ = MIN( N-J+1, KD+1 )
               CALL ZLARNV( 4, ISEED, LENJ, AB( 1, J ) )
               if ( J.NE.IY ) {
                  AB( 1, J ) = ZLARND( 5, ISEED )*TWO
               } else {
                  AB( 1, J ) = ZERO
               }
  290       CONTINUE
         }
         CALL ZLARNV( 2, ISEED, N, B )
         CALL ZDSCAL( N, TWO, B, 1 )

      } else if ( IMAT.EQ.16 ) {

         // Type 16:  Make the offdiagonal elements large to cause overflow
         // when adding a column of T.  In the non-transposed case, the
         // matrix is constructed to cause overflow when adding a column in
         // every other step.

         TSCAL = UNFL / ULP
         TSCAL = ( ONE-ULP ) / TSCAL
         DO 310 J = 1, N
            DO 300 I = 1, KD + 1
               AB( I, J ) = ZERO
  300       CONTINUE
  310    CONTINUE
         TEXP = ONE
         if ( KD.GT.0 ) {
            if ( UPPER ) {
               DO 330 J = N, 1, -KD
                  DO 320 I = J, MAX( 1, J-KD+1 ), -2
                     AB( 1+( J-I ), I ) = -TSCAL / DBLE( KD+2 )
                     AB( KD+1, I ) = ONE
                     B( I ) = TEXP*( ONE-ULP )
                     if ( I.GT.MAX( 1, J-KD+1 ) ) {
                        AB( 2+( J-I ), I-1 ) = -( TSCAL / DBLE( KD+2 ) ) / DBLE( KD+3 )
                        AB( KD+1, I-1 ) = ONE
                        B( I-1 ) = TEXP*DBLE( ( KD+1 )*( KD+1 )+KD )
                     }
                     TEXP = TEXP*TWO
  320             CONTINUE
                  B( MAX( 1, J-KD+1 ) ) = ( DBLE( KD+2 ) / DBLE( KD+3 ) )*TSCAL
  330          CONTINUE
            } else {
               DO 350 J = 1, N, KD
                  TEXP = ONE
                  LENJ = MIN( KD+1, N-J+1 )
                  DO 340 I = J, MIN( N, J+KD-1 ), 2
                     AB( LENJ-( I-J ), J ) = -TSCAL / DBLE( KD+2 )
                     AB( 1, J ) = ONE
                     B( J ) = TEXP*( ONE-ULP )
                     if ( I.LT.MIN( N, J+KD-1 ) ) {
                        AB( LENJ-( I-J+1 ), I+1 ) = -( TSCAL / DBLE( KD+2 ) ) / DBLE( KD+3 )
                        AB( 1, I+1 ) = ONE
                        B( I+1 ) = TEXP*DBLE( ( KD+1 )*( KD+1 )+KD )
                     }
                     TEXP = TEXP*TWO
  340             CONTINUE
                  B( MIN( N, J+KD-1 ) ) = ( DBLE( KD+2 ) / DBLE( KD+3 ) )*TSCAL
  350          CONTINUE
            }
         }

      } else if ( IMAT.EQ.17 ) {

         // Type 17:  Generate a unit triangular matrix with elements
         // between -1 and 1, and make the right hand side large so that it
         // requires scaling.

         if ( UPPER ) {
            DO 360 J = 1, N
               LENJ = MIN( J-1, KD )
               CALL ZLARNV( 4, ISEED, LENJ, AB( KD+1-LENJ, J ) )
               AB( KD+1, J ) = DBLE( J )
  360       CONTINUE
         } else {
            DO 370 J = 1, N
               LENJ = MIN( N-J, KD )
               IF( LENJ.GT.0 ) CALL ZLARNV( 4, ISEED, LENJ, AB( 2, J ) )
               AB( 1, J ) = DBLE( J )
  370       CONTINUE
         }

         // Set the right hand side so that the largest value is BIGNUM.

         CALL ZLARNV( 2, ISEED, N, B )
         IY = IZAMAX( N, B, 1 )
         BNORM = ABS( B( IY ) )
         BSCAL = BIGNUM / MAX( ONE, BNORM )
         CALL ZDSCAL( N, BSCAL, B, 1 )

      } else if ( IMAT.EQ.18 ) {

         // Type 18:  Generate a triangular matrix with elements between
         // BIGNUM/(KD+1) and BIGNUM so that at least one of the column
         // norms will exceed BIGNUM.
         // 1/3/91:  ZLATBS no longer can handle this case

         TLEFT = BIGNUM / DBLE( KD+1 )
         TSCAL = BIGNUM*( DBLE( KD+1 ) / DBLE( KD+2 ) )
         if ( UPPER ) {
            DO 390 J = 1, N
               LENJ = MIN( J, KD+1 )
               CALL ZLARNV( 5, ISEED, LENJ, AB( KD+2-LENJ, J ) )
               CALL DLARNV( 1, ISEED, LENJ, RWORK( KD+2-LENJ ) )
               DO 380 I = KD + 2 - LENJ, KD + 1
                  AB( I, J ) = AB( I, J )*( TLEFT+RWORK( I )*TSCAL )
  380          CONTINUE
  390       CONTINUE
         } else {
            DO 410 J = 1, N
               LENJ = MIN( N-J+1, KD+1 )
               CALL ZLARNV( 5, ISEED, LENJ, AB( 1, J ) )
               CALL DLARNV( 1, ISEED, LENJ, RWORK )
               DO 400 I = 1, LENJ
                  AB( I, J ) = AB( I, J )*( TLEFT+RWORK( I )*TSCAL )
  400          CONTINUE
  410       CONTINUE
         }
         CALL ZLARNV( 2, ISEED, N, B )
         CALL ZDSCAL( N, TWO, B, 1 )
      }

      // Flip the matrix if the transpose will be used.

      if ( .NOT.LSAME( TRANS, 'N' ) ) {
         if ( UPPER ) {
            DO 420 J = 1, N / 2
               LENJ = MIN( N-2*J+1, KD+1 )
               CALL ZSWAP( LENJ, AB( KD+1, J ), LDAB-1, AB( KD+2-LENJ, N-J+1 ), -1 )
  420       CONTINUE
         } else {
            DO 430 J = 1, N / 2
               LENJ = MIN( N-2*J+1, KD+1 )
               CALL ZSWAP( LENJ, AB( 1, J ), 1, AB( LENJ, N-J+2-LENJ ), -LDAB+1 )
  430       CONTINUE
         }
      }

      RETURN

      // End of ZLATTB

      }
