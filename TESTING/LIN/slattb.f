      SUBROUTINE SLATTB( IMAT, UPLO, TRANS, DIAG, ISEED, N, KD, AB, LDAB, B, WORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                IMAT, INFO, KD, LDAB, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      REAL               AB( LDAB, * ), B( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, TWO, ZERO
      PARAMETER          ( ONE = 1.0E+0, TWO = 2.0E+0, ZERO = 0.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      String             DIST, PACKIT, TYPE;
      String             PATH;
      int                I, IOFF, IY, J, JCOUNT, KL, KU, LENJ, MODE;
      REAL               ANORM, BIGNUM, BNORM, BSCAL, CNDNUM, PLUS1, PLUS2, REXP, SFAC, SMLNUM, STAR1, TEXP, TLEFT, TNORM, TSCAL, ULP, UNFL
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SLAMCH, SLARND
      // EXTERNAL LSAME, ISAMAX, SLAMCH, SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLARNV, SLATB4, SLATMS, SSCAL, SSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'TB'
      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      SMLNUM = UNFL
      BIGNUM = ( ONE-ULP ) / SMLNUM
      IF( ( IMAT.GE.6 .AND. IMAT.LE.9 ) .OR. IMAT.EQ.17 ) THEN
         DIAG = 'U'
      ELSE
         DIAG = 'N'
      END IF
      INFO = 0

      // Quick return if N.LE.0.

      IF( N.LE.0 ) RETURN

      // Call SLATB4 to set parameters for SLATMS.

      UPPER = LSAME( UPLO, 'U' )
      IF( UPPER ) THEN
         CALL SLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
         KU = KD
         IOFF = 1 + MAX( 0, KD-N+1 )
         KL = 0
         PACKIT = 'Q'
      ELSE
         CALL SLATB4( PATH, -IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
         KL = KD
         IOFF = 1
         KU = 0
         PACKIT = 'B'
      END IF

      // IMAT <= 5:  Non-unit triangular matrix

      IF( IMAT.LE.5 ) THEN
         CALL SLATMS( N, N, DIST, ISEED, TYPE, B, MODE, CNDNUM, ANORM, KL, KU, PACKIT, AB( IOFF, 1 ), LDAB, WORK, INFO )

      // IMAT > 5:  Unit triangular matrix
      // The diagonal is deliberately set to something other than 1.

      // IMAT = 6:  Matrix is the identity

      ELSE IF( IMAT.EQ.6 ) THEN
         IF( UPPER ) THEN
            DO 20 J = 1, N
               DO 10 I = MAX( 1, KD+2-J ), KD
                  AB( I, J ) = ZERO
   10          CONTINUE
               AB( KD+1, J ) = J
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               AB( 1, J ) = J
               DO 30 I = 2, MIN( KD+1, N-J+1 )
                  AB( I, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
         END IF

      // IMAT > 6:  Non-trivial unit triangular matrix

      // A unit triangular matrix T with condition CNDNUM is formed.
      // In this version, T only has bandwidth 2, the rest of it is zero.

      ELSE IF( IMAT.LE.9 ) THEN
         TNORM = SQRT( CNDNUM )

         // Initialize AB to zero.

         IF( UPPER ) THEN
            DO 60 J = 1, N
               DO 50 I = MAX( 1, KD+2-J ), KD
                  AB( I, J ) = ZERO
   50          CONTINUE
               AB( KD+1, J ) = REAL( J )
   60       CONTINUE
         ELSE
            DO 80 J = 1, N
               DO 70 I = 2, MIN( KD+1, N-J+1 )
                  AB( I, J ) = ZERO
   70          CONTINUE
               AB( 1, J ) = REAL( J )
   80       CONTINUE
         END IF

         // Special case:  T is tridiagonal.  Set every other offdiagonal
         // so that the matrix has norm TNORM+1.

         IF( KD.EQ.1 ) THEN
            IF( UPPER ) THEN
               AB( 1, 2 ) = SIGN( TNORM, SLARND( 2, ISEED ) )
               LENJ = ( N-3 ) / 2
               CALL SLARNV( 2, ISEED, LENJ, WORK )
               DO 90 J = 1, LENJ
                  AB( 1, 2*( J+1 ) ) = TNORM*WORK( J )
   90          CONTINUE
            ELSE
               AB( 2, 1 ) = SIGN( TNORM, SLARND( 2, ISEED ) )
               LENJ = ( N-3 ) / 2
               CALL SLARNV( 2, ISEED, LENJ, WORK )
               DO 100 J = 1, LENJ
                  AB( 2, 2*J+1 ) = TNORM*WORK( J )
  100          CONTINUE
            END IF
         ELSE IF( KD.GT.1 ) THEN

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

            STAR1 = SIGN( TNORM, SLARND( 2, ISEED ) )
            SFAC = SQRT( TNORM )
            PLUS1 = SIGN( SFAC, SLARND( 2, ISEED ) )
            DO 110 J = 1, N, 2
               PLUS2 = STAR1 / PLUS1
               WORK( J ) = PLUS1
               WORK( N+J ) = STAR1
               IF( J+1.LE.N ) THEN
                  WORK( J+1 ) = PLUS2
                  WORK( N+J+1 ) = ZERO
                  PLUS1 = STAR1 / PLUS2

                  // Generate a new *-value with norm between sqrt(TNORM)
                  // and TNORM.

                  REXP = SLARND( 2, ISEED )
                  IF( REXP.LT.ZERO ) THEN
                     STAR1 = -SFAC**( ONE-REXP )
                  ELSE
                     STAR1 = SFAC**( ONE+REXP )
                  END IF
               END IF
  110       CONTINUE

            // Copy the tridiagonal T to AB.

            IF( UPPER ) THEN
               CALL SCOPY( N-1, WORK, 1, AB( KD, 2 ), LDAB )
               CALL SCOPY( N-2, WORK( N+1 ), 1, AB( KD-1, 3 ), LDAB )
            ELSE
               CALL SCOPY( N-1, WORK, 1, AB( 2, 1 ), LDAB )
               CALL SCOPY( N-2, WORK( N+1 ), 1, AB( 3, 1 ), LDAB )
            END IF
         END IF

      // IMAT > 9:  Pathological test cases.  These triangular matrices
      // are badly scaled or badly conditioned, so when used in solving a
     t // riangular system they may cause overflow in the solution vector.

      ELSE IF( IMAT.EQ.10 ) THEN

         // Type 10:  Generate a triangular matrix with elements between
         // -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
         // Make the right hand side large so that it requires scaling.

         IF( UPPER ) THEN
            DO 120 J = 1, N
               LENJ = MIN( J, KD+1 )
               CALL SLARNV( 2, ISEED, LENJ, AB( KD+2-LENJ, J ) )
               AB( KD+1, J ) = SIGN( TWO, AB( KD+1, J ) )
  120       CONTINUE
         ELSE
            DO 130 J = 1, N
               LENJ = MIN( N-J+1, KD+1 )
               IF( LENJ.GT.0 ) CALL SLARNV( 2, ISEED, LENJ, AB( 1, J ) )
               AB( 1, J ) = SIGN( TWO, AB( 1, J ) )
  130       CONTINUE
         END IF

         // Set the right hand side so that the largest value is BIGNUM.

         CALL SLARNV( 2, ISEED, N, B )
         IY = ISAMAX( N, B, 1 )
         BNORM = ABS( B( IY ) )
         BSCAL = BIGNUM / MAX( ONE, BNORM )
         CALL SSCAL( N, BSCAL, B, 1 )

      ELSE IF( IMAT.EQ.11 ) THEN

         // Type 11:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 11, the offdiagonal elements are small (CNORM(j) < 1).

         CALL SLARNV( 2, ISEED, N, B )
         TSCAL = ONE / REAL( KD+1 )
         IF( UPPER ) THEN
            DO 140 J = 1, N
               LENJ = MIN( J, KD+1 )
               CALL SLARNV( 2, ISEED, LENJ, AB( KD+2-LENJ, J ) )
               CALL SSCAL( LENJ-1, TSCAL, AB( KD+2-LENJ, J ), 1 )
               AB( KD+1, J ) = SIGN( ONE, AB( KD+1, J ) )
  140       CONTINUE
            AB( KD+1, N ) = SMLNUM*AB( KD+1, N )
         ELSE
            DO 150 J = 1, N
               LENJ = MIN( N-J+1, KD+1 )
               CALL SLARNV( 2, ISEED, LENJ, AB( 1, J ) )
               IF( LENJ.GT.1 ) CALL SSCAL( LENJ-1, TSCAL, AB( 2, J ), 1 )
               AB( 1, J ) = SIGN( ONE, AB( 1, J ) )
  150       CONTINUE
            AB( 1, 1 ) = SMLNUM*AB( 1, 1 )
         END IF

      ELSE IF( IMAT.EQ.12 ) THEN

         // Type 12:  Make the first diagonal element in the solve small to
         // cause immediate overflow when dividing by T(j,j).
         // In type 12, the offdiagonal elements are O(1) (CNORM(j) > 1).

         CALL SLARNV( 2, ISEED, N, B )
         IF( UPPER ) THEN
            DO 160 J = 1, N
               LENJ = MIN( J, KD+1 )
               CALL SLARNV( 2, ISEED, LENJ, AB( KD+2-LENJ, J ) )
               AB( KD+1, J ) = SIGN( ONE, AB( KD+1, J ) )
  160       CONTINUE
            AB( KD+1, N ) = SMLNUM*AB( KD+1, N )
         ELSE
            DO 170 J = 1, N
               LENJ = MIN( N-J+1, KD+1 )
               CALL SLARNV( 2, ISEED, LENJ, AB( 1, J ) )
               AB( 1, J ) = SIGN( ONE, AB( 1, J ) )
  170       CONTINUE
            AB( 1, 1 ) = SMLNUM*AB( 1, 1 )
         END IF

      ELSE IF( IMAT.EQ.13 ) THEN

         // Type 13:  T is diagonal with small numbers on the diagonal to
         // make the growth factor underflow, but a small right hand side
         // chosen so that the solution does not overflow.

         IF( UPPER ) THEN
            JCOUNT = 1
            DO 190 J = N, 1, -1
               DO 180 I = MAX( 1, KD+1-( J-1 ) ), KD
                  AB( I, J ) = ZERO
  180          CONTINUE
               IF( JCOUNT.LE.2 ) THEN
                  AB( KD+1, J ) = SMLNUM
               ELSE
                  AB( KD+1, J ) = ONE
               END IF
               JCOUNT = JCOUNT + 1
               IF( JCOUNT.GT.4 ) JCOUNT = 1
  190       CONTINUE
         ELSE
            JCOUNT = 1
            DO 210 J = 1, N
               DO 200 I = 2, MIN( N-J+1, KD+1 )
                  AB( I, J ) = ZERO
  200          CONTINUE
               IF( JCOUNT.LE.2 ) THEN
                  AB( 1, J ) = SMLNUM
               ELSE
                  AB( 1, J ) = ONE
               END IF
               JCOUNT = JCOUNT + 1
               IF( JCOUNT.GT.4 ) JCOUNT = 1
  210       CONTINUE
         END IF

         // Set the right hand side alternately zero and small.

         IF( UPPER ) THEN
            B( 1 ) = ZERO
            DO 220 I = N, 2, -2
               B( I ) = ZERO
               B( I-1 ) = SMLNUM
  220       CONTINUE
         ELSE
            B( N ) = ZERO
            DO 230 I = 1, N - 1, 2
               B( I ) = ZERO
               B( I+1 ) = SMLNUM
  230       CONTINUE
         END IF

      ELSE IF( IMAT.EQ.14 ) THEN

         // Type 14:  Make the diagonal elements small to cause gradual
         // overflow when dividing by T(j,j).  To control the amount of
         // scaling needed, the matrix is bidiagonal.

         TEXP = ONE / REAL( KD+1 )
         TSCAL = SMLNUM**TEXP
         CALL SLARNV( 2, ISEED, N, B )
         IF( UPPER ) THEN
            DO 250 J = 1, N
               DO 240 I = MAX( 1, KD+2-J ), KD
                  AB( I, J ) = ZERO
  240          CONTINUE
               IF( J.GT.1 .AND. KD.GT.0 ) AB( KD, J ) = -ONE
               AB( KD+1, J ) = TSCAL
  250       CONTINUE
            B( N ) = ONE
         ELSE
            DO 270 J = 1, N
               DO 260 I = 3, MIN( N-J+1, KD+1 )
                  AB( I, J ) = ZERO
  260          CONTINUE
               IF( J.LT.N .AND. KD.GT.0 ) AB( 2, J ) = -ONE
               AB( 1, J ) = TSCAL
  270       CONTINUE
            B( 1 ) = ONE
         END IF

      ELSE IF( IMAT.EQ.15 ) THEN

         // Type 15:  One zero diagonal element.

         IY = N / 2 + 1
         IF( UPPER ) THEN
            DO 280 J = 1, N
               LENJ = MIN( J, KD+1 )
               CALL SLARNV( 2, ISEED, LENJ, AB( KD+2-LENJ, J ) )
               IF( J.NE.IY ) THEN
                  AB( KD+1, J ) = SIGN( TWO, AB( KD+1, J ) )
               ELSE
                  AB( KD+1, J ) = ZERO
               END IF
  280       CONTINUE
         ELSE
            DO 290 J = 1, N
               LENJ = MIN( N-J+1, KD+1 )
               CALL SLARNV( 2, ISEED, LENJ, AB( 1, J ) )
               IF( J.NE.IY ) THEN
                  AB( 1, J ) = SIGN( TWO, AB( 1, J ) )
               ELSE
                  AB( 1, J ) = ZERO
               END IF
  290       CONTINUE
         END IF
         CALL SLARNV( 2, ISEED, N, B )
         CALL SSCAL( N, TWO, B, 1 )

      ELSE IF( IMAT.EQ.16 ) THEN

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
         IF( KD.GT.0 ) THEN
            IF( UPPER ) THEN
               DO 330 J = N, 1, -KD
                  DO 320 I = J, MAX( 1, J-KD+1 ), -2
                     AB( 1+( J-I ), I ) = -TSCAL / REAL( KD+2 )
                     AB( KD+1, I ) = ONE
                     B( I ) = TEXP*( ONE-ULP )
                     IF( I.GT.MAX( 1, J-KD+1 ) ) THEN
                        AB( 2+( J-I ), I-1 ) = -( TSCAL / REAL( KD+2 ) ) / REAL( KD+3 )
                        AB( KD+1, I-1 ) = ONE
                        B( I-1 ) = TEXP*REAL( ( KD+1 )*( KD+1 )+KD )
                     END IF
                     TEXP = TEXP*TWO
  320             CONTINUE
                  B( MAX( 1, J-KD+1 ) ) = ( REAL( KD+2 ) / REAL( KD+3 ) )*TSCAL
  330          CONTINUE
            ELSE
               DO 350 J = 1, N, KD
                  TEXP = ONE
                  LENJ = MIN( KD+1, N-J+1 )
                  DO 340 I = J, MIN( N, J+KD-1 ), 2
                     AB( LENJ-( I-J ), J ) = -TSCAL / REAL( KD+2 )
                     AB( 1, J ) = ONE
                     B( J ) = TEXP*( ONE-ULP )
                     IF( I.LT.MIN( N, J+KD-1 ) ) THEN
                        AB( LENJ-( I-J+1 ), I+1 ) = -( TSCAL / REAL( KD+2 ) ) / REAL( KD+3 )
                        AB( 1, I+1 ) = ONE
                        B( I+1 ) = TEXP*REAL( ( KD+1 )*( KD+1 )+KD )
                     END IF
                     TEXP = TEXP*TWO
  340             CONTINUE
                  B( MIN( N, J+KD-1 ) ) = ( REAL( KD+2 ) / REAL( KD+3 ) )*TSCAL
  350          CONTINUE
            END IF
         ELSE
            DO 360 J = 1, N
               AB( 1, J ) = ONE
               B( J ) = REAL( J )
  360       CONTINUE
         END IF

      ELSE IF( IMAT.EQ.17 ) THEN

         // Type 17:  Generate a unit triangular matrix with elements
         // between -1 and 1, and make the right hand side large so that it
         // requires scaling.

         IF( UPPER ) THEN
            DO 370 J = 1, N
               LENJ = MIN( J-1, KD )
               CALL SLARNV( 2, ISEED, LENJ, AB( KD+1-LENJ, J ) )
               AB( KD+1, J ) = REAL( J )
  370       CONTINUE
         ELSE
            DO 380 J = 1, N
               LENJ = MIN( N-J, KD )
               IF( LENJ.GT.0 ) CALL SLARNV( 2, ISEED, LENJ, AB( 2, J ) )
               AB( 1, J ) = REAL( J )
  380       CONTINUE
         END IF

         // Set the right hand side so that the largest value is BIGNUM.

         CALL SLARNV( 2, ISEED, N, B )
         IY = ISAMAX( N, B, 1 )
         BNORM = ABS( B( IY ) )
         BSCAL = BIGNUM / MAX( ONE, BNORM )
         CALL SSCAL( N, BSCAL, B, 1 )

      ELSE IF( IMAT.EQ.18 ) THEN

         // Type 18:  Generate a triangular matrix with elements between
         // BIGNUM/KD and BIGNUM so that at least one of the column
         // norms will exceed BIGNUM.

         TLEFT = BIGNUM / MAX( ONE, REAL( KD ) )
         TSCAL = BIGNUM*( REAL( KD ) / REAL( KD+1 ) )
         IF( UPPER ) THEN
            DO 400 J = 1, N
               LENJ = MIN( J, KD+1 )
               CALL SLARNV( 2, ISEED, LENJ, AB( KD+2-LENJ, J ) )
               DO 390 I = KD + 2 - LENJ, KD + 1
                  AB( I, J ) = SIGN( TLEFT, AB( I, J ) ) + TSCAL*AB( I, J )
  390          CONTINUE
  400       CONTINUE
         ELSE
            DO 420 J = 1, N
               LENJ = MIN( N-J+1, KD+1 )
               CALL SLARNV( 2, ISEED, LENJ, AB( 1, J ) )
               DO 410 I = 1, LENJ
                  AB( I, J ) = SIGN( TLEFT, AB( I, J ) ) + TSCAL*AB( I, J )
  410          CONTINUE
  420       CONTINUE
         END IF
         CALL SLARNV( 2, ISEED, N, B )
         CALL SSCAL( N, TWO, B, 1 )
      END IF

      // Flip the matrix if the transpose will be used.

      IF( .NOT.LSAME( TRANS, 'N' ) ) THEN
         IF( UPPER ) THEN
            DO 430 J = 1, N / 2
               LENJ = MIN( N-2*J+1, KD+1 )
               CALL SSWAP( LENJ, AB( KD+1, J ), LDAB-1, AB( KD+2-LENJ, N-J+1 ), -1 )
  430       CONTINUE
         ELSE
            DO 440 J = 1, N / 2
               LENJ = MIN( N-2*J+1, KD+1 )
               CALL SSWAP( LENJ, AB( 1, J ), 1, AB( LENJ, N-J+2-LENJ ), -LDAB+1 )
  440       CONTINUE
         END IF
      END IF

      RETURN

      // End of SLATTB

      END
