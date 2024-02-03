      SUBROUTINE SGEQPF( M, N, A, LDA, JPVT, TAU, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      REAL               A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I, ITEMP, J, MA, MN, PVT;
      REAL               AII, TEMP, TEMP2, TOL3Z
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQR2, SLARF, SLARFG, SORM2R, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               SLAMCH, SNRM2
      // EXTERNAL ISAMAX, SLAMCH, SNRM2
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQPF', -INFO )
         RETURN
      END IF

      MN = MIN( M, N )
      TOL3Z = SQRT(SLAMCH('Epsilon'))

      // Move initial columns up front

      ITEMP = 1
      DO 10 I = 1, N
         IF( JPVT( I ).NE.0 ) THEN
            IF( I.NE.ITEMP ) THEN
               CALL SSWAP( M, A( 1, I ), 1, A( 1, ITEMP ), 1 )
               JPVT( I ) = JPVT( ITEMP )
               JPVT( ITEMP ) = I
            ELSE
               JPVT( I ) = I
            END IF
            ITEMP = ITEMP + 1
         ELSE
            JPVT( I ) = I
         END IF
   10 CONTINUE
      ITEMP = ITEMP - 1

      // Compute the QR factorization and update remaining columns

      IF( ITEMP.GT.0 ) THEN
         MA = MIN( ITEMP, M )
         CALL SGEQR2( M, MA, A, LDA, TAU, WORK, INFO )
         IF( MA.LT.N ) THEN
            CALL SORM2R( 'Left', 'Transpose', M, N-MA, MA, A, LDA, TAU, A( 1, MA+1 ), LDA, WORK, INFO )
         END IF
      END IF

      IF( ITEMP.LT.MN ) THEN

         // Initialize partial column norms. The first n elements of
         // work store the exact column norms.

         DO 20 I = ITEMP + 1, N
            WORK( I ) = SNRM2( M-ITEMP, A( ITEMP+1, I ), 1 )
            WORK( N+I ) = WORK( I )
   20    CONTINUE

         // Compute factorization

         DO 40 I = ITEMP + 1, MN

            // Determine ith pivot column and swap if necessary

            PVT = ( I-1 ) + ISAMAX( N-I+1, WORK( I ), 1 )

            IF( PVT.NE.I ) THEN
               CALL SSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
               ITEMP = JPVT( PVT )
               JPVT( PVT ) = JPVT( I )
               JPVT( I ) = ITEMP
               WORK( PVT ) = WORK( I )
               WORK( N+PVT ) = WORK( N+I )
            END IF

            // Generate elementary reflector H(i)

            IF( I.LT.M ) THEN
               CALL SLARFG( M-I+1, A( I, I ), A( I+1, I ), 1, TAU( I ) )
            ELSE
               CALL SLARFG( 1, A( M, M ), A( M, M ), 1, TAU( M ) )
            END IF

            IF( I.LT.N ) THEN

               // Apply H(i) to A(i:m,i+1:n) from the left

               AII = A( I, I )
               A( I, I ) = ONE
               CALL SLARF( 'LEFT', M-I+1, N-I, A( I, I ), 1, TAU( I ), A( I, I+1 ), LDA, WORK( 2*N+1 ) )
               A( I, I ) = AII
            END IF

            // Update partial column norms

            DO 30 J = I + 1, N
               IF( WORK( J ).NE.ZERO ) THEN

                  // NOTE: The following 4 lines follow from the analysis in
                  // Lapack Working Note 176.

                  TEMP = ABS( A( I, J ) ) / WORK( J )
                  TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                  TEMP2 = TEMP*( WORK( J ) / WORK( N+J ) )**2
                  IF( TEMP2 .LE. TOL3Z ) THEN
                     IF( M-I.GT.0 ) THEN
                        WORK( J ) = SNRM2( M-I, A( I+1, J ), 1 )
                        WORK( N+J ) = WORK( J )
                     ELSE
                        WORK( J ) = ZERO
                        WORK( N+J ) = ZERO
                     END IF
                  ELSE
                     WORK( J ) = WORK( J )*SQRT( TEMP )
                  END IF
               END IF
   30       CONTINUE

   40    CONTINUE
      END IF
      RETURN

      // End of SGEQPF

      }
