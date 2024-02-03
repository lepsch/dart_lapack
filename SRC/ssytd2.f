      SUBROUTINE SSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), D( * ), E( * ), TAU( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0, ZERO = 0.0, HALF = 1.0 / 2.0 )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I;
      REAL               ALPHA, TAUI
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SLARFG, SSYMV, SSYR2, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SDOT
      // EXTERNAL LSAME, SDOT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYTD2', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.LE.0 ) RETURN

      IF( UPPER ) THEN

         // Reduce the upper triangle of A

         DO 10 I = N - 1, 1, -1

            // Generate elementary reflector H(i) = I - tau * v * v**T
           t // o annihilate A(1:i-1,i+1)

            CALL SLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI )
            E( I ) = A( I, I+1 )

            IF( TAUI.NE.ZERO ) THEN

               // Apply H(i) from both sides to A(1:i,1:i)

               A( I, I+1 ) = ONE

               // Compute  x := tau * A * v  storing x in TAU(1:i)

               CALL SSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO, TAU, 1 )

               // Compute  w := x - 1/2 * tau * (x**T * v) * v

               ALPHA = -HALF*TAUI*SDOT( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL SAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**T - w * v**T

               CALL SSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A, LDA )

               A( I, I+1 ) = E( I )
            END IF
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE

         // Reduce the lower triangle of A

         DO 20 I = 1, N - 1

            // Generate elementary reflector H(i) = I - tau * v * v**T
           t // o annihilate A(i+2:n,i)

            CALL SLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1, TAUI )
            E( I ) = A( I+1, I )

            IF( TAUI.NE.ZERO ) THEN

               // Apply H(i) from both sides to A(i+1:n,i+1:n)

               A( I+1, I ) = ONE

               // Compute  x := tau * A * v  storing y in TAU(i:n-1)

               CALL SSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, TAU( I ), 1 )

               // Compute  w := x - 1/2 * tau * (x**T * v) * v

               ALPHA = -HALF*TAUI*SDOT( N-I, TAU( I ), 1, A( I+1, I ), 1 )
               CALL SAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**T - w * v**T

               CALL SSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1, A( I+1, I+1 ), LDA )

               A( I+1, I ) = E( I )
            END IF
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF

      RETURN

      // End of SSYTD2

      }
