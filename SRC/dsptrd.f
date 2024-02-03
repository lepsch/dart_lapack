      SUBROUTINE DSPTRD( UPLO, N, AP, D, E, TAU, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      double             AP( * ), D( * ), E( * ), TAU( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO, HALF;
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, HALF = 1.0D0 / 2.0D0 )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, I1, I1I1, II;
      double             ALPHA, TAUI;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DLARFG, DSPMV, DSPR2, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DDOT;
      // EXTERNAL LSAME, DDOT
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSPTRD', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.LE.0 ) RETURN

      IF( UPPER ) THEN

         // Reduce the upper triangle of A.
         // I1 is the index in AP of A(1,I+1).

         I1 = N*( N-1 ) / 2 + 1
         DO 10 I = N - 1, 1, -1

            // Generate elementary reflector H(i) = I - tau * v * v**T
           t // o annihilate A(1:i-1,i+1)

            CALL DLARFG( I, AP( I1+I-1 ), AP( I1 ), 1, TAUI )
            E( I ) = AP( I1+I-1 )

            IF( TAUI.NE.ZERO ) THEN

               // Apply H(i) from both sides to A(1:i,1:i)

               AP( I1+I-1 ) = ONE

               // Compute  y := tau * A * v  storing y in TAU(1:i)

               CALL DSPMV( UPLO, I, TAUI, AP, AP( I1 ), 1, ZERO, TAU, 1 )

               // Compute  w := y - 1/2 * tau * (y**T *v) * v

               ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, AP( I1 ), 1 )
               CALL DAXPY( I, ALPHA, AP( I1 ), 1, TAU, 1 )

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**T - w * v**T

               CALL DSPR2( UPLO, I, -ONE, AP( I1 ), 1, TAU, 1, AP )

               AP( I1+I-1 ) = E( I )
            END IF
            D( I+1 ) = AP( I1+I )
            TAU( I ) = TAUI
            I1 = I1 - I
   10    CONTINUE
         D( 1 ) = AP( 1 )
      ELSE

         // Reduce the lower triangle of A. II is the index in AP of
         // A(i,i) and I1I1 is the index of A(i+1,i+1).

         II = 1
         DO 20 I = 1, N - 1
            I1I1 = II + N - I + 1

            // Generate elementary reflector H(i) = I - tau * v * v**T
           t // o annihilate A(i+2:n,i)

            CALL DLARFG( N-I, AP( II+1 ), AP( II+2 ), 1, TAUI )
            E( I ) = AP( II+1 )

            IF( TAUI.NE.ZERO ) THEN

               // Apply H(i) from both sides to A(i+1:n,i+1:n)

               AP( II+1 ) = ONE

               // Compute  y := tau * A * v  storing y in TAU(i:n-1)

               CALL DSPMV( UPLO, N-I, TAUI, AP( I1I1 ), AP( II+1 ), 1, ZERO, TAU( I ), 1 )

               // Compute  w := y - 1/2 * tau * (y**T *v) * v

               ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, AP( II+1 ), 1 )
               CALL DAXPY( N-I, ALPHA, AP( II+1 ), 1, TAU( I ), 1 )

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**T - w * v**T

               CALL DSPR2( UPLO, N-I, -ONE, AP( II+1 ), 1, TAU( I ), 1, AP( I1I1 ) )

               AP( II+1 ) = E( I )
            END IF
            D( I ) = AP( II )
            TAU( I ) = TAUI
            II = I1I1
   20    CONTINUE
         D( N ) = AP( II )
      END IF

      RETURN

      // End of DSPTRD

      END
