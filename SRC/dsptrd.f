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
      const              ONE = 1.0D0, ZERO = 0.0D0, HALF = 1.0D0 / 2.0D0 ;
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
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      }
      if ( INFO.NE.0 ) {
         xerbla('DSPTRD', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.LE.0 ) RETURN

      if ( UPPER ) {

         // Reduce the upper triangle of A.
         // I1 is the index in AP of A(1,I+1).

         I1 = N*( N-1 ) / 2 + 1
         DO 10 I = N - 1, 1, -1

            // Generate elementary reflector H(i) = I - tau * v * v**T
            // to annihilate A(1:i-1,i+1)

            dlarfg(I, AP( I1+I-1 ), AP( I1 ), 1, TAUI );
            E( I ) = AP( I1+I-1 )

            if ( TAUI.NE.ZERO ) {

               // Apply H(i) from both sides to A(1:i,1:i)

               AP( I1+I-1 ) = ONE

               // Compute  y := tau * A * v  storing y in TAU(1:i)

               dspmv(UPLO, I, TAUI, AP, AP( I1 ), 1, ZERO, TAU, 1 );

               // Compute  w := y - 1/2 * tau * (y**T *v) * v

               ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, AP( I1 ), 1 )
               daxpy(I, ALPHA, AP( I1 ), 1, TAU, 1 );

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**T - w * v**T

               dspr2(UPLO, I, -ONE, AP( I1 ), 1, TAU, 1, AP );

               AP( I1+I-1 ) = E( I )
            }
            D( I+1 ) = AP( I1+I )
            TAU( I ) = TAUI
            I1 = I1 - I
   10    CONTINUE
         D( 1 ) = AP( 1 )
      } else {

         // Reduce the lower triangle of A. II is the index in AP of
         // A(i,i) and I1I1 is the index of A(i+1,i+1).

         II = 1
         DO 20 I = 1, N - 1
            I1I1 = II + N - I + 1

            // Generate elementary reflector H(i) = I - tau * v * v**T
            // to annihilate A(i+2:n,i)

            dlarfg(N-I, AP( II+1 ), AP( II+2 ), 1, TAUI );
            E( I ) = AP( II+1 )

            if ( TAUI.NE.ZERO ) {

               // Apply H(i) from both sides to A(i+1:n,i+1:n)

               AP( II+1 ) = ONE

               // Compute  y := tau * A * v  storing y in TAU(i:n-1)

               dspmv(UPLO, N-I, TAUI, AP( I1I1 ), AP( II+1 ), 1, ZERO, TAU( I ), 1 );

               // Compute  w := y - 1/2 * tau * (y**T *v) * v

               ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, AP( II+1 ), 1 )
               daxpy(N-I, ALPHA, AP( II+1 ), 1, TAU( I ), 1 );

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**T - w * v**T

               dspr2(UPLO, N-I, -ONE, AP( II+1 ), 1, TAU( I ), 1, AP( I1I1 ) );

               AP( II+1 ) = E( I )
            }
            D( I ) = AP( II )
            TAU( I ) = TAUI
            II = I1I1
   20    CONTINUE
         D( N ) = AP( II )
      }

      RETURN

      // End of DSPTRD

      }
