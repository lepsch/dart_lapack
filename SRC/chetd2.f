      SUBROUTINE CHETD2( UPLO, N, A, LDA, D, E, TAU, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * )
      COMPLEX            A( LDA, * ), TAU( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, ZERO, HALF
      const              ONE = ( 1.0E+0, 0.0E+0 ), ZERO = ( 0.0E+0, 0.0E+0 ), HALF = ( 0.5E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I;
      COMPLEX            ALPHA, TAUI
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CHEMV, CHER2, CLARFG, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      COMPLEX            CDOTC
      // EXTERNAL LSAME, CDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CHETD2', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.LE.0 ) RETURN

      if ( UPPER ) {

         // Reduce the upper triangle of A

         A( N, N ) = REAL( A( N, N ) )
         DO 10 I = N - 1, 1, -1

            // Generate elementary reflector H(i) = I - tau * v * v**H
           t // o annihilate A(1:i-1,i+1)

            ALPHA = A( I, I+1 )
            CALL CLARFG( I, ALPHA, A( 1, I+1 ), 1, TAUI )
            E( I ) = REAL( ALPHA )

            if ( TAUI.NE.ZERO ) {

               // Apply H(i) from both sides to A(1:i,1:i)

               A( I, I+1 ) = ONE

               // Compute  x := tau * A * v  storing x in TAU(1:i)

               CALL CHEMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO, TAU, 1 )

               // Compute  w := x - 1/2 * tau * (x**H * v) * v

               ALPHA = -HALF*TAUI*CDOTC( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL CAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**H - w * v**H

               CALL CHER2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A, LDA )

            } else {
               A( I, I ) = REAL( A( I, I ) )
            }
            A( I, I+1 ) = E( I )
            D( I+1 ) = REAL( A( I+1, I+1 ) )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = REAL( A( 1, 1 ) )
      } else {

         // Reduce the lower triangle of A

         A( 1, 1 ) = REAL( A( 1, 1 ) )
         DO 20 I = 1, N - 1

            // Generate elementary reflector H(i) = I - tau * v * v**H
           t // o annihilate A(i+2:n,i)

            ALPHA = A( I+1, I )
            CALL CLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAUI )
            E( I ) = REAL( ALPHA )

            if ( TAUI.NE.ZERO ) {

               // Apply H(i) from both sides to A(i+1:n,i+1:n)

               A( I+1, I ) = ONE

               // Compute  x := tau * A * v  storing y in TAU(i:n-1)

               CALL CHEMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, TAU( I ), 1 )

               // Compute  w := x - 1/2 * tau * (x**H * v) * v

               ALPHA = -HALF*TAUI*CDOTC( N-I, TAU( I ), 1, A( I+1, I ), 1 )
               CALL CAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**H - w * v**H

               CALL CHER2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1, A( I+1, I+1 ), LDA )

            } else {
               A( I+1, I+1 ) = REAL( A( I+1, I+1 ) )
            }
            A( I+1, I ) = E( I )
            D( I ) = REAL( A( I, I ) )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = REAL( A( N, N ) )
      }

      RETURN

      // End of CHETD2

      }
