      SUBROUTINE DLAGSY( N, K, D, A, LDA, ISEED, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDA, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             A( LDA, * ), D( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, HALF;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             ALPHA, TAU, WA, WB, WN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DGEMV, DGER, DLARNV, DSCAL, DSYMV, DSYR2, XERBLA
      // ..
      // .. External Functions ..
      double             DDOT, DNRM2;
      // EXTERNAL DDOT, DNRM2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -1
      } else if ( K.LT.0 .OR. K.GT.N-1 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      }
      if ( INFO.LT.0 ) {
         xerbla('DLAGSY', -INFO );
         RETURN
      }

      // initialize lower triangle of A to diagonal matrix

      for (J = 1; J <= N; J++) { // 20
         DO 10 I = J + 1, N
            A( I, J ) = ZERO
         } // 10
      } // 20
      for (I = 1; I <= N; I++) { // 30
         A( I, I ) = D( I )
      } // 30

      // Generate lower triangle of symmetric matrix

      DO 40 I = N - 1, 1, -1

         // generate random reflection

         dlarnv(3, ISEED, N-I+1, WORK );
         WN = DNRM2( N-I+1, WORK, 1 )
         WA = SIGN( WN, WORK( 1 ) )
         if ( WN.EQ.ZERO ) {
            TAU = ZERO
         } else {
            WB = WORK( 1 ) + WA
            dscal(N-I, ONE / WB, WORK( 2 ), 1 );
            WORK( 1 ) = ONE
            TAU = WB / WA
         }

         // apply random reflection to A(i:n,i:n) from the left
         // and the right

         // compute  y := tau * A * u

         dsymv('Lower', N-I+1, TAU, A( I, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 );

         // compute  v := y - 1/2 * tau * ( y, u ) * u

         ALPHA = -HALF*TAU*DDOT( N-I+1, WORK( N+1 ), 1, WORK, 1 )
         daxpy(N-I+1, ALPHA, WORK, 1, WORK( N+1 ), 1 );

         // apply the transformation as a rank-2 update to A(i:n,i:n)

         dsyr2('Lower', N-I+1, -ONE, WORK, 1, WORK( N+1 ), 1, A( I, I ), LDA );
      } // 40

      // Reduce number of subdiagonals to K

      DO 60 I = 1, N - 1 - K

         // generate reflection to annihilate A(k+i+1:n,i)

         WN = DNRM2( N-K-I+1, A( K+I, I ), 1 )
         WA = SIGN( WN, A( K+I, I ) )
         if ( WN.EQ.ZERO ) {
            TAU = ZERO
         } else {
            WB = A( K+I, I ) + WA
            dscal(N-K-I, ONE / WB, A( K+I+1, I ), 1 );
            A( K+I, I ) = ONE
            TAU = WB / WA
         }

         // apply reflection to A(k+i:n,i+1:k+i-1) from the left

         dgemv('Transpose', N-K-I+1, K-1, ONE, A( K+I, I+1 ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 )          CALL DGER( N-K-I+1, K-1, -TAU, A( K+I, I ), 1, WORK, 1, A( K+I, I+1 ), LDA );

         // apply reflection to A(k+i:n,k+i:n) from the left and the right

         // compute  y := tau * A * u

         dsymv('Lower', N-K-I+1, TAU, A( K+I, K+I ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 );

         // compute  v := y - 1/2 * tau * ( y, u ) * u

         ALPHA = -HALF*TAU*DDOT( N-K-I+1, WORK, 1, A( K+I, I ), 1 )
         daxpy(N-K-I+1, ALPHA, A( K+I, I ), 1, WORK, 1 );

         // apply symmetric rank-2 update to A(k+i:n,k+i:n)

         dsyr2('Lower', N-K-I+1, -ONE, A( K+I, I ), 1, WORK, 1, A( K+I, K+I ), LDA );

         A( K+I, I ) = -WA
         DO 50 J = K + I + 1, N
            A( J, I ) = ZERO
         } // 50
      } // 60

      // Store full symmetric matrix

      for (J = 1; J <= N; J++) { // 80
         DO 70 I = J + 1, N
            A( J, I ) = A( I, J )
         } // 70
      } // 80
      RETURN

      // End of DLAGSY

      }
