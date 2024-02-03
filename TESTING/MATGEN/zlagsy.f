      SUBROUTINE ZLAGSY( N, K, D, A, LDA, ISEED, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDA, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             D( * );
      COMPLEX*16         A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO, ONE, HALF
      const              ZERO = ( 0.0D+0, 0.0D+0 ), ONE = ( 1.0D+0, 0.0D+0 ), HALF = ( 0.5D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, II, J, JJ;
      double             WN;
      COMPLEX*16         ALPHA, TAU, WA, WB
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZGEMV, ZGERC, ZLACGV, ZLARNV, ZSCAL, ZSYMV
      // ..
      // .. External Functions ..
      double             DZNRM2;
      COMPLEX*16         ZDOTC
      // EXTERNAL DZNRM2, ZDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX
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
         xerbla('ZLAGSY', -INFO );
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

      DO 60 I = N - 1, 1, -1

         // generate random reflection

         zlarnv(3, ISEED, N-I+1, WORK );
         WN = DZNRM2( N-I+1, WORK, 1 )
         WA = ( WN / ABS( WORK( 1 ) ) )*WORK( 1 )
         if ( WN.EQ.ZERO ) {
            TAU = ZERO
         } else {
            WB = WORK( 1 ) + WA
            zscal(N-I, ONE / WB, WORK( 2 ), 1 );
            WORK( 1 ) = ONE
            TAU = DBLE( WB / WA )
         }

         // apply random reflection to A(i:n,i:n) from the left
         // and the right

         // compute  y := tau * A * conjg(u)

         zlacgv(N-I+1, WORK, 1 );
         zsymv('Lower', N-I+1, TAU, A( I, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 );
         zlacgv(N-I+1, WORK, 1 );

         // compute  v := y - 1/2 * tau * ( u, y ) * u

         ALPHA = -HALF*TAU*ZDOTC( N-I+1, WORK, 1, WORK( N+1 ), 1 )
         zaxpy(N-I+1, ALPHA, WORK, 1, WORK( N+1 ), 1 );

         // apply the transformation as a rank-2 update to A(i:n,i:n)

         // CALL ZSYR2( 'Lower', N-I+1, -ONE, WORK, 1, WORK( N+1 ), 1,
         // $               A( I, I ), LDA )

         for (JJ = I; JJ <= N; JJ++) { // 50
            for (II = JJ; II <= N; II++) { // 40
               A( II, JJ ) = A( II, JJ ) - WORK( II-I+1 )*WORK( N+JJ-I+1 ) - WORK( N+II-I+1 )*WORK( JJ-I+1 )
            } // 40
         } // 50
      } // 60

      // Reduce number of subdiagonals to K

      DO 100 I = 1, N - 1 - K

         // generate reflection to annihilate A(k+i+1:n,i)

         WN = DZNRM2( N-K-I+1, A( K+I, I ), 1 )
         WA = ( WN / ABS( A( K+I, I ) ) )*A( K+I, I )
         if ( WN.EQ.ZERO ) {
            TAU = ZERO
         } else {
            WB = A( K+I, I ) + WA
            zscal(N-K-I, ONE / WB, A( K+I+1, I ), 1 );
            A( K+I, I ) = ONE
            TAU = DBLE( WB / WA )
         }

         // apply reflection to A(k+i:n,i+1:k+i-1) from the left

         zgemv('Conjugate transpose', N-K-I+1, K-1, ONE, A( K+I, I+1 ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 )          CALL ZGERC( N-K-I+1, K-1, -TAU, A( K+I, I ), 1, WORK, 1, A( K+I, I+1 ), LDA );

         // apply reflection to A(k+i:n,k+i:n) from the left and the right

         // compute  y := tau * A * conjg(u)

         zlacgv(N-K-I+1, A( K+I, I ), 1 );
         zsymv('Lower', N-K-I+1, TAU, A( K+I, K+I ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 );
         zlacgv(N-K-I+1, A( K+I, I ), 1 );

         // compute  v := y - 1/2 * tau * ( u, y ) * u

         ALPHA = -HALF*TAU*ZDOTC( N-K-I+1, A( K+I, I ), 1, WORK, 1 )
         zaxpy(N-K-I+1, ALPHA, A( K+I, I ), 1, WORK, 1 );

         // apply symmetric rank-2 update to A(k+i:n,k+i:n)

         // CALL ZSYR2( 'Lower', N-K-I+1, -ONE, A( K+I, I ), 1, WORK, 1,
         // $               A( K+I, K+I ), LDA )

         DO 80 JJ = K + I, N
            for (II = JJ; II <= N; II++) { // 70
               A( II, JJ ) = A( II, JJ ) - A( II, I )*WORK( JJ-K-I+1 ) - WORK( II-K-I+1 )*A( JJ, I )
            } // 70
         } // 80

         A( K+I, I ) = -WA
         DO 90 J = K + I + 1, N
            A( J, I ) = ZERO
         } // 90
      } // 100

      // Store full symmetric matrix

      for (J = 1; J <= N; J++) { // 120
         DO 110 I = J + 1, N
            A( J, I ) = A( I, J )
         } // 110
      } // 120
      RETURN

      // End of ZLAGSY

      }
