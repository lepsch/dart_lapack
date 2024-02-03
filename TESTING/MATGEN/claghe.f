      SUBROUTINE CLAGHE( N, K, D, A, LDA, ISEED, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDA, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      REAL               D( * )
      COMPLEX            A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO, ONE, HALF
      const              ZERO = ( 0.0E+0, 0.0E+0 ), ONE = ( 1.0E+0, 0.0E+0 ), HALF = ( 0.5E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               WN
      COMPLEX            ALPHA, TAU, WA, WB
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CGEMV, CGERC, CHEMV, CHER2, CLARNV, CSCAL, XERBLA
      // ..
      // .. External Functions ..
      REAL               SCNRM2
      COMPLEX            CDOTC
      // EXTERNAL SCNRM2, CDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, REAL
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
         CALL XERBLA( 'CLAGHE', -INFO )
         RETURN
      }

      // initialize lower triangle of A to diagonal matrix

      DO 20 J = 1, N
         DO 10 I = J + 1, N
            A( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, N
         A( I, I ) = D( I )
   30 CONTINUE

      // Generate lower triangle of hermitian matrix

      DO 40 I = N - 1, 1, -1

         // generate random reflection

         CALL CLARNV( 3, ISEED, N-I+1, WORK )
         WN = SCNRM2( N-I+1, WORK, 1 )
         WA = ( WN / ABS( WORK( 1 ) ) )*WORK( 1 )
         if ( WN.EQ.ZERO ) {
            TAU = ZERO
         } else {
            WB = WORK( 1 ) + WA
            CALL CSCAL( N-I, ONE / WB, WORK( 2 ), 1 )
            WORK( 1 ) = ONE
            TAU = REAL( WB / WA )
         }

         // apply random reflection to A(i:n,i:n) from the left
         // and the right

         // compute  y := tau * A * u

         CALL CHEMV( 'Lower', N-I+1, TAU, A( I, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 )

         // compute  v := y - 1/2 * tau * ( y, u ) * u

         ALPHA = -HALF*TAU*CDOTC( N-I+1, WORK( N+1 ), 1, WORK, 1 )
         CALL CAXPY( N-I+1, ALPHA, WORK, 1, WORK( N+1 ), 1 )

         // apply the transformation as a rank-2 update to A(i:n,i:n)

         CALL CHER2( 'Lower', N-I+1, -ONE, WORK, 1, WORK( N+1 ), 1, A( I, I ), LDA )
   40 CONTINUE

      // Reduce number of subdiagonals to K

      DO 60 I = 1, N - 1 - K

         // generate reflection to annihilate A(k+i+1:n,i)

         WN = SCNRM2( N-K-I+1, A( K+I, I ), 1 )
         WA = ( WN / ABS( A( K+I, I ) ) )*A( K+I, I )
         if ( WN.EQ.ZERO ) {
            TAU = ZERO
         } else {
            WB = A( K+I, I ) + WA
            CALL CSCAL( N-K-I, ONE / WB, A( K+I+1, I ), 1 )
            A( K+I, I ) = ONE
            TAU = REAL( WB / WA )
         }

         // apply reflection to A(k+i:n,i+1:k+i-1) from the left

         CALL CGEMV( 'Conjugate transpose', N-K-I+1, K-1, ONE, A( K+I, I+1 ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 )          CALL CGERC( N-K-I+1, K-1, -TAU, A( K+I, I ), 1, WORK, 1, A( K+I, I+1 ), LDA )

         // apply reflection to A(k+i:n,k+i:n) from the left and the right

         // compute  y := tau * A * u

         CALL CHEMV( 'Lower', N-K-I+1, TAU, A( K+I, K+I ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 )

         // compute  v := y - 1/2 * tau * ( y, u ) * u

         ALPHA = -HALF*TAU*CDOTC( N-K-I+1, WORK, 1, A( K+I, I ), 1 )
         CALL CAXPY( N-K-I+1, ALPHA, A( K+I, I ), 1, WORK, 1 )

         // apply hermitian rank-2 update to A(k+i:n,k+i:n)

         CALL CHER2( 'Lower', N-K-I+1, -ONE, A( K+I, I ), 1, WORK, 1, A( K+I, K+I ), LDA )

         A( K+I, I ) = -WA
         DO 50 J = K + I + 1, N
            A( J, I ) = ZERO
   50    CONTINUE
   60 CONTINUE

      // Store full hermitian matrix

      DO 80 J = 1, N
         DO 70 I = J + 1, N
            A( J, I ) = CONJG( A( I, J ) )
   70    CONTINUE
   80 CONTINUE
      RETURN

      // End of CLAGHE

      }
