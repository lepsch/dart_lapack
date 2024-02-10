      void ctzrqf(M, N, A, LDA, TAU, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      Complex            A( LDA, * ), TAU( * );
      // ..

      Complex            CONE, CZERO;
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      int                I, K, M1;
      Complex            ALPHA;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGEMV, CGERC, CLACGV, CLARFG, XERBLA

      // Test the input parameters.

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < M ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('CTZRQF', -INFO );
         return;
      }

      // Perform the factorization.

      if (M == 0) return;
      if ( M == N ) {
         for (I = 1; I <= N; I++) { // 10
            TAU[I] = CZERO;
         } // 10
      } else {
         M1 = min( M+1, N );
         for (K = M; K >= 1; K--) { // 20

            // Use a Householder reflection to zero the kth row of A.
            // First set up the reflection.

            A[K][K] = CONJG( A( K, K ) );
            clacgv(N-M, A( K, M1 ), LDA );
            ALPHA = A( K, K );
            clarfg(N-M+1, ALPHA, A( K, M1 ), LDA, TAU( K ) );
            A[K][K] = ALPHA;
            TAU[K] = CONJG( TAU( K ) );

            if ( TAU( K ) != CZERO && K > 1 ) {

               // We now perform the operation  A := A*P( k )**H.

               // Use the first ( k - 1 ) elements of TAU to store  a( k ),
               // where  a( k ) consists of the first ( k - 1 ) elements of
               // the  kth column  of  A.  Also  let  B  denote  the  first
               // ( k - 1 ) rows of the last ( n - m ) columns of A.

               ccopy(K-1, A( 1, K ), 1, TAU, 1 );

               // Form   w = a( k ) + B*z( k )  in TAU.

               cgemv('No transpose', K-1, N-M, CONE, A( 1, M1 ), LDA, A( K, M1 ), LDA, CONE, TAU, 1 );

               // Now form  a( k ) := a( k ) - conjg(tau)*w
               // and       B      := B      - conjg(tau)*w*z( k )**H.

               caxpy(K-1, -CONJG( TAU( K ) ), TAU, 1, A( 1, K ), 1 );
               cgerc(K-1, N-M, -CONJG( TAU( K ) ), TAU, 1, A( K, M1 ), LDA, A( 1, M1 ), LDA );
            }
         } // 20
      }

      }
