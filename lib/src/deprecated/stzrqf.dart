      void stzrqf(final int M, final int N, final Matrix<double> A, final int LDA, final int TAU, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      double               A( LDA, * ), TAU( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, K, M1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SGEMV, SGER, SLARFG, XERBLA

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
         xerbla('STZRQF', -INFO );
         return;
      }

      // Perform the factorization.

      if (M == 0) return;
      if ( M == N ) {
         for (I = 1; I <= N; I++) { // 10
            TAU[I] = ZERO;
         } // 10
      } else {
         M1 = min( M+1, N );
         for (K = M; K >= 1; K--) { // 20

            // Use a Householder reflection to zero the kth row of A.
            // First set up the reflection.

            slarfg(N-M+1, A( K, K ), A( K, M1 ), LDA, TAU( K ) );

            if ( ( TAU( K ) != ZERO ) && ( K > 1 ) ) {

               // We now perform the operation  A := A*P( k ).

               // Use the first ( k - 1 ) elements of TAU to store  a( k ),
               // where  a( k ) consists of the first ( k - 1 ) elements of
               // the  kth column  of  A.  Also  let  B  denote  the  first
               // ( k - 1 ) rows of the last ( n - m ) columns of A.

               scopy(K-1, A( 1, K ), 1, TAU, 1 );

               // Form   w = a( k ) + B*z( k )  in TAU.

               sgemv('No transpose', K-1, N-M, ONE, A( 1, M1 ), LDA, A( K, M1 ), LDA, ONE, TAU, 1 );

               // Now form  a( k ) := a( k ) - tau*w
               // and       B      := B      - tau*w*z( k )**T.

               saxpy(K-1, -TAU( K ), TAU, 1, A( 1, K ), 1 );
               sger(K-1, N-M, -TAU( K ), TAU, 1, A( K, M1 ), LDA, A( 1, M1 ), LDA );
            }
         } // 20
      }

      }
