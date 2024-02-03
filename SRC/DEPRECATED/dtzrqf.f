      void dtzrqf(M, N, A, LDA, TAU, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), TAU( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, K, M1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGEMV, DGER, DLARFG, XERBLA
      // ..
      // .. Executable Statements ..

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
         xerbla('DTZRQF', -INFO );
         return;
      }

      // Perform the factorization.

      if (M == 0) return;
      if ( M == N ) {
         for (I = 1; I <= N; I++) { // 10
            TAU( I ) = ZERO;
         } // 10
      } else {
         M1 = min( M+1, N );
         DO 20 K = M, 1, -1;

            // Use a Householder reflection to zero the kth row of A.
            // First set up the reflection.

            dlarfg(N-M+1, A( K, K ), A( K, M1 ), LDA, TAU( K ) );

            if ( ( TAU( K ) != ZERO ) && ( K > 1 ) ) {

               // We now perform the operation  A := A*P( k ).

               // Use the first ( k - 1 ) elements of TAU to store  a( k ),
               // where  a( k ) consists of the first ( k - 1 ) elements of
               // the  kth column  of  A.  Also  let  B  denote  the  first
               // ( k - 1 ) rows of the last ( n - m ) columns of A.

               dcopy(K-1, A( 1, K ), 1, TAU, 1 );

               // Form   w = a( k ) + B*z( k )  in TAU.

               dgemv('No transpose', K-1, N-M, ONE, A( 1, M1 ), LDA, A( K, M1 ), LDA, ONE, TAU, 1 );

               // Now form  a( k ) := a( k ) - tau*w
               // and       B      := B      - tau*w*z( k )**T.

               daxpy(K-1, -TAU( K ), TAU, 1, A( 1, K ), 1 );
               dger(K-1, N-M, -TAU( K ), TAU, 1, A( K, M1 ), LDA, A( 1, M1 ), LDA );
            }
         } // 20
      }

      return;
      }
