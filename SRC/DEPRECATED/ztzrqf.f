      SUBROUTINE ZTZRQF( M, N, A, LDA, TAU, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * )
      // ..

* =====================================================================

      // .. Parameters ..
      COMPLEX*16         CONE, CZERO
      const              CONE = ( 1.0D+0, 0.0D+0 ), CZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, K, M1;
      COMPLEX*16         ALPHA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX, MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZCOPY, ZGEMV, ZGERC, ZLACGV, ZLARFG
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.M ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('ZTZRQF', -INFO );
         RETURN
      }

      // Perform the factorization.

      IF( M.EQ.0 ) RETURN
      if ( M.EQ.N ) {
         for (I = 1; I <= N; I++) { // 10
            TAU( I ) = CZERO
   10    CONTINUE
      } else {
         M1 = MIN( M+1, N )
         DO 20 K = M, 1, -1

            // Use a Householder reflection to zero the kth row of A.
            // First set up the reflection.

            A( K, K ) = DCONJG( A( K, K ) )
            zlacgv(N-M, A( K, M1 ), LDA );
            ALPHA = A( K, K )
            zlarfg(N-M+1, ALPHA, A( K, M1 ), LDA, TAU( K ) );
            A( K, K ) = ALPHA
            TAU( K ) = DCONJG( TAU( K ) )

            if ( TAU( K ).NE.CZERO .AND. K.GT.1 ) {

               // We now perform the operation  A := A*P( k )**H.

               // Use the first ( k - 1 ) elements of TAU to store  a( k ),
               // where  a( k ) consists of the first ( k - 1 ) elements of
               // the  kth column  of  A.  Also  let  B  denote  the  first
               // ( k - 1 ) rows of the last ( n - m ) columns of A.

               zcopy(K-1, A( 1, K ), 1, TAU, 1 );

               // Form   w = a( k ) + B*z( k )  in TAU.

               zgemv('No transpose', K-1, N-M, CONE, A( 1, M1 ), LDA, A( K, M1 ), LDA, CONE, TAU, 1 );

               // Now form  a( k ) := a( k ) - conjg(tau)*w
               // and       B      := B      - conjg(tau)*w*z( k )**H.

               zaxpy(K-1, -DCONJG( TAU( K ) ), TAU, 1, A( 1, K ), 1 )                CALL ZGERC( K-1, N-M, -DCONJG( TAU( K ) ), TAU, 1, A( K, M1 ), LDA, A( 1, M1 ), LDA );
            }
   20    CONTINUE
      }

      RETURN

      // End of ZTZRQF

      }
