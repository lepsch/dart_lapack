      SUBROUTINE CTZRQF( M, N, A, LDA, TAU, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), TAU( * )
      // ..

* =====================================================================

      // .. Parameters ..
      COMPLEX            CONE, CZERO
      const              CONE = ( 1.0E+0, 0.0E+0 ), CZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, K, M1;
      COMPLEX            ALPHA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGEMV, CGERC, CLACGV, CLARFG, XERBLA
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
         xerbla('CTZRQF', -INFO );
         RETURN
      }

      // Perform the factorization.

      IF( M.EQ.0 ) RETURN
      if ( M.EQ.N ) {
         DO 10 I = 1, N
            TAU( I ) = CZERO
   10    CONTINUE
      } else {
         M1 = MIN( M+1, N )
         DO 20 K = M, 1, -1

            // Use a Householder reflection to zero the kth row of A.
            // First set up the reflection.

            A( K, K ) = CONJG( A( K, K ) )
            clacgv(N-M, A( K, M1 ), LDA );
            ALPHA = A( K, K )
            clarfg(N-M+1, ALPHA, A( K, M1 ), LDA, TAU( K ) );
            A( K, K ) = ALPHA
            TAU( K ) = CONJG( TAU( K ) )

            if ( TAU( K ).NE.CZERO .AND. K.GT.1 ) {

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

               caxpy(K-1, -CONJG( TAU( K ) ), TAU, 1, A( 1, K ), 1 )                CALL CGERC( K-1, N-M, -CONJG( TAU( K ) ), TAU, 1, A( K, M1 ), LDA, A( 1, M1 ), LDA );
            }
   20    CONTINUE
      }

      RETURN

      // End of CTZRQF

      }
