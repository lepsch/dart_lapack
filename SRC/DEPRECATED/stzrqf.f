      SUBROUTINE STZRQF( M, N, A, LDA, TAU, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, K, M1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SGEMV, SGER, SLARFG, XERBLA
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
         CALL XERBLA( 'STZRQF', -INFO )
         RETURN
      }

      // Perform the factorization.

      IF( M.EQ.0 ) RETURN
      if ( M.EQ.N ) {
         DO 10 I = 1, N
            TAU( I ) = ZERO
   10    CONTINUE
      } else {
         M1 = MIN( M+1, N )
         DO 20 K = M, 1, -1

            // Use a Householder reflection to zero the kth row of A.
            // First set up the reflection.

            CALL SLARFG( N-M+1, A( K, K ), A( K, M1 ), LDA, TAU( K ) )

            if ( ( TAU( K ).NE.ZERO ) .AND. ( K.GT.1 ) ) {

               // We now perform the operation  A := A*P( k ).

               // Use the first ( k - 1 ) elements of TAU to store  a( k ),
               // where  a( k ) consists of the first ( k - 1 ) elements of
               // the  kth column  of  A.  Also  let  B  denote  the  first
               // ( k - 1 ) rows of the last ( n - m ) columns of A.

               CALL SCOPY( K-1, A( 1, K ), 1, TAU, 1 )

               // Form   w = a( k ) + B*z( k )  in TAU.

               CALL SGEMV( 'No transpose', K-1, N-M, ONE, A( 1, M1 ), LDA, A( K, M1 ), LDA, ONE, TAU, 1 )

               // Now form  a( k ) := a( k ) - tau*w
               // and       B      := B      - tau*w*z( k )**T.

               CALL SAXPY( K-1, -TAU( K ), TAU, 1, A( 1, K ), 1 )
               CALL SGER( K-1, N-M, -TAU( K ), TAU, 1, A( K, M1 ), LDA, A( 1, M1 ), LDA )
            }
   20    CONTINUE
      }

      RETURN

      // End of STZRQF

      }
