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
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ), CZERO = ( 0.0D+0, 0.0D+0 ) )
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
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTZRQF', -INFO )
         RETURN
      END IF

      // Perform the factorization.

      IF( M.EQ.0 ) RETURN
      IF( M.EQ.N ) THEN
         DO 10 I = 1, N
            TAU( I ) = CZERO
   10    CONTINUE
      ELSE
         M1 = MIN( M+1, N )
         DO 20 K = M, 1, -1

            // Use a Householder reflection to zero the kth row of A.
            // First set up the reflection.

            A( K, K ) = DCONJG( A( K, K ) )
            CALL ZLACGV( N-M, A( K, M1 ), LDA )
            ALPHA = A( K, K )
            CALL ZLARFG( N-M+1, ALPHA, A( K, M1 ), LDA, TAU( K ) )
            A( K, K ) = ALPHA
            TAU( K ) = DCONJG( TAU( K ) )

            IF( TAU( K ).NE.CZERO .AND. K.GT.1 ) THEN

               // We now perform the operation  A := A*P( k )**H.

               // Use the first ( k - 1 ) elements of TAU to store  a( k ),
               // where  a( k ) consists of the first ( k - 1 ) elements of
              t // he  kth column  of  A.  Also  let  B  denote  the  first
               // ( k - 1 ) rows of the last ( n - m ) columns of A.

               CALL ZCOPY( K-1, A( 1, K ), 1, TAU, 1 )

               // Form   w = a( k ) + B*z( k )  in TAU.

               CALL ZGEMV( 'No transpose', K-1, N-M, CONE, A( 1, M1 ), LDA, A( K, M1 ), LDA, CONE, TAU, 1 )

               // Now form  a( k ) := a( k ) - conjg(tau)*w
               // and       B      := B      - conjg(tau)*w*z( k )**H.

               CALL ZAXPY( K-1, -DCONJG( TAU( K ) ), TAU, 1, A( 1, K ), 1 )                CALL ZGERC( K-1, N-M, -DCONJG( TAU( K ) ), TAU, 1, A( K, M1 ), LDA, A( 1, M1 ), LDA )
            END IF
   20    CONTINUE
      END IF

      RETURN

      // End of ZTZRQF

      }
