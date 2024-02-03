      SUBROUTINE DLAGSY( N, K, D, A, LDA, ISEED, WORK, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, K, LDA, N;
*     ..
*     .. Array Arguments ..
      int                ISEED( 4 );
      double             A( LDA, * ), D( * ), WORK( * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE, HALF;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      int                I, J;
      double             ALPHA, TAU, WA, WB, WN;
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DGEMV, DGER, DLARNV, DSCAL, DSYMV, DSYR2, XERBLA
*     ..
*     .. External Functions ..
      double             DDOT, DNRM2;
      EXTERNAL           DDOT, DNRM2
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( K.LT.0 .OR. K.GT.N-1 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.LT.0 ) THEN
         CALL XERBLA( 'DLAGSY', -INFO )
         RETURN
      END IF
*
*     initialize lower triangle of A to diagonal matrix
*
      DO 20 J = 1, N
         DO 10 I = J + 1, N
            A( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
      DO 30 I = 1, N
         A( I, I ) = D( I )
   30 CONTINUE
*
*     Generate lower triangle of symmetric matrix
*
      DO 40 I = N - 1, 1, -1
*
*        generate random reflection
*
         CALL DLARNV( 3, ISEED, N-I+1, WORK )
         WN = DNRM2( N-I+1, WORK, 1 )
         WA = SIGN( WN, WORK( 1 ) )
         IF( WN.EQ.ZERO ) THEN
            TAU = ZERO
         ELSE
            WB = WORK( 1 ) + WA
            CALL DSCAL( N-I, ONE / WB, WORK( 2 ), 1 )
            WORK( 1 ) = ONE
            TAU = WB / WA
         END IF
*
*        apply random reflection to A(i:n,i:n) from the left
*        and the right
*
*        compute  y := tau * A * u
*
         CALL DSYMV( 'Lower', N-I+1, TAU, A( I, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 )
*
*        compute  v := y - 1/2 * tau * ( y, u ) * u
*
         ALPHA = -HALF*TAU*DDOT( N-I+1, WORK( N+1 ), 1, WORK, 1 )
         CALL DAXPY( N-I+1, ALPHA, WORK, 1, WORK( N+1 ), 1 )
*
*        apply the transformation as a rank-2 update to A(i:n,i:n)
*
         CALL DSYR2( 'Lower', N-I+1, -ONE, WORK, 1, WORK( N+1 ), 1, A( I, I ), LDA )
   40 CONTINUE
*
*     Reduce number of subdiagonals to K
*
      DO 60 I = 1, N - 1 - K
*
*        generate reflection to annihilate A(k+i+1:n,i)
*
         WN = DNRM2( N-K-I+1, A( K+I, I ), 1 )
         WA = SIGN( WN, A( K+I, I ) )
         IF( WN.EQ.ZERO ) THEN
            TAU = ZERO
         ELSE
            WB = A( K+I, I ) + WA
            CALL DSCAL( N-K-I, ONE / WB, A( K+I+1, I ), 1 )
            A( K+I, I ) = ONE
            TAU = WB / WA
         END IF
*
*        apply reflection to A(k+i:n,i+1:k+i-1) from the left
*
         CALL DGEMV( 'Transpose', N-K-I+1, K-1, ONE, A( K+I, I+1 ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 )          CALL DGER( N-K-I+1, K-1, -TAU, A( K+I, I ), 1, WORK, 1, A( K+I, I+1 ), LDA )
*
*        apply reflection to A(k+i:n,k+i:n) from the left and the right
*
*        compute  y := tau * A * u
*
         CALL DSYMV( 'Lower', N-K-I+1, TAU, A( K+I, K+I ), LDA, A( K+I, I ), 1, ZERO, WORK, 1 )
*
*        compute  v := y - 1/2 * tau * ( y, u ) * u
*
         ALPHA = -HALF*TAU*DDOT( N-K-I+1, WORK, 1, A( K+I, I ), 1 )
         CALL DAXPY( N-K-I+1, ALPHA, A( K+I, I ), 1, WORK, 1 )
*
*        apply symmetric rank-2 update to A(k+i:n,k+i:n)
*
         CALL DSYR2( 'Lower', N-K-I+1, -ONE, A( K+I, I ), 1, WORK, 1, A( K+I, K+I ), LDA )
*
         A( K+I, I ) = -WA
         DO 50 J = K + I + 1, N
            A( J, I ) = ZERO
   50    CONTINUE
   60 CONTINUE
*
*     Store full symmetric matrix
*
      DO 80 J = 1, N
         DO 70 I = J + 1, N
            A( J, I ) = A( I, J )
   70    CONTINUE
   80 CONTINUE
      RETURN
*
*     End of DLAGSY
*
      END
