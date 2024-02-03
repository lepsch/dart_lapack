      SUBROUTINE ZGEQRT2( M, N, A, LDA, T, LDT, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int       INFO, LDA, LDT, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16   A( LDA, * ), T( LDT, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      COMPLEX*16  ONE, ZERO
      PARAMETER( ONE = (1.0D+00,0.0D+00), ZERO = (0.0D+00,0.0D+00) )
      // ..
      // .. Local Scalars ..
      int       I, K;
      COMPLEX*16   AII, ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLARFG, ZGEMV, ZGERC, ZTRMV, XERBLA
      // ..
      // .. Executable Statements ..
*
      // Test the input arguments
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( M.LT.N ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEQRT2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO I = 1, K
*
         // Generate elem. refl. H(i) to annihilate A(i+1:m,i), tau(I) -> T(I,1)
*
         CALL ZLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, T( I, 1 ) )
         IF( I.LT.N ) THEN
*
            // Apply H(i) to A(I:M,I+1:N) from the left
*
            AII = A( I, I )
            A( I, I ) = ONE
*
            // W(1:N-I) := A(I:M,I+1:N)^H * A(I:M,I) [W = T(:,N)]
*
            CALL ZGEMV( 'C',M-I+1, N-I, ONE, A( I, I+1 ), LDA, A( I, I ), 1, ZERO, T( 1, N ), 1 )
*
            // A(I:M,I+1:N) = A(I:m,I+1:N) + alpha*A(I:M,I)*W(1:N-1)^H
*
            ALPHA = -CONJG(T( I, 1 ))
            CALL ZGERC( M-I+1, N-I, ALPHA, A( I, I ), 1, T( 1, N ), 1, A( I, I+1 ), LDA )
            A( I, I ) = AII
         END IF
      END DO
*
      DO I = 2, N
         AII = A( I, I )
         A( I, I ) = ONE
*
         // T(1:I-1,I) := alpha * A(I:M,1:I-1)**H * A(I:M,I)
*
         ALPHA = -T( I, 1 )
         CALL ZGEMV( 'C', M-I+1, I-1, ALPHA, A( I, 1 ), LDA, A( I, I ), 1, ZERO, T( 1, I ), 1 )
         A( I, I ) = AII
*
         // T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)
*
         CALL ZTRMV( 'U', 'N', 'N', I-1, T, LDT, T( 1, I ), 1 )
*
            // T(I,I) = tau(I)
*
            T( I, I ) = T( I, 1 )
            T( I, 1) = ZERO
      END DO

*
      // End of ZGEQRT2
*
      END
