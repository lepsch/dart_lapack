      SUBROUTINE CTPQRT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int       INFO, LDA, LDB, LDT, N, M, L;
*     ..
*     .. Array Arguments ..
      COMPLEX   A( LDA, * ), B( LDB, * ), T( LDT, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX  ONE, ZERO
      PARAMETER( ONE = (1.0,0.0), ZERO = (0.0,0.0) )
*     ..
*     .. Local Scalars ..
      int       I, J, P, MP, NP;
      COMPLEX   ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL  CLARFG, CGEMV, CGERC, CTRMV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( L.LT.0 .OR. L.GT.MIN(M,N) ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -7
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTPQRT2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 ) RETURN
*
      DO I = 1, N
*
*        Generate elementary reflector H(I) to annihilate B(:,I)
*
         P = M-L+MIN( L, I )
         CALL CLARFG( P+1, A( I, I ), B( 1, I ), 1, T( I, 1 ) )
         IF( I.LT.N ) THEN
*
*           W(1:N-I) := C(I:M,I+1:N)**H * C(I:M,I) [use W = T(:,N)]
*
            DO J = 1, N-I
               T( J, N ) = CONJG(A( I, I+J ))
            END DO
            CALL CGEMV( 'C', P, N-I, ONE, B( 1, I+1 ), LDB, B( 1, I ), 1, ONE, T( 1, N ), 1 )
*
*           C(I:M,I+1:N) = C(I:m,I+1:N) + alpha*C(I:M,I)*W(1:N-1)**H
*
            ALPHA = -CONJG(T( I, 1 ))
            DO J = 1, N-I
               A( I, I+J ) = A( I, I+J ) + ALPHA*CONJG(T( J, N ))
            END DO
            CALL CGERC( P, N-I, ALPHA, B( 1, I ), 1, T( 1, N ), 1, B( 1, I+1 ), LDB )
         END IF
      END DO
*
      DO I = 2, N
*
*        T(1:I-1,I) := C(I:M,1:I-1)**H * (alpha * C(I:M,I))
*
         ALPHA = -T( I, 1 )

         DO J = 1, I-1
            T( J, I ) = ZERO
         END DO
         P = MIN( I-1, L )
         MP = MIN( M-L+1, M )
         NP = MIN( P+1, N )
*
*        Triangular part of B2
*
         DO J = 1, P
            T( J, I ) = ALPHA*B( M-L+J, I )
         END DO
         CALL CTRMV( 'U', 'C', 'N', P, B( MP, 1 ), LDB, T( 1, I ), 1 )
*
*        Rectangular part of B2
*
         CALL CGEMV( 'C', L, I-1-P, ALPHA, B( MP, NP ), LDB, B( MP, I ), 1, ZERO, T( NP, I ), 1 )
*
*        B1
*
         CALL CGEMV( 'C', M-L, I-1, ALPHA, B, LDB, B( 1, I ), 1, ONE, T( 1, I ), 1 )
*
*        T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)
*
         CALL CTRMV( 'U', 'N', 'N', I-1, T, LDT, T( 1, I ), 1 )
*
*        T(I,I) = tau(I)
*
         T( I, I ) = T( I, 1 )
         T( I, 1 ) = ZERO
      END DO

*
*     End of CTPQRT2
*
      END
