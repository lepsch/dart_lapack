      SUBROUTINE CTPLQT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int            INFO, LDA, LDB, LDT, N, M, L;
*     ..
*     .. Array Arguments ..
      COMPLEX     A( LDA, * ), B( LDB, * ), T( LDT, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX  ONE, ZERO
      PARAMETER( ZERO = ( 0.0E+0, 0.0E+0 ),ONE  = ( 1.0E+0, 0.0E+0 ) )
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
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -7
      ELSE IF( LDT.LT.MAX( 1, M ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTPLQT2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 ) RETURN
*
      DO I = 1, M
*
*        Generate elementary reflector H(I) to annihilate B(I,:)
*
         P = N-L+MIN( L, I )
         CALL CLARFG( P+1, A( I, I ), B( I, 1 ), LDB, T( 1, I ) )
         T(1,I)=CONJG(T(1,I))
         IF( I.LT.M ) THEN
            DO J = 1, P
               B( I, J ) = CONJG(B(I,J))
            END DO
*
*           W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)]
*
            DO J = 1, M-I
               T( M, J ) = (A( I+J, I ))
            END DO
            CALL CGEMV( 'N', M-I, P, ONE, B( I+1, 1 ), LDB, B( I, 1 ), LDB, ONE, T( M, 1 ), LDT )
*
*           C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H
*
            ALPHA = -(T( 1, I ))
            DO J = 1, M-I
               A( I+J, I ) = A( I+J, I ) + ALPHA*(T( M, J ))
            END DO
            CALL CGERC( M-I, P, (ALPHA),  T( M, 1 ), LDT, B( I, 1 ), LDB, B( I+1, 1 ), LDB )
            DO J = 1, P
               B( I, J ) = CONJG(B(I,J))
            END DO
         END IF
      END DO
*
      DO I = 2, M
*
*        T(I,1:I-1) := C(I:I-1,1:N)**H * (alpha * C(I,I:N))
*
         ALPHA = -(T( 1, I ))
         DO J = 1, I-1
            T( I, J ) = ZERO
         END DO
         P = MIN( I-1, L )
         NP = MIN( N-L+1, N )
         MP = MIN( P+1, M )
         DO J = 1, N-L+P
           B(I,J)=CONJG(B(I,J))
         END DO
*
*        Triangular part of B2
*
         DO J = 1, P
            T( I, J ) = (ALPHA*B( I, N-L+J ))
         END DO
         CALL CTRMV( 'L', 'N', 'N', P, B( 1, NP ), LDB, T( I, 1 ), LDT )
*
*        Rectangular part of B2
*
         CALL CGEMV( 'N', I-1-P, L,  ALPHA, B( MP, NP ), LDB, B( I, NP ), LDB, ZERO, T( I,MP ), LDT )
*
*        B1

*
         CALL CGEMV( 'N', I-1, N-L, ALPHA, B, LDB, B( I, 1 ), LDB, ONE, T( I, 1 ), LDT )
*

*
*        T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1)
*
         DO J = 1, I-1
            T(I,J)=CONJG(T(I,J))
         END DO
         CALL CTRMV( 'L', 'C', 'N', I-1, T, LDT, T( I, 1 ), LDT )
         DO J = 1, I-1
            T(I,J)=CONJG(T(I,J))
         END DO
         DO J = 1, N-L+P
            B(I,J)=CONJG(B(I,J))
         END DO
*
*        T(I,I) = tau(I)
*
         T( I, I ) = T( 1, I )
         T( 1, I ) = ZERO
      END DO
      DO I=1,M
         DO J= I+1,M
            T(I,J)=(T(J,I))
            T(J,I)=ZERO
         END DO
      END DO

*
*     End of CTPLQT2
*
      END
