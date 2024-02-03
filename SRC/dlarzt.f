      SUBROUTINE DLARZT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DTRMV, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Check for currently supported options
*
      INFO = 0
      IF( .NOT.LSAME( DIRECT, 'B' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( STOREV, 'R' ) ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLARZT', -INFO )
         RETURN
      END IF
*
      DO 20 I = K, 1, -1
         IF( TAU( I ).EQ.ZERO ) THEN
*
*           H(i)  =  I
*
            DO 10 J = I, K
               T( J, I ) = ZERO
   10       CONTINUE
         ELSE
*
*           general case
*
            IF( I.LT.K ) THEN
*
*              T(i+1:k,i) = - tau(i) * V(i+1:k,1:n) * V(i,1:n)**T
*
               CALL DGEMV( 'No transpose', K-I, N, -TAU( I ), V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, T( I+1, I ), 1 )
*
*              T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)
*
               CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
            END IF
            T( I, I ) = TAU( I )
         END IF
   20 CONTINUE
      RETURN
*
*     End of DLARZT
*
      END
