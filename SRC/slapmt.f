      SUBROUTINE SLAPMT( FORWRD, M, N, X, LDX, K )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            FORWRD
      INTEGER            LDX, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            K( * )
      REAL               X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, II, J, IN
      REAL               TEMP
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.1 ) RETURN
*
      DO 10 I = 1, N
         K( I ) = -K( I )
   10 CONTINUE
*
      IF( FORWRD ) THEN
*
*        Forward permutation
*
         DO 60 I = 1, N
*
            IF( K( I ).GT.0 ) GO TO 40
*
            J = I
            K( J ) = -K( J )
            IN = K( J )
*
   20       CONTINUE
            IF( K( IN ).GT.0 ) GO TO 40
*
            DO 30 II = 1, M
               TEMP = X( II, J )
               X( II, J ) = X( II, IN )
               X( II, IN ) = TEMP
   30       CONTINUE
*
            K( IN ) = -K( IN )
            J = IN
            IN = K( IN )
            GO TO 20
*
   40       CONTINUE
*
   60    CONTINUE
*
      ELSE
*
*        Backward permutation
*
         DO 110 I = 1, N
*
            IF( K( I ).GT.0 ) GO TO 100
*
            K( I ) = -K( I )
            J = K( I )
   80       CONTINUE
            IF( J.EQ.I ) GO TO 100
*
            DO 90 II = 1, M
               TEMP = X( II, I )
               X( II, I ) = X( II, J )
               X( II, J ) = TEMP
   90       CONTINUE
*
            K( J ) = -K( J )
            J = K( J )
            GO TO 80
*
  100       CONTINUE

  110    CONTINUE
*
      END IF
*
      RETURN
*
*     End of SLAPMT
*
      END
