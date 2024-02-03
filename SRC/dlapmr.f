      SUBROUTINE DLAPMR( FORWRD, M, N, X, LDX, K )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      bool               FORWRD;
      int                LDX, M, N;
      // ..
      // .. Array Arguments ..
      int                K( * );
      double             X( LDX, * );
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      int                I, IN, J, JJ;
      double             TEMP;
      // ..
      // .. Executable Statements ..
*
      IF( M.LE.1 ) RETURN
*
      DO 10 I = 1, M
         K( I ) = -K( I )
   10 CONTINUE
*
      IF( FORWRD ) THEN
*
         // Forward permutation
*
         DO 50 I = 1, M
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
            DO 30 JJ = 1, N
               TEMP = X( J, JJ )
               X( J, JJ ) = X( IN, JJ )
               X( IN, JJ ) = TEMP
   30       CONTINUE
*
            K( IN ) = -K( IN )
            J = IN
            IN = K( IN )
            GO TO 20
*
   40       CONTINUE
*
   50    CONTINUE
*
      ELSE
*
         // Backward permutation
*
         DO 90 I = 1, M
*
            IF( K( I ).GT.0 ) GO TO 80
*
            K( I ) = -K( I )
            J = K( I )
   60       CONTINUE
            IF( J.EQ.I ) GO TO 80
*
            DO 70 JJ = 1, N
               TEMP = X( I, JJ )
               X( I, JJ ) = X( J, JJ )
               X( J, JJ ) = TEMP
   70       CONTINUE
*
            K( J ) = -K( J )
            J = K( J )
            GO TO 60
*
   80       CONTINUE
*
   90    CONTINUE
*
      END IF
*
      RETURN
*
      // End of DLAPMR
*
      END
