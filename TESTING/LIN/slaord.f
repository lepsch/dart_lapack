      SUBROUTINE SLAORD( JOB, N, X, INCX )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOB;
      int                INCX, N;
      // ..
      // .. Array Arguments ..
      REAL               X( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, INC, IX, IXNEXT;
      REAL               TEMP
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      INC = ABS( INCX )
      IF( LSAME( JOB, 'I' ) ) THEN

         // Sort in increasing order

         DO 20 I = 2, N
            IX = 1 + ( I-1 )*INC
   10       CONTINUE
            IF( IX.EQ.1 ) GO TO 20
            IXNEXT = IX - INC
            IF( X( IX ).GT.X( IXNEXT ) ) THEN
               GO TO 20
            ELSE
               TEMP = X( IX )
               X( IX ) = X( IXNEXT )
               X( IXNEXT ) = TEMP
            END IF
            IX = IXNEXT
            GO TO 10
   20    CONTINUE

      ELSE IF( LSAME( JOB, 'D' ) ) THEN

         // Sort in decreasing order

         DO 40 I = 2, N
            IX = 1 + ( I-1 )*INC
   30       CONTINUE
            IF( IX.EQ.1 ) GO TO 40
            IXNEXT = IX - INC
            IF( X( IX ).LT.X( IXNEXT ) ) THEN
               GO TO 40
            ELSE
               TEMP = X( IX )
               X( IX ) = X( IXNEXT )
               X( IXNEXT ) = TEMP
            END IF
            IX = IXNEXT
            GO TO 30
   40    CONTINUE
      END IF
      RETURN

      // End of SLAORD

      END
