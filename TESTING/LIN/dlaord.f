      SUBROUTINE DLAORD( JOB, N, X, INCX )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOB;
      int                INCX, N;
      // ..
      // .. Array Arguments ..
      double             X( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, INC, IX, IXNEXT;
      double             TEMP;
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
      if ( LSAME( JOB, 'I' ) ) {

         // Sort in increasing order

         DO 20 I = 2, N
            IX = 1 + ( I-1 )*INC
   10       CONTINUE
            IF( IX.EQ.1 ) GO TO 20
            IXNEXT = IX - INC
            if ( X( IX ).GT.X( IXNEXT ) ) {
               GO TO 20
            } else {
               TEMP = X( IX )
               X( IX ) = X( IXNEXT )
               X( IXNEXT ) = TEMP
            }
            IX = IXNEXT
            GO TO 10
   20    CONTINUE

      } else if ( LSAME( JOB, 'D' ) ) {

         // Sort in decreasing order

         DO 40 I = 2, N
            IX = 1 + ( I-1 )*INC
   30       CONTINUE
            IF( IX.EQ.1 ) GO TO 40
            IXNEXT = IX - INC
            if ( X( IX ).LT.X( IXNEXT ) ) {
               GO TO 40
            } else {
               TEMP = X( IX )
               X( IX ) = X( IXNEXT )
               X( IXNEXT ) = TEMP
            }
            IX = IXNEXT
            GO TO 30
   40    CONTINUE
      }
      RETURN

      // End of DLAORD

      }
