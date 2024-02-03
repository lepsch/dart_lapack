      SUBROUTINE DLAORD( JOB, N, X, INCX );

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

      INC = ABS( INCX );
      if ( LSAME( JOB, 'I' ) ) {

         // Sort in increasing order

         for (I = 2; I <= N; I++) { // 20
            IX = 1 + ( I-1 )*INC;
            } // 10
            if (IX == 1) GO TO 20;
            IXNEXT = IX - INC;
            if ( X( IX ) > X( IXNEXT ) ) {
               GO TO 20;
            } else {
               TEMP = X( IX );
               X( IX ) = X( IXNEXT );
               X( IXNEXT ) = TEMP;
            }
            IX = IXNEXT;
            GO TO 10;
         } // 20

      } else if ( LSAME( JOB, 'D' ) ) {

         // Sort in decreasing order

         for (I = 2; I <= N; I++) { // 40
            IX = 1 + ( I-1 )*INC;
            } // 30
            if (IX == 1) GO TO 40;
            IXNEXT = IX - INC;
            if ( X( IX ) < X( IXNEXT ) ) {
               GO TO 40;
            } else {
               TEMP = X( IX );
               X( IX ) = X( IXNEXT );
               X( IXNEXT ) = TEMP;
            }
            IX = IXNEXT;
            GO TO 30;
         } // 40
      }
      RETURN;

      // End of DLAORD

      }
