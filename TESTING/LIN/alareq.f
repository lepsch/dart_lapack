      SUBROUTINE ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NIN, NMATS, NOUT, NTYPES;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               FIRSTT;
      String             C1;
      String             INTSTR;
      String             LINE;
      int                I, I1, IC, J, K, LENP, NT;
      // ..
      // .. Local Arrays ..
      int                NREQ( 100 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LEN
      // ..
      // .. Data statements ..
      DATA               INTSTR / '0123456789' /
      // ..
      // .. Executable Statements ..

      if ( NMATS.GE.NTYPES ) {

         // Test everything if NMATS >= NTYPES.

         DO 10 I = 1, NTYPES
            DOTYPE( I ) = .TRUE.
   10    CONTINUE
      } else {
         DO 20 I = 1, NTYPES
            DOTYPE( I ) = .FALSE.
   20    CONTINUE
         FIRSTT = .TRUE.

         // Read a line of matrix types if 0 < NMATS < NTYPES.

         if ( NMATS.GT.0 ) {
            READ( NIN, FMT = '(A80)', END = 90 )LINE
            LENP = LEN( LINE )
            I = 0
            DO 60 J = 1, NMATS
               NREQ( J ) = 0
               I1 = 0
   30          CONTINUE
               I = I + 1
               if ( I.GT.LENP ) {
                  if ( J.EQ.NMATS .AND. I1.GT.0 ) {
                     GO TO 60
                  } else {
                     WRITE( NOUT, FMT = 9995 )LINE
                     WRITE( NOUT, FMT = 9994 )NMATS
                     GO TO 80
                  }
               }
               if ( LINE( I: I ).NE.' ' .AND. LINE( I: I ).NE.',' ) {
                  I1 = I
                  C1 = LINE( I1: I1 )

               // Check that a valid integer was read

                  DO 40 K = 1, 10
                     if ( C1.EQ.INTSTR( K: K ) ) {
                        IC = K - 1
                        GO TO 50
                     }
   40             CONTINUE
                  WRITE( NOUT, FMT = 9996 )I, LINE
                  WRITE( NOUT, FMT = 9994 )NMATS
                  GO TO 80
   50             CONTINUE
                  NREQ( J ) = 10*NREQ( J ) + IC
                  GO TO 30
               } else if ( I1.GT.0 ) {
                  GO TO 60
               } else {
                  GO TO 30
               }
   60       CONTINUE
         }
         DO 70 I = 1, NMATS
            NT = NREQ( I )
            if ( NT.GT.0 .AND. NT.LE.NTYPES ) {
               if ( DOTYPE( NT ) ) {
                  IF( FIRSTT ) WRITE( NOUT, FMT = * )
                  FIRSTT = .FALSE.
                  WRITE( NOUT, FMT = 9997 )NT, PATH
               }
               DOTYPE( NT ) = .TRUE.
            } else {
               WRITE( NOUT, FMT = 9999 )PATH, NT, NTYPES
 9999          FORMAT( ' *** Invalid type request for ', A3, ', type  ', I4, ': must satisfy  1 <= type <= ', I2 )
            }
   70    CONTINUE
   80    CONTINUE
      }
      RETURN

   90 CONTINUE
      WRITE( NOUT, FMT = 9998 )PATH
 9998 FORMAT( /' *** End of file reached when trying to read matrix ', 'types for ', A3, /' *** Check that you are requesting the', ' right number of types for each path', / )
 9997 FORMAT( ' *** Warning:  duplicate request of matrix type ', I2, ' for ', A3 )
 9996 FORMAT( //' *** Invalid int     value in column ', I2,; ' of input', ' line:', /A79 )
 9995 FORMAT( //' *** Not enough matrix types on input line', /A79 )
 9994 FORMAT( ' ==> Specify ', I4, ' matrix types on this line or ', 'adjust NTYPES on previous line' )
      WRITE( NOUT, FMT = * )
      STOP

      // End of ALAREQ

      }
