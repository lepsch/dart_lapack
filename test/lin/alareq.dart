      void alareq(final int PATH, final int NMATS, final Array<bool> DOTYPE, final int NTYPES, final int NIN, final int NOUT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NIN, NMATS, NOUT, NTYPES;
      bool               DOTYPE( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               FIRSTT;
      String             C1;
      String             INTSTR;
      String             LINE;
      int                I, I1, IC, J, K, LENP, NT;
      int                NREQ( 100 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LEN
      // ..
      // .. Data statements ..
      const INTSTR = '0123456789';

      if ( NMATS >= NTYPES ) {

         // Test everything if NMATS >= NTYPES.

         for (I = 1; I <= NTYPES; I++) { // 10
            DOTYPE[I] = true;
         } // 10
      } else {
         for (I = 1; I <= NTYPES; I++) { // 20
            DOTYPE[I] = false;
         } // 20
         FIRSTT = true;

         // Read a line of matrix types if 0 < NMATS < NTYPES.

         if ( NMATS > 0 ) {
            READ( NIN, FMT = '(A80)', END = 90 )LINE;
            LENP = LINE.length;
            I = 0;
            for (J = 1; J <= NMATS; J++) { // 60
               NREQ[J] = 0;
               I1 = 0;
               } // 30
               I = I + 1;
               if ( I > LENP ) {
                  if ( J == NMATS && I1 > 0 ) {
                     GO TO 60;
                  } else {
                     WRITE( NOUT, FMT = 9995 )LINE;
                     WRITE( NOUT, FMT = 9994 )NMATS;
                     GO TO 80;
                  }
               }
               if ( LINE( I: I ) != ' ' && LINE( I: I ) != ',' ) {
                  I1 = I;
                  C1 = LINE( I1: I1 );

               // Check that a valid integer was read

                  for (K = 1; K <= 10; K++) { // 40
                     if ( C1 == INTSTR( K: K ) ) {
                        IC = K - 1;
                        GO TO 50;
                     }
                  } // 40
                  WRITE( NOUT, FMT = 9996 )I, LINE;
                  WRITE( NOUT, FMT = 9994 )NMATS;
                  GO TO 80;
                  } // 50
                  NREQ[J] = 10*NREQ( J ) + IC;
                  GO TO 30;
               } else if ( I1 > 0 ) {
                  GO TO 60;
               } else {
                  GO TO 30;
               }
            } // 60
         }
         for (I = 1; I <= NMATS; I++) { // 70
            NT = NREQ( I );
            if ( NT > 0 && NT <= NTYPES ) {
               if ( DOTYPE( NT ) ) {
                  if (FIRSTT) WRITE( NOUT, FMT = * );
                  FIRSTT = false;
                  WRITE( NOUT, FMT = 9997 )NT, PATH;
               }
               DOTYPE[NT] = true;
            } else {
               WRITE( NOUT, FMT = 9999 )PATH, NT, NTYPES;
 9999          FORMAT( ' *** Invalid type request for ${.a3}, type  ${.i4}: must satisfy  1 <= type <= ${.i2}');
            }
         } // 70
         } // 80
      }
      return;

      } // 90
      WRITE( NOUT, FMT = 9998 )PATH;
 9998 FORMAT('\n *** End of file reached when trying to read matrix types for ${.a3}\n *** Check that you are requesting the right number of types for each path\n');
 9997 FORMAT( ' *** Warning:  duplicate request of matrix type ${.i2} for ${.a3}');
 9996 FORMAT( '\n\n *** Invalid integer value in column ', I2,; ' of input line:', /A79 )
 9995 FORMAT( '\n\n *** Not enough matrix types on input line', /A79 );
 9994 FORMAT( ' ==> Specify ${.i4} matrix types on this line or adjust NTYPES on previous line' );
      WRITE( NOUT, FMT = * );
      STOP;
      }
