      void alaesm(PATH, OK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               OK;
      String             PATH;
      int                NOUT;
      // ..

// =====================================================================

      // .. Executable Statements ..

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH;
      } else {
         WRITE( NOUT, FMT = 9998 )PATH;
      }

 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits' );
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ', 'exits ***' );
      return;

      // End of ALAESM

      }
