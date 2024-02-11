      void alaesm(final int PATH, final int OK, final int NOUT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
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

 9999 FORMAT(' ${.a3} routines passed the tests of the error exits' );
 9998 FORMAT( ' *** ${.a3} routines failed the tests of the error exits ***' );
      }
