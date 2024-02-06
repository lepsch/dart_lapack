      void chkxer(SRNAMT, INFOT, NOUT, LERR, OK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               LERR, OK;
      List<String>        srnamc.SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LEN_TRIM
      if ( !LERR ) {
         WRITE( NOUT, FMT = 9999 )INFOT,srnamc.SRNAMT( 1:LEN_TRIM(srnamc.SRNAMT ) );
         OK = false;
      }
      LERR = false;
      return;

 9999 FORMAT( ' *** Illegal value of parameter number ${.i2} not detected by ${} ***' );
      }
