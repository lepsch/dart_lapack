      SUBROUTINE CHKXER( SRNAMT, INFOT, NOUT, LERR, OK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LERR, OK;
      List<String>         SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LEN_TRIM
      // ..
      // .. Executable Statements ..
      if ( !LERR ) {
         WRITE( NOUT, FMT = 9999 )INFOT, SRNAMT( 1:LEN_TRIM( SRNAMT ) );
         OK = false;
      }
      LERR = false;
      return;

 9999 FORMAT( ' *** Illegal value of parameter number ', I2, ' not detected by ', A6, ' ***' );

      // End of CHKXER

      }
