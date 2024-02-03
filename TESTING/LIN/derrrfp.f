      SUBROUTINE DERRRFP( NUNIT );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NUNIT;
      // ..

// =====================================================================

      // ..
      // .. Local Scalars ..
      int                INFO;
      double             ALPHA, BETA;
      // ..
      // .. Local Arrays ..
      double             A( 1, 1), B( 1, 1);
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, DTFSM, DTFTRI, DSFRK, DTFTTP, DTFTTR, DPFTRI, DPFTRF, DPFTRS, DTPTTF, DTPTTR, DTRTTF, DTRTTP
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      OK = true;
      A( 1, 1 ) = 1.0;
      B( 1, 1 ) = 1.0;
      ALPHA     = 1.0;
      BETA      = 1.0;

      SRNAMT = 'DPFTRF';
      INFOT = 1;
      dpftrf('/', 'U', 0, A, INFO );
      chkxer('DPFTRF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dpftrf('N', '/', 0, A, INFO );
      chkxer('DPFTRF', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dpftrf('N', 'U', -1, A, INFO );
      chkxer('DPFTRF', INFOT, NOUT, LERR, OK );

      SRNAMT = 'DPFTRS';
      INFOT = 1;
      dpftrs('/', 'U', 0, 0, A, B, 1, INFO );
      chkxer('DPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dpftrs('N', '/', 0, 0, A, B, 1, INFO );
      chkxer('DPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dpftrs('N', 'U', -1, 0, A, B, 1, INFO );
      chkxer('DPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dpftrs('N', 'U', 0, -1, A, B, 1, INFO );
      chkxer('DPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      dpftrs('N', 'U', 0, 0, A, B, 0, INFO );
      chkxer('DPFTRS', INFOT, NOUT, LERR, OK );

      SRNAMT = 'DPFTRI';
      INFOT = 1;
      dpftri('/', 'U', 0, A, INFO );
      chkxer('DPFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dpftri('N', '/', 0, A, INFO );
      chkxer('DPFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dpftri('N', 'U', -1, A, INFO );
      chkxer('DPFTRI', INFOT, NOUT, LERR, OK );

      SRNAMT = 'DTFSM ';
      INFOT = 1;
      dtfsm('/', 'L', 'U', 'T', 'U', 0, 0, ALPHA, A, B, 1 );
      chkxer('DTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtfsm('N', '/', 'U', 'T', 'U', 0, 0, ALPHA, A, B, 1 );
      chkxer('DTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtfsm('N', 'L', '/', 'T', 'U', 0, 0, ALPHA, A, B, 1 );
      chkxer('DTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dtfsm('N', 'L', 'U', '/', 'U', 0, 0, ALPHA, A, B, 1 );
      chkxer('DTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      dtfsm('N', 'L', 'U', 'T', '/', 0, 0, ALPHA, A, B, 1 );
      chkxer('DTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      dtfsm('N', 'L', 'U', 'T', 'U', -1, 0, ALPHA, A, B, 1 );
      chkxer('DTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      dtfsm('N', 'L', 'U', 'T', 'U', 0, -1, ALPHA, A, B, 1 );
      chkxer('DTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      dtfsm('N', 'L', 'U', 'T', 'U', 0, 0, ALPHA, A, B, 0 );
      chkxer('DTFSM ', INFOT, NOUT, LERR, OK );

      SRNAMT = 'DTFTRI';
      INFOT = 1;
      dtftri('/', 'L', 'N', 0, A, INFO );
      chkxer('DTFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtftri('N', '/', 'N', 0, A, INFO );
      chkxer('DTFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtftri('N', 'L', '/', 0, A, INFO );
      chkxer('DTFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dtftri('N', 'L', 'N', -1, A, INFO );
      chkxer('DTFTRI', INFOT, NOUT, LERR, OK );

      SRNAMT = 'DTFTTR';
      INFOT = 1;
      dtfttr('/', 'U', 0, A, B, 1, INFO );
      chkxer('DTFTTR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtfttr('N', '/', 0, A, B, 1, INFO );
      chkxer('DTFTTR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtfttr('N', 'U', -1, A, B, 1, INFO );
      chkxer('DTFTTR', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      dtfttr('N', 'U', 0, A, B, 0, INFO );
      chkxer('DTFTTR', INFOT, NOUT, LERR, OK );

      SRNAMT = 'DTRTTF';
      INFOT = 1;
      dtrttf('/', 'U', 0, A, 1, B, INFO );
      chkxer('DTRTTF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtrttf('N', '/', 0, A, 1, B, INFO );
      chkxer('DTRTTF', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtrttf('N', 'U', -1, A, 1, B, INFO );
      chkxer('DTRTTF', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      dtrttf('N', 'U', 0, A, 0, B, INFO );
      chkxer('DTRTTF', INFOT, NOUT, LERR, OK );

      SRNAMT = 'DTFTTP';
      INFOT = 1;
      dtfttp('/', 'U', 0, A, B, INFO );
      chkxer('DTFTTP', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtfttp('N', '/', 0, A, B, INFO );
      chkxer('DTFTTP', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtfttp('N', 'U', -1, A, B, INFO );
      chkxer('DTFTTP', INFOT, NOUT, LERR, OK );

      SRNAMT = 'DTPTTF';
      INFOT = 1;
      dtpttf('/', 'U', 0, A, B, INFO );
      chkxer('DTPTTF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtpttf('N', '/', 0, A, B, INFO );
      chkxer('DTPTTF', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtpttf('N', 'U', -1, A, B, INFO );
      chkxer('DTPTTF', INFOT, NOUT, LERR, OK );

      SRNAMT = 'DTRTTP';
      INFOT = 1;
      dtrttp('/', 0, A, 1,  B, INFO );
      chkxer('DTRTTP', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtrttp('U', -1, A, 1,  B, INFO );
      chkxer('DTRTTP', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dtrttp('U', 0, A, 0,  B, INFO );
      chkxer('DTRTTP', INFOT, NOUT, LERR, OK );

      SRNAMT = 'DTPTTR';
      INFOT = 1;
      dtpttr('/', 0, A, B, 1,  INFO );
      chkxer('DTPTTR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtpttr('U', -1, A, B, 1,  INFO );
      chkxer('DTPTTR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      dtpttr('U', 0, A, B, 0, INFO );
      chkxer('DTPTTR', INFOT, NOUT, LERR, OK );

      SRNAMT = 'DSFRK ';
      INFOT = 1;
      dsfrk('/', 'U', 'N', 0, 0, ALPHA, A, 1, BETA, B );
      chkxer('DSFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dsfrk('N', '/', 'N', 0, 0, ALPHA, A, 1, BETA, B );
      chkxer('DSFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dsfrk('N', 'U', '/', 0, 0, ALPHA, A, 1, BETA, B );
      chkxer('DSFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dsfrk('N', 'U', 'N', -1, 0, ALPHA, A, 1, BETA, B );
      chkxer('DSFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      dsfrk('N', 'U', 'N', 0, -1, ALPHA, A, 1, BETA, B );
      chkxer('DSFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      dsfrk('N', 'U', 'N', 0, 0, ALPHA, A, 0, BETA, B );
      chkxer('DSFRK ', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 );
      } else {
         WRITE( NOUT, FMT = 9998 );
      }

 9999 FORMAT( 1X, 'double           RFP routines passed the tests of ',; 'the error exits' )
 9998 FORMAT( ' *** RFP routines failed the tests of the error ', 'exits ***' );
      return;

      // End of DERRRFP

      }
