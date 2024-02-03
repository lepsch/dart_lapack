      void serrrfp(NUNIT ) {

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
      REAL               ALPHA, BETA;
      // ..
      // .. Local Arrays ..
      REAL               A( 1, 1), B( 1, 1);
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, STFSM, STFTRI, SSFRK, STFTTP, STFTTR, SPFTRI, SPFTRF, SPFTRS, STPTTF, STPTTR, STRTTF, STRTTP
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

      SRNAMT = 'SPFTRF';
      INFOT = 1;
      spftrf('/', 'U', 0, A, INFO );
      chkxer('SPFTRF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      spftrf('N', '/', 0, A, INFO );
      chkxer('SPFTRF', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      spftrf('N', 'U', -1, A, INFO );
      chkxer('SPFTRF', INFOT, NOUT, LERR, OK );

      SRNAMT = 'SPFTRS';
      INFOT = 1;
      spftrs('/', 'U', 0, 0, A, B, 1, INFO );
      chkxer('SPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      spftrs('N', '/', 0, 0, A, B, 1, INFO );
      chkxer('SPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      spftrs('N', 'U', -1, 0, A, B, 1, INFO );
      chkxer('SPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      spftrs('N', 'U', 0, -1, A, B, 1, INFO );
      chkxer('SPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      spftrs('N', 'U', 0, 0, A, B, 0, INFO );
      chkxer('SPFTRS', INFOT, NOUT, LERR, OK );

      SRNAMT = 'SPFTRI';
      INFOT = 1;
      spftri('/', 'U', 0, A, INFO );
      chkxer('SPFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      spftri('N', '/', 0, A, INFO );
      chkxer('SPFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      spftri('N', 'U', -1, A, INFO );
      chkxer('SPFTRI', INFOT, NOUT, LERR, OK );

      SRNAMT = 'STFSM ';
      INFOT = 1;
      stfsm('/', 'L', 'U', 'T', 'U', 0, 0, ALPHA, A, B, 1 );
      chkxer('STFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      stfsm('N', '/', 'U', 'T', 'U', 0, 0, ALPHA, A, B, 1 );
      chkxer('STFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stfsm('N', 'L', '/', 'T', 'U', 0, 0, ALPHA, A, B, 1 );
      chkxer('STFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      stfsm('N', 'L', 'U', '/', 'U', 0, 0, ALPHA, A, B, 1 );
      chkxer('STFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      stfsm('N', 'L', 'U', 'T', '/', 0, 0, ALPHA, A, B, 1 );
      chkxer('STFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      stfsm('N', 'L', 'U', 'T', 'U', -1, 0, ALPHA, A, B, 1 );
      chkxer('STFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      stfsm('N', 'L', 'U', 'T', 'U', 0, -1, ALPHA, A, B, 1 );
      chkxer('STFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      stfsm('N', 'L', 'U', 'T', 'U', 0, 0, ALPHA, A, B, 0 );
      chkxer('STFSM ', INFOT, NOUT, LERR, OK );

      SRNAMT = 'STFTRI';
      INFOT = 1;
      stftri('/', 'L', 'N', 0, A, INFO );
      chkxer('STFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      stftri('N', '/', 'N', 0, A, INFO );
      chkxer('STFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stftri('N', 'L', '/', 0, A, INFO );
      chkxer('STFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      stftri('N', 'L', 'N', -1, A, INFO );
      chkxer('STFTRI', INFOT, NOUT, LERR, OK );

      SRNAMT = 'STFTTR';
      INFOT = 1;
      stfttr('/', 'U', 0, A, B, 1, INFO );
      chkxer('STFTTR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      stfttr('N', '/', 0, A, B, 1, INFO );
      chkxer('STFTTR', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stfttr('N', 'U', -1, A, B, 1, INFO );
      chkxer('STFTTR', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      stfttr('N', 'U', 0, A, B, 0, INFO );
      chkxer('STFTTR', INFOT, NOUT, LERR, OK );

      SRNAMT = 'STRTTF';
      INFOT = 1;
      strttf('/', 'U', 0, A, 1, B, INFO );
      chkxer('STRTTF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      strttf('N', '/', 0, A, 1, B, INFO );
      chkxer('STRTTF', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      strttf('N', 'U', -1, A, 1, B, INFO );
      chkxer('STRTTF', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      strttf('N', 'U', 0, A, 0, B, INFO );
      chkxer('STRTTF', INFOT, NOUT, LERR, OK );

      SRNAMT = 'STFTTP';
      INFOT = 1;
      stfttp('/', 'U', 0, A, B, INFO );
      chkxer('STFTTP', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      stfttp('N', '/', 0, A, B, INFO );
      chkxer('STFTTP', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stfttp('N', 'U', -1, A, B, INFO );
      chkxer('STFTTP', INFOT, NOUT, LERR, OK );

      SRNAMT = 'STPTTF';
      INFOT = 1;
      stpttf('/', 'U', 0, A, B, INFO );
      chkxer('STPTTF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      stpttf('N', '/', 0, A, B, INFO );
      chkxer('STPTTF', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stpttf('N', 'U', -1, A, B, INFO );
      chkxer('STPTTF', INFOT, NOUT, LERR, OK );

      SRNAMT = 'STRTTP';
      INFOT = 1;
      strttp('/', 0, A, 1,  B, INFO );
      chkxer('STRTTP', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      strttp('U', -1, A, 1,  B, INFO );
      chkxer('STRTTP', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      strttp('U', 0, A, 0,  B, INFO );
      chkxer('STRTTP', INFOT, NOUT, LERR, OK );

      SRNAMT = 'STPTTR';
      INFOT = 1;
      stpttr('/', 0, A, B, 1,  INFO );
      chkxer('STPTTR', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      stpttr('U', -1, A, B, 1,  INFO );
      chkxer('STPTTR', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      stpttr('U', 0, A, B, 0, INFO );
      chkxer('STPTTR', INFOT, NOUT, LERR, OK );

      SRNAMT = 'SSFRK ';
      INFOT = 1;
      ssfrk('/', 'U', 'N', 0, 0, ALPHA, A, 1, BETA, B );
      chkxer('SSFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      ssfrk('N', '/', 'N', 0, 0, ALPHA, A, 1, BETA, B );
      chkxer('SSFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      ssfrk('N', 'U', '/', 0, 0, ALPHA, A, 1, BETA, B );
      chkxer('SSFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      ssfrk('N', 'U', 'N', -1, 0, ALPHA, A, 1, BETA, B );
      chkxer('SSFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      ssfrk('N', 'U', 'N', 0, -1, ALPHA, A, 1, BETA, B );
      chkxer('SSFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      ssfrk('N', 'U', 'N', 0, 0, ALPHA, A, 0, BETA, B );
      chkxer('SSFRK ', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 );
      } else {
         WRITE( NOUT, FMT = 9998 );
      }

 9999 FORMAT( 1X, 'REAL RFP routines passed the tests of ', 'the error exits' );
 9998 FORMAT( ' *** RFP routines failed the tests of the error ', 'exits ***' );
      return;
      }
