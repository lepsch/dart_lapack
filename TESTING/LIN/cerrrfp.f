      SUBROUTINE CERRRFP( NUNIT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NUNIT;
      // ..

*  =====================================================================

      // ..
      // .. Local Scalars ..
      int                INFO;
      COMPLEX            ALPHACMPLX
      REAL               ALPHA, BETA
      // ..
      // .. Local Arrays ..
      COMPLEX            A( 1, 1), B( 1, 1)
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, CTFSM, CTFTRI, CHFRK, CTFTTP, CTFTTR, CPFTRI, CPFTRF, CPFTRS, CTPTTF, CTPTTR, CTRTTF, CTRTTP
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      OK = true;
      A( 1, 1 ) = CMPLX( 1.0 , 1.0 )
      B( 1, 1 ) = CMPLX( 1.0 , 1.0  )
      ALPHACMPLX = CMPLX( 1.0 , 1.0  )
      ALPHA = 1.0
      BETA = 1.0

      SRNAMT = 'CPFTRF'
      INFOT = 1
      cpftrf('/', 'U', 0, A, INFO );
      chkxer('CPFTRF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      cpftrf('N', '/', 0, A, INFO );
      chkxer('CPFTRF', INFOT, NOUT, LERR, OK );
      INFOT = 3
      cpftrf('N', 'U', -1, A, INFO );
      chkxer('CPFTRF', INFOT, NOUT, LERR, OK );

      SRNAMT = 'CPFTRS'
      INFOT = 1
      cpftrs('/', 'U', 0, 0, A, B, 1, INFO );
      chkxer('CPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 2
      cpftrs('N', '/', 0, 0, A, B, 1, INFO );
      chkxer('CPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 3
      cpftrs('N', 'U', -1, 0, A, B, 1, INFO );
      chkxer('CPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 4
      cpftrs('N', 'U', 0, -1, A, B, 1, INFO );
      chkxer('CPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 7
      cpftrs('N', 'U', 0, 0, A, B, 0, INFO );
      chkxer('CPFTRS', INFOT, NOUT, LERR, OK );

      SRNAMT = 'CPFTRI'
      INFOT = 1
      cpftri('/', 'U', 0, A, INFO );
      chkxer('CPFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 2
      cpftri('N', '/', 0, A, INFO );
      chkxer('CPFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 3
      cpftri('N', 'U', -1, A, INFO );
      chkxer('CPFTRI', INFOT, NOUT, LERR, OK );

      SRNAMT = 'CTFSM '
      INFOT = 1
      ctfsm('/', 'L', 'U', 'C', 'U', 0, 0, ALPHACMPLX, A, B, 1 );
      chkxer('CTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ctfsm('N', '/', 'U', 'C', 'U', 0, 0, ALPHACMPLX, A, B, 1 );
      chkxer('CTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      ctfsm('N', 'L', '/', 'C', 'U', 0, 0, ALPHACMPLX, A, B, 1 );
      chkxer('CTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 4
      ctfsm('N', 'L', 'U', '/', 'U', 0, 0, ALPHACMPLX, A, B, 1 );
      chkxer('CTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      ctfsm('N', 'L', 'U', 'C', '/', 0, 0, ALPHACMPLX, A, B, 1 );
      chkxer('CTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 6
      ctfsm('N', 'L', 'U', 'C', 'U', -1, 0, ALPHACMPLX, A, B, 1 );
      chkxer('CTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      ctfsm('N', 'L', 'U', 'C', 'U', 0, -1, ALPHACMPLX, A, B, 1 );
      chkxer('CTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 11
      ctfsm('N', 'L', 'U', 'C', 'U', 0, 0, ALPHACMPLX, A, B, 0 );
      chkxer('CTFSM ', INFOT, NOUT, LERR, OK );

      SRNAMT = 'CTFTRI'
      INFOT = 1
      ctftri('/', 'L', 'N', 0, A, INFO );
      chkxer('CTFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ctftri('N', '/', 'N', 0, A, INFO );
      chkxer('CTFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 3
      ctftri('N', 'L', '/', 0, A, INFO );
      chkxer('CTFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 4
      ctftri('N', 'L', 'N', -1, A, INFO );
      chkxer('CTFTRI', INFOT, NOUT, LERR, OK );

      SRNAMT = 'CTFTTR'
      INFOT = 1
      ctfttr('/', 'U', 0, A, B, 1, INFO );
      chkxer('CTFTTR', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ctfttr('N', '/', 0, A, B, 1, INFO );
      chkxer('CTFTTR', INFOT, NOUT, LERR, OK );
      INFOT = 3
      ctfttr('N', 'U', -1, A, B, 1, INFO );
      chkxer('CTFTTR', INFOT, NOUT, LERR, OK );
      INFOT = 6
      ctfttr('N', 'U', 0, A, B, 0, INFO );
      chkxer('CTFTTR', INFOT, NOUT, LERR, OK );

      SRNAMT = 'CTRTTF'
      INFOT = 1
      ctrttf('/', 'U', 0, A, 1, B, INFO );
      chkxer('CTRTTF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ctrttf('N', '/', 0, A, 1, B, INFO );
      chkxer('CTRTTF', INFOT, NOUT, LERR, OK );
      INFOT = 3
      ctrttf('N', 'U', -1, A, 1, B, INFO );
      chkxer('CTRTTF', INFOT, NOUT, LERR, OK );
      INFOT = 5
      ctrttf('N', 'U', 0, A, 0, B, INFO );
      chkxer('CTRTTF', INFOT, NOUT, LERR, OK );

      SRNAMT = 'CTFTTP'
      INFOT = 1
      ctfttp('/', 'U', 0, A, B, INFO );
      chkxer('CTFTTP', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ctfttp('N', '/', 0, A, B, INFO );
      chkxer('CTFTTP', INFOT, NOUT, LERR, OK );
      INFOT = 3
      ctfttp('N', 'U', -1, A, B, INFO );
      chkxer('CTFTTP', INFOT, NOUT, LERR, OK );

      SRNAMT = 'CTPTTF'
      INFOT = 1
      ctpttf('/', 'U', 0, A, B, INFO );
      chkxer('CTPTTF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ctpttf('N', '/', 0, A, B, INFO );
      chkxer('CTPTTF', INFOT, NOUT, LERR, OK );
      INFOT = 3
      ctpttf('N', 'U', -1, A, B, INFO );
      chkxer('CTPTTF', INFOT, NOUT, LERR, OK );

      SRNAMT = 'CTRTTP'
      INFOT = 1
      ctrttp('/', 0, A, 1,  B, INFO );
      chkxer('CTRTTP', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ctrttp('U', -1, A, 1,  B, INFO );
      chkxer('CTRTTP', INFOT, NOUT, LERR, OK );
      INFOT = 4
      ctrttp('U', 0, A, 0,  B, INFO );
      chkxer('CTRTTP', INFOT, NOUT, LERR, OK );

      SRNAMT = 'CTPTTR'
      INFOT = 1
      ctpttr('/', 0, A, B, 1,  INFO );
      chkxer('CTPTTR', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ctpttr('U', -1, A, B, 1,  INFO );
      chkxer('CTPTTR', INFOT, NOUT, LERR, OK );
      INFOT = 5
      ctpttr('U', 0, A, B, 0, INFO );
      chkxer('CTPTTR', INFOT, NOUT, LERR, OK );

      SRNAMT = 'CHFRK '
      INFOT = 1
      chfrk('/', 'U', 'N', 0, 0, ALPHA, A, 1, BETA, B );
      chkxer('CHFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      chfrk('N', '/', 'N', 0, 0, ALPHA, A, 1, BETA, B );
      chkxer('CHFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      chfrk('N', 'U', '/', 0, 0, ALPHA, A, 1, BETA, B );
      chkxer('CHFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 4
      chfrk('N', 'U', 'N', -1, 0, ALPHA, A, 1, BETA, B );
      chkxer('CHFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      chfrk('N', 'U', 'N', 0, -1, ALPHA, A, 1, BETA, B );
      chkxer('CHFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 8
      chfrk('N', 'U', 'N', 0, 0, ALPHA, A, 0, BETA, B );
      chkxer('CHFRK ', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )
      } else {
         WRITE( NOUT, FMT = 9998 )
      }

 9999 FORMAT( 1X, 'COMPLEX RFP routines passed the tests of the ', 'error exits' )
 9998 FORMAT( ' *** RFP routines failed the tests of the error ', 'exits ***' )
      RETURN

      // End of CERRRFP

      }
