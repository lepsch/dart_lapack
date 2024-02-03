      SUBROUTINE ZERRRFP( NUNIT )

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
      double             ALPHA, BETA;
      COMPLEX*16         CALPHA
      // ..
      // .. Local Arrays ..
      COMPLEX*16         A( 1, 1), B( 1, 1)
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZTFSM, ZTFTRI, ZHFRK, ZTFTTP, ZTFTTR, ZPFTRI, ZPFTRF, ZPFTRS, ZTPTTF, ZTPTTR, ZTRTTF, ZTRTTP
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      OK = .TRUE.
      A( 1, 1 ) = DCMPLX( 1.0D0 , 1.0D0  )
      B( 1, 1 ) = DCMPLX( 1.0D0 , 1.0D0  )
      ALPHA     = 1.0D0
      CALPHA    = DCMPLX( 1.0D0 , 1.0D0  )
      BETA      = 1.0D0

      SRNAMT = 'ZPFTRF'
      INFOT = 1
      zpftrf('/', 'U', 0, A, INFO );
      chkxer('ZPFTRF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zpftrf('N', '/', 0, A, INFO );
      chkxer('ZPFTRF', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zpftrf('N', 'U', -1, A, INFO );
      chkxer('ZPFTRF', INFOT, NOUT, LERR, OK );

      SRNAMT = 'ZPFTRS'
      INFOT = 1
      zpftrs('/', 'U', 0, 0, A, B, 1, INFO );
      chkxer('ZPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zpftrs('N', '/', 0, 0, A, B, 1, INFO );
      chkxer('ZPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zpftrs('N', 'U', -1, 0, A, B, 1, INFO );
      chkxer('ZPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zpftrs('N', 'U', 0, -1, A, B, 1, INFO );
      chkxer('ZPFTRS', INFOT, NOUT, LERR, OK );
      INFOT = 7
      zpftrs('N', 'U', 0, 0, A, B, 0, INFO );
      chkxer('ZPFTRS', INFOT, NOUT, LERR, OK );

      SRNAMT = 'ZPFTRI'
      INFOT = 1
      zpftri('/', 'U', 0, A, INFO );
      chkxer('ZPFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zpftri('N', '/', 0, A, INFO );
      chkxer('ZPFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zpftri('N', 'U', -1, A, INFO );
      chkxer('ZPFTRI', INFOT, NOUT, LERR, OK );

      SRNAMT = 'ZTFSM '
      INFOT = 1
      ztfsm('/', 'L', 'U', 'C', 'U', 0, 0, CALPHA, A, B, 1 );
      chkxer('ZTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztfsm('N', '/', 'U', 'C', 'U', 0, 0, CALPHA, A, B, 1 );
      chkxer('ZTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      ztfsm('N', 'L', '/', 'C', 'U', 0, 0, CALPHA, A, B, 1 );
      chkxer('ZTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 4
      ztfsm('N', 'L', 'U', '/', 'U', 0, 0, CALPHA, A, B, 1 );
      chkxer('ZTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      ztfsm('N', 'L', 'U', 'C', '/', 0, 0, CALPHA, A, B, 1 );
      chkxer('ZTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 6
      ztfsm('N', 'L', 'U', 'C', 'U', -1, 0, CALPHA, A, B, 1 );
      chkxer('ZTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 7
      ztfsm('N', 'L', 'U', 'C', 'U', 0, -1, CALPHA, A, B, 1 );
      chkxer('ZTFSM ', INFOT, NOUT, LERR, OK );
      INFOT = 11
      ztfsm('N', 'L', 'U', 'C', 'U', 0, 0, CALPHA, A, B, 0 );
      chkxer('ZTFSM ', INFOT, NOUT, LERR, OK );

      SRNAMT = 'ZTFTRI'
      INFOT = 1
      ztftri('/', 'L', 'N', 0, A, INFO );
      chkxer('ZTFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztftri('N', '/', 'N', 0, A, INFO );
      chkxer('ZTFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 3
      ztftri('N', 'L', '/', 0, A, INFO );
      chkxer('ZTFTRI', INFOT, NOUT, LERR, OK );
      INFOT = 4
      ztftri('N', 'L', 'N', -1, A, INFO );
      chkxer('ZTFTRI', INFOT, NOUT, LERR, OK );

      SRNAMT = 'ZTFTTR'
      INFOT = 1
      ztfttr('/', 'U', 0, A, B, 1, INFO );
      chkxer('ZTFTTR', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztfttr('N', '/', 0, A, B, 1, INFO );
      chkxer('ZTFTTR', INFOT, NOUT, LERR, OK );
      INFOT = 3
      ztfttr('N', 'U', -1, A, B, 1, INFO );
      chkxer('ZTFTTR', INFOT, NOUT, LERR, OK );
      INFOT = 6
      ztfttr('N', 'U', 0, A, B, 0, INFO );
      chkxer('ZTFTTR', INFOT, NOUT, LERR, OK );

      SRNAMT = 'ZTRTTF'
      INFOT = 1
      ztrttf('/', 'U', 0, A, 1, B, INFO );
      chkxer('ZTRTTF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztrttf('N', '/', 0, A, 1, B, INFO );
      chkxer('ZTRTTF', INFOT, NOUT, LERR, OK );
      INFOT = 3
      ztrttf('N', 'U', -1, A, 1, B, INFO );
      chkxer('ZTRTTF', INFOT, NOUT, LERR, OK );
      INFOT = 5
      ztrttf('N', 'U', 0, A, 0, B, INFO );
      chkxer('ZTRTTF', INFOT, NOUT, LERR, OK );

      SRNAMT = 'ZTFTTP'
      INFOT = 1
      ztfttp('/', 'U', 0, A, B, INFO );
      chkxer('ZTFTTP', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztfttp('N', '/', 0, A, B, INFO );
      chkxer('ZTFTTP', INFOT, NOUT, LERR, OK );
      INFOT = 3
      ztfttp('N', 'U', -1, A, B, INFO );
      chkxer('ZTFTTP', INFOT, NOUT, LERR, OK );

      SRNAMT = 'ZTPTTF'
      INFOT = 1
      ztpttf('/', 'U', 0, A, B, INFO );
      chkxer('ZTPTTF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztpttf('N', '/', 0, A, B, INFO );
      chkxer('ZTPTTF', INFOT, NOUT, LERR, OK );
      INFOT = 3
      ztpttf('N', 'U', -1, A, B, INFO );
      chkxer('ZTPTTF', INFOT, NOUT, LERR, OK );

      SRNAMT = 'ZTRTTP'
      INFOT = 1
      ztrttp('/', 0, A, 1,  B, INFO );
      chkxer('ZTRTTP', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztrttp('U', -1, A, 1,  B, INFO );
      chkxer('ZTRTTP', INFOT, NOUT, LERR, OK );
      INFOT = 4
      ztrttp('U', 0, A, 0,  B, INFO );
      chkxer('ZTRTTP', INFOT, NOUT, LERR, OK );

      SRNAMT = 'ZTPTTR'
      INFOT = 1
      ztpttr('/', 0, A, B, 1,  INFO );
      chkxer('ZTPTTR', INFOT, NOUT, LERR, OK );
      INFOT = 2
      ztpttr('U', -1, A, B, 1,  INFO );
      chkxer('ZTPTTR', INFOT, NOUT, LERR, OK );
      INFOT = 5
      ztpttr('U', 0, A, B, 0, INFO );
      chkxer('ZTPTTR', INFOT, NOUT, LERR, OK );

      SRNAMT = 'ZHFRK '
      INFOT = 1
      zhfrk('/', 'U', 'N', 0, 0, ALPHA, A, 1, BETA, B );
      chkxer('ZHFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 2
      zhfrk('N', '/', 'N', 0, 0, ALPHA, A, 1, BETA, B );
      chkxer('ZHFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 3
      zhfrk('N', 'U', '/', 0, 0, ALPHA, A, 1, BETA, B );
      chkxer('ZHFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 4
      zhfrk('N', 'U', 'N', -1, 0, ALPHA, A, 1, BETA, B );
      chkxer('ZHFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 5
      zhfrk('N', 'U', 'N', 0, -1, ALPHA, A, 1, BETA, B );
      chkxer('ZHFRK ', INFOT, NOUT, LERR, OK );
      INFOT = 8
      zhfrk('N', 'U', 'N', 0, 0, ALPHA, A, 0, BETA, B );
      chkxer('ZHFRK ', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )
      } else {
         WRITE( NOUT, FMT = 9998 )
      }

 9999 FORMAT( 1X, 'COMPLEX*16 RFP routines passed the tests of the ', 'error exits' )
 9998 FORMAT( ' *** RFP routines failed the tests of the error ', 'exits ***' )
      RETURN

      // End of ZERRRFP

      }
