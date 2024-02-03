      SUBROUTINE DERRLQTP( PATH, NUNIT );
      // IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 2 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      double             A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DTPLQT2, DTPLQT, DTPMLQT
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
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A( I, J ) = 1.0 / DBLE( I+J );
            C( I, J ) = 1.0 / DBLE( I+J );
            T( I, J ) = 1.0 / DBLE( I+J );
         }
         W( J ) = 0.0;
      }
      OK = true;

      // Error exits for TPLQT factorization

      // DTPLQT

      SRNAMT = 'DTPLQT';
      INFOT = 1;
      dtplqt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtplqt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtplqt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtplqt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dtplqt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dtplqt(1, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      dtplqt(2, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      dtplqt(2, 1, 0, 1, A, 2, B, 1, T, 1, W, INFO );
      chkxer('DTPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      dtplqt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO );
      chkxer('DTPLQT', INFOT, NOUT, LERR, OK );

      // DTPLQT2

      SRNAMT = 'DTPLQT2';
      INFOT = 1;
      dtplqt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('DTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtplqt2(0, -1, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('DTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtplqt2(0, 0, -1, A, 1, B, 1, T, 1, INFO );
      chkxer('DTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      dtplqt2(2, 2, 0, A, 1, B, 2, T, 2, INFO );
      chkxer('DTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      dtplqt2(2, 2, 0, A, 2, B, 1, T, 2, INFO );
      chkxer('DTPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      dtplqt2(2, 2, 0, A, 2, B, 2, T, 1, INFO );
      chkxer('DTPLQT2', INFOT, NOUT, LERR, OK );

      // DTPMLQT

      SRNAMT = 'DTPMLQT';
      INFOT = 1;
      dtpmlqt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      dtpmlqt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      dtpmlqt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      dtpmlqt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      dtpmlqt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      INFOT = 6;
      dtpmlqt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      dtpmlqt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      dtpmlqt('R', 'N', 2, 2, 2, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      dtpmlqt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO );
      chkxer('DTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      dtpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO );
      chkxer('DTPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 15;
      dtpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO );
      chkxer('DTPMLQT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN;

      // End of DERRLQTP

      }
