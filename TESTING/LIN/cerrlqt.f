      SUBROUTINE CERRLQT( PATH, NUNIT )
      IMPLICIT NONE

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
      COMPLEX            A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CGELQT3, CGELQT, CGEMLQT
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
      // INTRINSIC REAL, CMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A( I, J ) = 1.E0 / CMPLX( REAL( I+J ), 0.E0 )
            C( I, J ) = 1.E0 / CMPLX( REAL( I+J ), 0.E0 )
            T( I, J ) = 1.E0 / CMPLX( REAL( I+J ), 0.E0 )
         }
         W( J ) = 0.E0
      }
      OK = true;

      // Error exits for LQT factorization

      // CGELQT

      SRNAMT = 'CGELQT'
      INFOT = 1
      cgelqt(-1, 0, 1, A, 1, T, 1, W, INFO );
      chkxer('CGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 2
      cgelqt(0, -1, 1, A, 1, T, 1, W, INFO );
      chkxer('CGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 3
      cgelqt(0, 0, 0, A, 1, T, 1, W, INFO );
      chkxer('CGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 5
      cgelqt(2, 1, 1, A, 1, T, 1, W, INFO );
      chkxer('CGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 7
      cgelqt(2, 2, 2, A, 2, T, 1, W, INFO );
      chkxer('CGELQT', INFOT, NOUT, LERR, OK );

      // CGELQT3

      SRNAMT = 'CGELQT3'
      INFOT = 1
      cgelqt3(-1, 0, A, 1, T, 1, INFO );
      chkxer('CGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 2
      cgelqt3(0, -1, A, 1, T, 1, INFO );
      chkxer('CGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 4
      cgelqt3(2, 2, A, 1, T, 1, INFO );
      chkxer('CGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 6
      cgelqt3(2, 2, A, 2, T, 1, INFO );
      chkxer('CGELQT3', INFOT, NOUT, LERR, OK );

      // CGEMLQT

      SRNAMT = 'CGEMLQT'
      INFOT = 1
      cgemlqt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2
      cgemlqt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3
      cgemlqt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4
      cgemlqt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5
      cgemlqt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5
      cgemlqt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 6
      cgemlqt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8
      cgemlqt('R', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8
      cgemlqt('L', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('CGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 10
      cgemlqt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO );
      chkxer('CGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 12
      cgemlqt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO );
      chkxer('CGEMLQT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of CERRLQT

      }
