      SUBROUTINE SERRLQT( PATH, NUNIT )
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
      REAL               A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SGELQT3, SGELQT, SGEMLQT
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      DO J = 1, NMAX
         DO I = 1, NMAX
            A( I, J ) = 1. / REAL( I+J )
            C( I, J ) = 1. / REAL( I+J )
            T( I, J ) = 1. / REAL( I+J )
         END DO
         W( J ) = 0.
      END DO
      OK = .TRUE.

      // Error exits for LQT factorization

      // SGELQT

      SRNAMT = 'SGELQT'
      INFOT = 1
      sgelqt(-1, 0, 1, A, 1, T, 1, W, INFO );
      chkxer('SGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sgelqt(0, -1, 1, A, 1, T, 1, W, INFO );
      chkxer('SGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sgelqt(0, 0, 0, A, 1, T, 1, W, INFO );
      chkxer('SGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sgelqt(2, 1, 1, A, 1, T, 1, W, INFO );
      chkxer('SGELQT', INFOT, NOUT, LERR, OK );
      INFOT = 7
      sgelqt(2, 2, 2, A, 2, T, 1, W, INFO );
      chkxer('SGELQT', INFOT, NOUT, LERR, OK );

      // SGELQT3

      SRNAMT = 'SGELQT3'
      INFOT = 1
      sgelqt3(-1, 0, A, 1, T, 1, INFO );
      chkxer('SGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sgelqt3(0, -1, A, 1, T, 1, INFO );
      chkxer('SGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 4
      sgelqt3(2, 2, A, 1, T, 1, INFO );
      chkxer('SGELQT3', INFOT, NOUT, LERR, OK );
      INFOT = 6
      sgelqt3(2, 2, A, 2, T, 1, INFO );
      chkxer('SGELQT3', INFOT, NOUT, LERR, OK );

      // SGEMLQT

      SRNAMT = 'SGEMLQT'
      INFOT = 1
      sgemlqt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2
      sgemlqt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3
      sgemlqt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4
      sgemlqt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sgemlqt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5
      sgemlqt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 6
      sgemlqt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8
      sgemlqt('R', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8
      sgemlqt('L', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 10
      sgemlqt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 12
      sgemlqt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO );
      chkxer('SGEMLQT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of SERRLQT

      }
