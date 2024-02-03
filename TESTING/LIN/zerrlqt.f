      SUBROUTINE ZERRLQT( PATH, NUNIT )
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
      PARAMETER          ( NMAX = 2 )
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      COMPLEX*16   A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGELQT3, ZGELQT, ZGEMLQT
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
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      DO J = 1, NMAX
         DO I = 1, NMAX
            A( I, J ) = 1.D0 / DCMPLX( DBLE( I+J ), 0.D0 )
            C( I, J ) = 1.D0 / DCMPLX( DBLE( I+J ), 0.D0 )
            T( I, J ) = 1.D0 / DCMPLX( DBLE( I+J ), 0.D0 )
         END DO
         W( J ) = 0.D0
      END DO
      OK = .TRUE.

      // Error exits for LQT factorization

      // ZGELQT

      SRNAMT = 'ZGELQT'
      INFOT = 1
      CALL ZGELQT( -1, 0, 1, A, 1, T, 1, W, INFO )
      CALL CHKXER( 'ZGELQT', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZGELQT( 0, -1, 1, A, 1, T, 1, W, INFO )
      CALL CHKXER( 'ZGELQT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL ZGELQT( 0, 0, 0, A, 1, T, 1, W, INFO )
      CALL CHKXER( 'ZGELQT', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL ZGELQT( 2, 1, 1, A, 1, T, 1, W, INFO )
      CALL CHKXER( 'ZGELQT', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL ZGELQT( 2, 2, 2, A, 2, T, 1, W, INFO )
      CALL CHKXER( 'ZGELQT', INFOT, NOUT, LERR, OK )

      // ZGELQT3

      SRNAMT = 'ZGELQT3'
      INFOT = 1
      CALL ZGELQT3( -1, 0, A, 1, T, 1, INFO )
      CALL CHKXER( 'ZGELQT3', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZGELQT3( 0, -1, A, 1, T, 1, INFO )
      CALL CHKXER( 'ZGELQT3', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZGELQT3( 2, 2, A, 1, T, 1, INFO )
      CALL CHKXER( 'ZGELQT3', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL ZGELQT3( 2, 2, A, 2, T, 1, INFO )
      CALL CHKXER( 'ZGELQT3', INFOT, NOUT, LERR, OK )

      // ZGEMLQT

      SRNAMT = 'ZGEMLQT'
      INFOT = 1
      CALL ZGEMLQT( '/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZGEMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZGEMLQT( 'L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZGEMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL ZGEMLQT( 'L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZGEMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZGEMLQT( 'L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZGEMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL ZGEMLQT( 'L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZGEMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL ZGEMLQT( 'R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZGEMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL ZGEMLQT( 'L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZGEMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL ZGEMLQT( 'R', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZGEMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL ZGEMLQT( 'L', 'N', 2, 2, 2, 1, A, 1, T, 1, C, 1, W, INFO )
      CALL CHKXER( 'ZGEMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL ZGEMLQT( 'R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO )
      CALL CHKXER( 'ZGEMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 12
      CALL ZGEMLQT( 'L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO )
      CALL CHKXER( 'ZGEMLQT', INFOT, NOUT, LERR, OK )

      // Print a summary line.

      CALL ALAESM( PATH, OK, NOUT )

      RETURN

      // End of ZERRLQT

      END
