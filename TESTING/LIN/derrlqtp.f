      SUBROUTINE DERRLQTP( PATH, NUNIT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      int                NMAX;
      PARAMETER          ( NMAX = 2 )
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
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
*
      // Set the variables to innocuous values.
*
      DO J = 1, NMAX
         DO I = 1, NMAX
            A( I, J ) = 1.D0 / DBLE( I+J )
            C( I, J ) = 1.D0 / DBLE( I+J )
            T( I, J ) = 1.D0 / DBLE( I+J )
         END DO
         W( J ) = 0.0
      END DO
      OK = .TRUE.
*
      // Error exits for TPLQT factorization
*
      // DTPLQT
*
      SRNAMT = 'DTPLQT'
      INFOT = 1
      CALL DTPLQT( -1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'DTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DTPLQT( 1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'DTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DTPLQT( 0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'DTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DTPLQT( 0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'DTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DTPLQT( 0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'DTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DTPLQT( 1, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'DTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL DTPLQT( 2, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'DTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL DTPLQT( 2, 1, 0, 1, A, 2, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'DTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL DTPLQT( 2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO )
      CALL CHKXER( 'DTPLQT', INFOT, NOUT, LERR, OK )
*
      // DTPLQT2
*
      SRNAMT = 'DTPLQT2'
      INFOT = 1
      CALL DTPLQT2( -1, 0, 0, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'DTPLQT2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DTPLQT2( 0, -1, 0, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'DTPLQT2', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DTPLQT2( 0, 0, -1, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'DTPLQT2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DTPLQT2( 2, 2, 0, A, 1, B, 2, T, 2, INFO )
      CALL CHKXER( 'DTPLQT2', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DTPLQT2( 2, 2, 0, A, 2, B, 1, T, 2, INFO )
      CALL CHKXER( 'DTPLQT2', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL DTPLQT2( 2, 2, 0, A, 2, B, 2, T, 1, INFO )
      CALL CHKXER( 'DTPLQT2', INFOT, NOUT, LERR, OK )
*
      // DTPMLQT
*
      SRNAMT = 'DTPMLQT'
      INFOT = 1
      CALL DTPMLQT( '/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'DTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DTPMLQT( 'L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'DTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DTPMLQT( 'L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'DTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DTPMLQT( 'L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'DTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DTPMLQT( 'L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      INFOT = 6
      CALL DTPMLQT( 'L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'DTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DTPMLQT( 'L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'DTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL DTPMLQT( 'R', 'N', 2, 2, 2, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'DTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 11
      CALL DTPMLQT( 'R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'DTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 13
      CALL DTPMLQT( 'L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO )
      CALL CHKXER( 'DTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 15
      CALL DTPMLQT( 'L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO )
      CALL CHKXER( 'DTPMLQT', INFOT, NOUT, LERR, OK )
*
      // Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
      // End of DERRLQTP
*
      END
