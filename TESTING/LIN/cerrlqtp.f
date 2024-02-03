      SUBROUTINE CERRLQTP( PATH, NUNIT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             PATH;
      int                NUNIT
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                NMAX
      PARAMETER          ( NMAX = 2 )
*     ..
*     .. Local Scalars ..
      int                I, INFO, J
*     ..
*     .. Local Arrays ..
      COMPLEX            A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, CTPLQT2, CTPLQT, CTPMLQT
*     ..
*     .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, CMPLX
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
*
*     Set the variables to innocuous values.
*
      DO J = 1, NMAX
         DO I = 1, NMAX
            A( I, J ) = 1.E0 / CMPLX( REAL( I+J ), 0.E0 )
            C( I, J ) = 1.E0 / CMPLX( REAL( I+J ), 0.E0 )
            T( I, J ) = 1.E0 / CMPLX( REAL( I+J ), 0.E0 )
         END DO
         W( J ) = 0.E0
      END DO
      OK = .TRUE.
*
*     Error exits for TPLQT factorization
*
*     CTPLQT
*
      SRNAMT = 'CTPLQT'
      INFOT = 1
      CALL CTPLQT( -1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CTPLQT( 1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CTPLQT( 0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CTPLQT( 0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CTPLQT( 0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CTPLQT( 1, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL CTPLQT( 2, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL CTPLQT( 2, 1, 0, 1, A, 2, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'CTPLQT', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL CTPLQT( 2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO )
      CALL CHKXER( 'CTPLQT', INFOT, NOUT, LERR, OK )
*
*     CTPLQT2
*
      SRNAMT = 'CTPLQT2'
      INFOT = 1
      CALL CTPLQT2( -1, 0, 0, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'CTPLQT2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CTPLQT2( 0, -1, 0, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'CTPLQT2', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CTPLQT2( 0, 0, -1, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'CTPLQT2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CTPLQT2( 2, 2, 0, A, 1, B, 2, T, 2, INFO )
      CALL CHKXER( 'CTPLQT2', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL CTPLQT2( 2, 2, 0, A, 2, B, 1, T, 2, INFO )
      CALL CHKXER( 'CTPLQT2', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL CTPLQT2( 2, 2, 0, A, 2, B, 2, T, 1, INFO )
      CALL CHKXER( 'CTPLQT2', INFOT, NOUT, LERR, OK )
*
*     CTPMLQT
*
      SRNAMT = 'CTPMLQT'
      INFOT = 1
      CALL CTPMLQT( '/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CTPMLQT( 'L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL CTPMLQT( 'L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CTPMLQT( 'L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL CTPMLQT( 'L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      INFOT = 6
      CALL CTPMLQT( 'L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL CTPMLQT( 'L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL CTPMLQT( 'R', 'N', 2, 2, 2, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 11
      CALL CTPMLQT( 'R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'CTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 13
      CALL CTPMLQT( 'L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO )
      CALL CHKXER( 'CTPMLQT', INFOT, NOUT, LERR, OK )
      INFOT = 15
      CALL CTPMLQT( 'L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO )
      CALL CHKXER( 'CTPMLQT', INFOT, NOUT, LERR, OK )
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of CERRLQTP
*
      END
