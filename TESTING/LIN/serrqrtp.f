      SUBROUTINE SERRQRTP( PATH, NUNIT )
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
      REAL               A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, STPQRT2, STPQRT, STPMQRT
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      String             SRNAMT;
      int                INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          FLOAT
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
            A( I, J ) = 1.0 / FLOAT( I+J )
            C( I, J ) = 1.0 / FLOAT( I+J )
            T( I, J ) = 1.0 / FLOAT( I+J )
         END DO
         W( J ) = 0.0
      END DO
      OK = .TRUE.
*
*     Error exits for TPQRT factorization
*
*     STPQRT
*
      SRNAMT = 'STPQRT'
      INFOT = 1
      CALL STPQRT( -1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'STPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL STPQRT( 1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'STPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL STPQRT( 0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'STPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL STPQRT( 0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'STPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL STPQRT( 0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'STPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL STPQRT( 0, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'STPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 6
      CALL STPQRT( 1, 2, 0, 2, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'STPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL STPQRT( 2, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO )
      CALL CHKXER( 'STPQRT', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL STPQRT( 2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO )
      CALL CHKXER( 'STPQRT', INFOT, NOUT, LERR, OK )
*
*     STPQRT2
*
      SRNAMT = 'STPQRT2'
      INFOT = 1
      CALL STPQRT2( -1, 0, 0, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'STPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL STPQRT2( 0, -1, 0, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'STPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL STPQRT2( 0, 0, -1, A, 1, B, 1, T, 1, INFO )
      CALL CHKXER( 'STPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL STPQRT2( 2, 2, 0, A, 1, B, 2, T, 2, INFO )
      CALL CHKXER( 'STPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL STPQRT2( 2, 2, 0, A, 2, B, 1, T, 2, INFO )
      CALL CHKXER( 'STPQRT2', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL STPQRT2( 2, 2, 0, A, 2, B, 2, T, 1, INFO )
      CALL CHKXER( 'STPQRT2', INFOT, NOUT, LERR, OK )
*
*     STPMQRT
*
      SRNAMT = 'STPMQRT'
      INFOT = 1
      CALL STPMQRT( '/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'STPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL STPMQRT( 'L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'STPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL STPMQRT( 'L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'STPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL STPMQRT( 'L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'STPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL STPMQRT( 'L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      INFOT = 6
      CALL STPMQRT( 'L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'STPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL STPMQRT( 'L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'STPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL STPMQRT( 'R', 'N', 1, 2, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'STPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 9
      CALL STPMQRT( 'L', 'N', 2, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'STPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 11
      CALL STPMQRT( 'R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO )
      CALL CHKXER( 'STPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 13
      CALL STPMQRT( 'L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO )
      CALL CHKXER( 'STPMQRT', INFOT, NOUT, LERR, OK )
      INFOT = 15
      CALL STPMQRT( 'L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO )
      CALL CHKXER( 'STPMQRT', INFOT, NOUT, LERR, OK )
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of SERRQRTP
*
      END
