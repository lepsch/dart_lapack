      SUBROUTINE SERRTZ( PATH, NUNIT )
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
      String             C2;
      int                INFO
*     ..
*     .. Local Arrays ..
      REAL               A( NMAX, NMAX ), TAU( NMAX ), W( NMAX )
*     ..
*     .. External Functions ..
      bool               LSAMEN;
      EXTERNAL           LSAMEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, STZRZF
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
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
      A( 1, 1 ) = 1.E+0
      A( 1, 2 ) = 2.E+0
      A( 2, 2 ) = 3.E+0
      A( 2, 1 ) = 4.E+0
      W( 1 ) = 0.0E+0
      W( 2 ) = 0.0E+0
      OK = .TRUE.
*
      IF( LSAMEN( 2, C2, 'TZ' ) ) THEN
*
*        Test error exits for the trapezoidal routines.
*
*        STZRZF
*
         SRNAMT = 'STZRZF'
         INFOT = 1
         CALL STZRZF( -1, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'STZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL STZRZF( 1, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'STZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL STZRZF( 2, 2, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'STZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL STZRZF( 2, 2, A, 2, TAU, W, 0, INFO )
         CALL CHKXER( 'STZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL STZRZF( 2, 3, A, 2, TAU, W, 1, INFO )
         CALL CHKXER( 'STZRZF', INFOT, NOUT, LERR, OK )
      END IF
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of SERRTZ
*
      END
