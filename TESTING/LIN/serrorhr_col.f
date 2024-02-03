      SUBROUTINE SERRORHR_COL( PATH, NUNIT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER(LEN=3)   PATH
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
      REAL               A( NMAX, NMAX ), T( NMAX, NMAX ), D(NMAX)
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, SORHR_COL
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER(LEN=32)  SRNAMT
      int                INFOT, NOUT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
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
            A( I, J ) = 1.E+0 / REAL( I+J )
            T( I, J ) = 1.E+0 / REAL( I+J )
         END DO
         D( J ) = 0.E+0
      END DO
      OK = .TRUE.
*
*     Error exits for Householder reconstruction
*
*     SORHR_COL
*
      SRNAMT = 'SORHR_COL'
*
      INFOT = 1
      CALL SORHR_COL( -1, 0, 1, A, 1, T, 1, D, INFO )
      CALL CHKXER( 'SORHR_COL', INFOT, NOUT, LERR, OK )
*
      INFOT = 2
      CALL SORHR_COL( 0, -1, 1, A, 1, T, 1, D, INFO )
      CALL CHKXER( 'SORHR_COL', INFOT, NOUT, LERR, OK )
      CALL SORHR_COL( 1, 2, 1, A, 1, T, 1, D, INFO )
      CALL CHKXER( 'SORHR_COL', INFOT, NOUT, LERR, OK )
*
      INFOT = 3
      CALL SORHR_COL( 0, 0, -1, A, 1, T, 1, D, INFO )
      CALL CHKXER( 'SORHR_COL', INFOT, NOUT, LERR, OK )
*
      CALL SORHR_COL( 0, 0, 0, A, 1, T, 1, D, INFO )
      CALL CHKXER( 'SORHR_COL', INFOT, NOUT, LERR, OK )
*
      INFOT = 5
      CALL SORHR_COL( 0, 0, 1, A, -1, T, 1, D, INFO )
      CALL CHKXER( 'SORHR_COL', INFOT, NOUT, LERR, OK )
*
      CALL SORHR_COL( 0, 0, 1, A, 0, T, 1, D, INFO )
      CALL CHKXER( 'SORHR_COL', INFOT, NOUT, LERR, OK )
*
      CALL SORHR_COL( 2, 0, 1, A, 1, T, 1, D, INFO )
      CALL CHKXER( 'SORHR_COL', INFOT, NOUT, LERR, OK )
*
      INFOT = 7
      CALL SORHR_COL( 0, 0, 1, A, 1, T, -1, D, INFO )
      CALL CHKXER( 'SORHR_COL', INFOT, NOUT, LERR, OK )
*
      CALL SORHR_COL( 0, 0, 1, A, 1, T, 0, D, INFO )
      CALL CHKXER( 'SORHR_COL', INFOT, NOUT, LERR, OK )
*
      CALL SORHR_COL( 4, 3, 2, A, 4, T, 1, D, INFO )
      CALL CHKXER( 'SORHR_COL', INFOT, NOUT, LERR, OK )
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of SERRORHR_COL
*
      END
