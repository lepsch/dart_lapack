      SUBROUTINE ZERRUNHR_COL( PATH, NUNIT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String   (LEN=3)   PATH;
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
      COMPLEX*16         A( NMAX, NMAX ), T( NMAX, NMAX ), D(NMAX)
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZUNHR_COL
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String   (LEN=32)  SRNAMT;
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
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
*
      // Set the variables to innocuous values.
*
      DO J = 1, NMAX
         DO I = 1, NMAX
            A( I, J ) = DCMPLX( 1.D+0 / DBLE( I+J ) )
            T( I, J ) = DCMPLX( 1.D+0 / DBLE( I+J ) )
         END DO
         D( J ) = ( 0.D+0, 0.D+0 )
      END DO
      OK = .TRUE.
*
      // Error exits for Householder reconstruction
*
      // ZUNHR_COL
*
      SRNAMT = 'ZUNHR_COL'
*
      INFOT = 1
      CALL ZUNHR_COL( -1, 0, 1, A, 1, T, 1, D, INFO )
      CALL CHKXER( 'ZUNHR_COL', INFOT, NOUT, LERR, OK )
*
      INFOT = 2
      CALL ZUNHR_COL( 0, -1, 1, A, 1, T, 1, D, INFO )
      CALL CHKXER( 'ZUNHR_COL', INFOT, NOUT, LERR, OK )
      CALL ZUNHR_COL( 1, 2, 1, A, 1, T, 1, D, INFO )
      CALL CHKXER( 'ZUNHR_COL', INFOT, NOUT, LERR, OK )
*
      INFOT = 3
      CALL ZUNHR_COL( 0, 0, -1, A, 1, T, 1, D, INFO )
      CALL CHKXER( 'ZUNHR_COL', INFOT, NOUT, LERR, OK )
*
      CALL ZUNHR_COL( 0, 0, 0, A, 1, T, 1, D, INFO )
      CALL CHKXER( 'ZUNHR_COL', INFOT, NOUT, LERR, OK )
*
      INFOT = 5
      CALL ZUNHR_COL( 0, 0, 1, A, -1, T, 1, D, INFO )
      CALL CHKXER( 'ZUNHR_COL', INFOT, NOUT, LERR, OK )
*
      CALL ZUNHR_COL( 0, 0, 1, A, 0, T, 1, D, INFO )
      CALL CHKXER( 'ZUNHR_COL', INFOT, NOUT, LERR, OK )
*
      CALL ZUNHR_COL( 2, 0, 1, A, 1, T, 1, D, INFO )
      CALL CHKXER( 'ZUNHR_COL', INFOT, NOUT, LERR, OK )
*
      INFOT = 7
      CALL ZUNHR_COL( 0, 0, 1, A, 1, T, -1, D, INFO )
      CALL CHKXER( 'ZUNHR_COL', INFOT, NOUT, LERR, OK )
*
      CALL ZUNHR_COL( 0, 0, 1, A, 1, T, 0, D, INFO )
      CALL CHKXER( 'ZUNHR_COL', INFOT, NOUT, LERR, OK )
*
      CALL ZUNHR_COL( 4, 3, 2, A, 4, T, 1, D, INFO )
      CALL CHKXER( 'ZUNHR_COL', INFOT, NOUT, LERR, OK )
*
      // Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
      // End of ZERRUNHR_COL
*
      END
