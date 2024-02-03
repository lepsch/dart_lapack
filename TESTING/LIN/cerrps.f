      SUBROUTINE CERRPS( PATH, NUNIT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                NUNIT;
      String             PATH;
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      int                NMAX;
      PARAMETER          ( NMAX = 4 )
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, RANK;
      // ..
      // .. Local Arrays ..
      COMPLEX            A( NMAX, NMAX )
      REAL               RWORK( 2*NMAX )
      int                PIV( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CPSTF2, CPSTRF
      // ..
      // .. Scalars in Common ..
      int                INFOT, NOUT;
      bool               LERR, OK;
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
*
      // Set the variables to innocuous values.
*
      DO 110 J = 1, NMAX
         DO 100 I = 1, NMAX
            A( I, J ) = 1.0 / REAL( I+J )
*
  100    CONTINUE
         PIV( J ) = J
         RWORK( J ) = 0.
         RWORK( NMAX+J ) = 0.
*
  110 CONTINUE
      OK = .TRUE.
*
*
         // Test error exits of the routines that use the Cholesky
         // decomposition of an Hermitian positive semidefinite matrix.
*
         // CPSTRF
*
      SRNAMT = 'CPSTRF'
      INFOT = 1
      CALL CPSTRF( '/', 0, A, 1, PIV, RANK, -1.0, RWORK, INFO )
      CALL CHKXER( 'CPSTRF', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CPSTRF( 'U', -1, A, 1, PIV, RANK, -1.0, RWORK, INFO )
      CALL CHKXER( 'CPSTRF', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CPSTRF( 'U', 2, A, 1, PIV, RANK, -1.0, RWORK, INFO )
      CALL CHKXER( 'CPSTRF', INFOT, NOUT, LERR, OK )
*
         // CPSTF2
*
      SRNAMT = 'CPSTF2'
      INFOT = 1
      CALL CPSTF2( '/', 0, A, 1, PIV, RANK, -1.0, RWORK, INFO )
      CALL CHKXER( 'CPSTF2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL CPSTF2( 'U', -1, A, 1, PIV, RANK, -1.0, RWORK, INFO )
      CALL CHKXER( 'CPSTF2', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL CPSTF2( 'U', 2, A, 1, PIV, RANK, -1.0, RWORK, INFO )
      CALL CHKXER( 'CPSTF2', INFOT, NOUT, LERR, OK )
*
*
      // Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
      // End of CERRPS
*
      END
