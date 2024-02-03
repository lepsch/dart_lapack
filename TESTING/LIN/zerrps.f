      SUBROUTINE ZERRPS( PATH, NUNIT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NUNIT;
      String             PATH;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      PARAMETER          ( NMAX = 4 )
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, RANK;
      // ..
      // .. Local Arrays ..
      COMPLEX*16         A( NMAX, NMAX )
      double             RWORK( 2*NMAX );
      int                PIV( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZPSTF2, ZPSTRF
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
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      DO 110 J = 1, NMAX
         DO 100 I = 1, NMAX
            A( I, J ) = 1.D0 / DBLE( I+J )

  100    CONTINUE
         PIV( J ) = J
         RWORK( J ) = 0.D0
         RWORK( NMAX+J ) = 0.D0

  110 CONTINUE
      OK = .TRUE.


         // Test error exits of the routines that use the Cholesky
         // decomposition of an Hermitian positive semidefinite matrix.

         // ZPSTRF

      SRNAMT = 'ZPSTRF'
      INFOT = 1
      CALL ZPSTRF( '/', 0, A, 1, PIV, RANK, -1.D0, RWORK, INFO )
      CALL CHKXER( 'ZPSTRF', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZPSTRF( 'U', -1, A, 1, PIV, RANK, -1.D0, RWORK, INFO )
      CALL CHKXER( 'ZPSTRF', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZPSTRF( 'U', 2, A, 1, PIV, RANK, -1.D0, RWORK, INFO )
      CALL CHKXER( 'ZPSTRF', INFOT, NOUT, LERR, OK )

         // ZPSTF2

      SRNAMT = 'ZPSTF2'
      INFOT = 1
      CALL ZPSTF2( '/', 0, A, 1, PIV, RANK, -1.D0, RWORK, INFO )
      CALL CHKXER( 'ZPSTF2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL ZPSTF2( 'U', -1, A, 1, PIV, RANK, -1.D0, RWORK, INFO )
      CALL CHKXER( 'ZPSTF2', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL ZPSTF2( 'U', 2, A, 1, PIV, RANK, -1.D0, RWORK, INFO )
      CALL CHKXER( 'ZPSTF2', INFOT, NOUT, LERR, OK )


      // Print a summary line.

      CALL ALAESM( PATH, OK, NOUT )

      RETURN

      // End of ZERRPS

      }
