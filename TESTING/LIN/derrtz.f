      SUBROUTINE DERRTZ( PATH, NUNIT )

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
      String             C2;
      int                INFO;
      // ..
      // .. Local Arrays ..
      double             A( NMAX, NMAX ), TAU( NMAX ), W( NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DTZRZF
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
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
      A( 1, 1 ) = 1.D+0
      A( 1, 2 ) = 2.D+0
      A( 2, 2 ) = 3.D+0
      A( 2, 1 ) = 4.D+0
      W( 1 ) = 0.0D+0
      W( 2 ) = 0.0D+0
      OK = .TRUE.

      IF( LSAMEN( 2, C2, 'TZ' ) ) THEN

         // Test error exits for the trapezoidal routines.

         // DTZRZF

         SRNAMT = 'DTZRZF'
         INFOT = 1
         CALL DTZRZF( -1, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'DTZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 2
         CALL DTZRZF( 1, 0, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'DTZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 4
         CALL DTZRZF( 2, 2, A, 1, TAU, W, 1, INFO )
         CALL CHKXER( 'DTZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DTZRZF( 2, 2, A, 2, TAU, W, 0, INFO )
         CALL CHKXER( 'DTZRZF', INFOT, NOUT, LERR, OK )
         INFOT = 7
         CALL DTZRZF( 2, 3, A, 2, TAU, W, 1, INFO )
         CALL CHKXER( 'DTZRZF', INFOT, NOUT, LERR, OK )
      END IF

      // Print a summary line.

      CALL ALAESM( PATH, OK, NOUT )

      RETURN

      // End of DERRTZ

      }
