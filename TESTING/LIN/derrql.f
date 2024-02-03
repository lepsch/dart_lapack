      SUBROUTINE DERRQL( PATH, NUNIT )
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
      double             A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGEQL2, DGEQLF, DGEQLS, DORG2L, DORGQL, DORM2L, DORMQL
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
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            A( I, J ) = 1.D0 / DBLE( I+J )
            AF( I, J ) = 1.D0 / DBLE( I+J )
   10    CONTINUE
         B( J ) = 0.D0
         W( J ) = 0.D0
         X( J ) = 0.D0
   20 CONTINUE
      OK = .TRUE.
*
      // Error exits for QL factorization
*
      // DGEQLF
*
      SRNAMT = 'DGEQLF'
      INFOT = 1
      CALL DGEQLF( -1, 0, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'DGEQLF', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DGEQLF( 0, -1, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'DGEQLF', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DGEQLF( 2, 1, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'DGEQLF', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DGEQLF( 1, 2, A, 1, B, W, 1, INFO )
      CALL CHKXER( 'DGEQLF', INFOT, NOUT, LERR, OK )
*
      // DGEQL2
*
      SRNAMT = 'DGEQL2'
      INFOT = 1
      CALL DGEQL2( -1, 0, A, 1, B, W, INFO )
      CALL CHKXER( 'DGEQL2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DGEQL2( 0, -1, A, 1, B, W, INFO )
      CALL CHKXER( 'DGEQL2', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DGEQL2( 2, 1, A, 1, B, W, INFO )
      CALL CHKXER( 'DGEQL2', INFOT, NOUT, LERR, OK )
*
      // DGEQLS
*
      SRNAMT = 'DGEQLS'
      INFOT = 1
      CALL DGEQLS( -1, 0, 0, A, 1, X, B, 1, W, 1, INFO )
      CALL CHKXER( 'DGEQLS', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DGEQLS( 0, -1, 0, A, 1, X, B, 1, W, 1, INFO )
      CALL CHKXER( 'DGEQLS', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DGEQLS( 1, 2, 0, A, 1, X, B, 1, W, 1, INFO )
      CALL CHKXER( 'DGEQLS', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DGEQLS( 0, 0, -1, A, 1, X, B, 1, W, 1, INFO )
      CALL CHKXER( 'DGEQLS', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DGEQLS( 2, 1, 0, A, 1, X, B, 2, W, 1, INFO )
      CALL CHKXER( 'DGEQLS', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL DGEQLS( 2, 1, 0, A, 2, X, B, 1, W, 1, INFO )
      CALL CHKXER( 'DGEQLS', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL DGEQLS( 1, 1, 2, A, 1, X, B, 1, W, 1, INFO )
      CALL CHKXER( 'DGEQLS', INFOT, NOUT, LERR, OK )
*
      // DORGQL
*
      SRNAMT = 'DORGQL'
      INFOT = 1
      CALL DORGQL( -1, 0, 0, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGQL', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORGQL( 0, -1, 0, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGQL', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORGQL( 1, 2, 0, A, 1, X, W, 2, INFO )
      CALL CHKXER( 'DORGQL', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORGQL( 0, 0, -1, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGQL', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORGQL( 1, 1, 2, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGQL', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORGQL( 2, 1, 0, A, 1, X, W, 1, INFO )
      CALL CHKXER( 'DORGQL', INFOT, NOUT, LERR, OK )
      INFOT = 8
      CALL DORGQL( 2, 2, 0, A, 2, X, W, 1, INFO )
      CALL CHKXER( 'DORGQL', INFOT, NOUT, LERR, OK )
*
      // DORG2L
*
      SRNAMT = 'DORG2L'
      INFOT = 1
      CALL DORG2L( -1, 0, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'DORG2L', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORG2L( 0, -1, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'DORG2L', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORG2L( 1, 2, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'DORG2L', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORG2L( 0, 0, -1, A, 1, X, W, INFO )
      CALL CHKXER( 'DORG2L', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORG2L( 2, 1, 2, A, 2, X, W, INFO )
      CALL CHKXER( 'DORG2L', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORG2L( 2, 1, 0, A, 1, X, W, INFO )
      CALL CHKXER( 'DORG2L', INFOT, NOUT, LERR, OK )
*
      // DORMQL
*
      SRNAMT = 'DORMQL'
      INFOT = 1
      CALL DORMQL( '/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMQL', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORMQL( 'L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMQL', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORMQL( 'L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMQL', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DORMQL( 'L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMQL', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORMQL( 'L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMQL', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORMQL( 'L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMQL', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORMQL( 'R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMQL', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DORMQL( 'L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO )
      CALL CHKXER( 'DORMQL', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DORMQL( 'R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMQL', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL DORMQL( 'L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMQL', INFOT, NOUT, LERR, OK )
      INFOT = 12
      CALL DORMQL( 'L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO )
      CALL CHKXER( 'DORMQL', INFOT, NOUT, LERR, OK )
      INFOT = 12
      CALL DORMQL( 'R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO )
      CALL CHKXER( 'DORMQL', INFOT, NOUT, LERR, OK )
*
      // DORM2L
*
      SRNAMT = 'DORM2L'
      INFOT = 1
      CALL DORM2L( '/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORM2L', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DORM2L( 'L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORM2L', INFOT, NOUT, LERR, OK )
      INFOT = 3
      CALL DORM2L( 'L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORM2L', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DORM2L( 'L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORM2L', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORM2L( 'L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORM2L', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORM2L( 'L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORM2L', INFOT, NOUT, LERR, OK )
      INFOT = 5
      CALL DORM2L( 'R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORM2L', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DORM2L( 'L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO )
      CALL CHKXER( 'DORM2L', INFOT, NOUT, LERR, OK )
      INFOT = 7
      CALL DORM2L( 'R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORM2L', INFOT, NOUT, LERR, OK )
      INFOT = 10
      CALL DORM2L( 'L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO )
      CALL CHKXER( 'DORM2L', INFOT, NOUT, LERR, OK )
*
      // Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
      // End of DERRQL
*
      END
