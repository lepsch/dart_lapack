      SUBROUTINE DERRQL( PATH, NUNIT )

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
      const              NMAX = 2 ;
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
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = 1.D0 / DBLE( I+J )
            AF( I, J ) = 1.D0 / DBLE( I+J )
         } // 10
         B( J ) = 0.D0
         W( J ) = 0.D0
         X( J ) = 0.D0
      } // 20
      OK = .TRUE.

      // Error exits for QL factorization

      // DGEQLF

      SRNAMT = 'DGEQLF'
      INFOT = 1
      dgeqlf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('DGEQLF', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgeqlf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('DGEQLF', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgeqlf(2, 1, A, 1, B, W, 1, INFO );
      chkxer('DGEQLF', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dgeqlf(1, 2, A, 1, B, W, 1, INFO );
      chkxer('DGEQLF', INFOT, NOUT, LERR, OK );

      // DGEQL2

      SRNAMT = 'DGEQL2'
      INFOT = 1
      dgeql2(-1, 0, A, 1, B, W, INFO );
      chkxer('DGEQL2', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgeql2(0, -1, A, 1, B, W, INFO );
      chkxer('DGEQL2', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dgeql2(2, 1, A, 1, B, W, INFO );
      chkxer('DGEQL2', INFOT, NOUT, LERR, OK );

      // DGEQLS

      SRNAMT = 'DGEQLS'
      INFOT = 1
      dgeqls(-1, 0, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgeqls(0, -1, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dgeqls(1, 2, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dgeqls(0, 0, -1, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dgeqls(2, 1, 0, A, 1, X, B, 2, W, 1, INFO );
      chkxer('DGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 8
      dgeqls(2, 1, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('DGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 10
      dgeqls(1, 1, 2, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGEQLS', INFOT, NOUT, LERR, OK );

      // DORGQL

      SRNAMT = 'DORGQL'
      INFOT = 1
      dorgql(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGQL', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorgql(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGQL', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorgql(1, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('DORGQL', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorgql(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('DORGQL', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorgql(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('DORGQL', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorgql(2, 1, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGQL', INFOT, NOUT, LERR, OK );
      INFOT = 8
      dorgql(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('DORGQL', INFOT, NOUT, LERR, OK );

      // DORG2L

      SRNAMT = 'DORG2L'
      INFOT = 1
      dorg2l(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('DORG2L', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorg2l(0, -1, 0, A, 1, X, W, INFO );
      chkxer('DORG2L', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorg2l(1, 2, 0, A, 1, X, W, INFO );
      chkxer('DORG2L', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorg2l(0, 0, -1, A, 1, X, W, INFO );
      chkxer('DORG2L', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorg2l(2, 1, 2, A, 2, X, W, INFO );
      chkxer('DORG2L', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorg2l(2, 1, 0, A, 1, X, W, INFO );
      chkxer('DORG2L', INFOT, NOUT, LERR, OK );

      // DORMQL

      SRNAMT = 'DORMQL'
      INFOT = 1
      dormql('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dormql('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dormql('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dormql('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormql('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormql('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dormql('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dormql('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dormql('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 10
      dormql('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 12
      dormql('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 12
      dormql('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMQL', INFOT, NOUT, LERR, OK );

      // DORM2L

      SRNAMT = 'DORM2L'
      INFOT = 1
      dorm2l('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 2
      dorm2l('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 3
      dorm2l('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 4
      dorm2l('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorm2l('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorm2l('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 5
      dorm2l('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dorm2l('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO );
      chkxer('DORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 7
      dorm2l('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 10
      dorm2l('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('DORM2L', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of DERRQL

      }
