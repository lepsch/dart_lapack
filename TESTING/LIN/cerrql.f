      void cerrql(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 2 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. Local Arrays ..
      COMPLEX            A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CGEQL2, CGEQLF, CGEQLS, CHKXER, CUNG2L, CUNGQL, CUNM2L, CUNMQL
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
      // INTRINSIC CMPLX, REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
            AF( I, J ) = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
         } // 10
         B( J ) = 0.;
         W( J ) = 0.;
         X( J ) = 0.;
      } // 20
      OK = true;

      // Error exits for QL factorization

      // CGEQLF

      SRNAMT = 'CGEQLF';
      INFOT = 1;
      cgeqlf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('CGEQLF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgeqlf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('CGEQLF', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgeqlf(2, 1, A, 1, B, W, 1, INFO );
      chkxer('CGEQLF', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cgeqlf(1, 2, A, 1, B, W, 1, INFO );
      chkxer('CGEQLF', INFOT, NOUT, LERR, OK );

      // CGEQL2

      SRNAMT = 'CGEQL2';
      INFOT = 1;
      cgeql2(-1, 0, A, 1, B, W, INFO );
      chkxer('CGEQL2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgeql2(0, -1, A, 1, B, W, INFO );
      chkxer('CGEQL2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgeql2(2, 1, A, 1, B, W, INFO );
      chkxer('CGEQL2', INFOT, NOUT, LERR, OK );

      // CGEQLS

      SRNAMT = 'CGEQLS';
      INFOT = 1;
      cgeqls(-1, 0, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('CGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgeqls(0, -1, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('CGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgeqls(1, 2, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('CGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgeqls(0, 0, -1, A, 1, X, B, 1, W, 1, INFO );
      chkxer('CGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgeqls(2, 1, 0, A, 1, X, B, 2, W, 1, INFO );
      chkxer('CGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgeqls(2, 1, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('CGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cgeqls(1, 1, 2, A, 1, X, B, 1, W, 1, INFO );
      chkxer('CGEQLS', INFOT, NOUT, LERR, OK );

      // CUNGQL

      SRNAMT = 'CUNGQL';
      INFOT = 1;
      cungql(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('CUNGQL', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cungql(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('CUNGQL', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cungql(1, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('CUNGQL', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cungql(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('CUNGQL', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cungql(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('CUNGQL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cungql(2, 1, 0, A, 1, X, W, 1, INFO );
      chkxer('CUNGQL', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cungql(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('CUNGQL', INFOT, NOUT, LERR, OK );

      // CUNG2L

      SRNAMT = 'CUNG2L';
      INFOT = 1;
      cung2l(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('CUNG2L', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cung2l(0, -1, 0, A, 1, X, W, INFO );
      chkxer('CUNG2L', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cung2l(1, 2, 0, A, 1, X, W, INFO );
      chkxer('CUNG2L', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cung2l(0, 0, -1, A, 1, X, W, INFO );
      chkxer('CUNG2L', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cung2l(2, 1, 2, A, 2, X, W, INFO );
      chkxer('CUNG2L', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cung2l(2, 1, 0, A, 1, X, W, INFO );
      chkxer('CUNG2L', INFOT, NOUT, LERR, OK );

      // CUNMQL

      SRNAMT = 'CUNMQL';
      INFOT = 1;
      cunmql('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cunmql('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cunmql('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cunmql('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmql('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmql('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmql('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunmql('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('CUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunmql('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cunmql('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cunmql('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cunmql('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('CUNMQL', INFOT, NOUT, LERR, OK );

      // CUNM2L

      SRNAMT = 'CUNM2L';
      INFOT = 1;
      cunm2l('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cunm2l('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cunm2l('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cunm2l('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunm2l('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunm2l('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunm2l('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunm2l('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO );
      chkxer('CUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunm2l('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cunm2l('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('CUNM2L', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of CERRQL

      }
