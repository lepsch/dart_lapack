      void serrql(PATH, NUNIT ) {

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
      REAL               A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SGEQL2, SGEQLF, SGEQLS, SORG2L, SORGQL, SORM2L, SORMQL
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
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I, J] = 1. / REAL( I+J );
            AF[I, J] = 1. / REAL( I+J );
         } // 10
         B[J] = 0.;
         W[J] = 0.;
         X[J] = 0.;
      } // 20
      OK = true;

      // Error exits for QL factorization

      // SGEQLF

      SRNAMT = 'SGEQLF';
      INFOT = 1;
      sgeqlf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('SGEQLF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgeqlf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('SGEQLF', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgeqlf(2, 1, A, 1, B, W, 1, INFO );
      chkxer('SGEQLF', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sgeqlf(1, 2, A, 1, B, W, 1, INFO );
      chkxer('SGEQLF', INFOT, NOUT, LERR, OK );

      // SGEQL2

      SRNAMT = 'SGEQL2';
      INFOT = 1;
      sgeql2(-1, 0, A, 1, B, W, INFO );
      chkxer('SGEQL2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgeql2(0, -1, A, 1, B, W, INFO );
      chkxer('SGEQL2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sgeql2(2, 1, A, 1, B, W, INFO );
      chkxer('SGEQL2', INFOT, NOUT, LERR, OK );

      // SGEQLS

      SRNAMT = 'SGEQLS';
      INFOT = 1;
      sgeqls(-1, 0, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('SGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgeqls(0, -1, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('SGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sgeqls(1, 2, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('SGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sgeqls(0, 0, -1, A, 1, X, B, 1, W, 1, INFO );
      chkxer('SGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sgeqls(2, 1, 0, A, 1, X, B, 2, W, 1, INFO );
      chkxer('SGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sgeqls(2, 1, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('SGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      sgeqls(1, 1, 2, A, 1, X, B, 1, W, 1, INFO );
      chkxer('SGEQLS', INFOT, NOUT, LERR, OK );

      // SORGQL

      SRNAMT = 'SORGQL';
      INFOT = 1;
      sorgql(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('SORGQL', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sorgql(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('SORGQL', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sorgql(1, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('SORGQL', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sorgql(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('SORGQL', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sorgql(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('SORGQL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sorgql(2, 1, 0, A, 1, X, W, 1, INFO );
      chkxer('SORGQL', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      sorgql(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('SORGQL', INFOT, NOUT, LERR, OK );

      // SORG2L

      SRNAMT = 'SORG2L';
      INFOT = 1;
      sorg2l(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('SORG2L', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sorg2l(0, -1, 0, A, 1, X, W, INFO );
      chkxer('SORG2L', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sorg2l(1, 2, 0, A, 1, X, W, INFO );
      chkxer('SORG2L', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sorg2l(0, 0, -1, A, 1, X, W, INFO );
      chkxer('SORG2L', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sorg2l(2, 1, 2, A, 2, X, W, INFO );
      chkxer('SORG2L', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sorg2l(2, 1, 0, A, 1, X, W, INFO );
      chkxer('SORG2L', INFOT, NOUT, LERR, OK );

      // SORMQL

      SRNAMT = 'SORMQL';
      INFOT = 1;
      sormql('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sormql('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sormql('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sormql('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sormql('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sormql('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sormql('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sormql('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('SORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sormql('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      sormql('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('SORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      sormql('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('SORMQL', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      sormql('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('SORMQL', INFOT, NOUT, LERR, OK );

      // SORM2L

      SRNAMT = 'SORM2L';
      INFOT = 1;
      sorm2l('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      sorm2l('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      sorm2l('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      sorm2l('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sorm2l('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sorm2l('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      sorm2l('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sorm2l('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO );
      chkxer('SORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      sorm2l('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('SORM2L', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      sorm2l('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('SORM2L', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
