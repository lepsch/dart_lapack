      void zerrql(PATH, NUNIT ) {

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
      Complex         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGEQL2, ZGEQLF, ZGEQLS, ZUNG2L, ZUNGQL, ZUNM2L, ZUNMQL
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
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) );
         } // 10
         B( J ) = 0.0;
         W( J ) = 0.0;
         X( J ) = 0.0;
      } // 20
      OK = true;

      // Error exits for QL factorization

      // ZGEQLF

      SRNAMT = 'ZGEQLF';
      INFOT = 1;
      zgeqlf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('ZGEQLF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgeqlf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('ZGEQLF', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zgeqlf(2, 1, A, 1, B, W, 1, INFO );
      chkxer('ZGEQLF', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zgeqlf(1, 2, A, 1, B, W, 1, INFO );
      chkxer('ZGEQLF', INFOT, NOUT, LERR, OK );

      // ZGEQL2

      SRNAMT = 'ZGEQL2';
      INFOT = 1;
      zgeql2(-1, 0, A, 1, B, W, INFO );
      chkxer('ZGEQL2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgeql2(0, -1, A, 1, B, W, INFO );
      chkxer('ZGEQL2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zgeql2(2, 1, A, 1, B, W, INFO );
      chkxer('ZGEQL2', INFOT, NOUT, LERR, OK );

      // ZGEQLS

      SRNAMT = 'ZGEQLS';
      INFOT = 1;
      zgeqls(-1, 0, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('ZGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgeqls(0, -1, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('ZGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zgeqls(1, 2, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('ZGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zgeqls(0, 0, -1, A, 1, X, B, 1, W, 1, INFO );
      chkxer('ZGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zgeqls(2, 1, 0, A, 1, X, B, 2, W, 1, INFO );
      chkxer('ZGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      zgeqls(2, 1, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('ZGEQLS', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      zgeqls(1, 1, 2, A, 1, X, B, 1, W, 1, INFO );
      chkxer('ZGEQLS', INFOT, NOUT, LERR, OK );

      // ZUNGQL

      SRNAMT = 'ZUNGQL';
      INFOT = 1;
      zungql(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQL', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zungql(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQL', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zungql(1, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('ZUNGQL', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zungql(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQL', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zungql(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zungql(2, 1, 0, A, 1, X, W, 1, INFO );
      chkxer('ZUNGQL', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      zungql(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('ZUNGQL', INFOT, NOUT, LERR, OK );

      // ZUNG2L

      SRNAMT = 'ZUNG2L';
      INFOT = 1;
      zung2l(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('ZUNG2L', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zung2l(0, -1, 0, A, 1, X, W, INFO );
      chkxer('ZUNG2L', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zung2l(1, 2, 0, A, 1, X, W, INFO );
      chkxer('ZUNG2L', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zung2l(0, 0, -1, A, 1, X, W, INFO );
      chkxer('ZUNG2L', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zung2l(2, 1, 2, A, 2, X, W, INFO );
      chkxer('ZUNG2L', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zung2l(2, 1, 0, A, 1, X, W, INFO );
      chkxer('ZUNG2L', INFOT, NOUT, LERR, OK );

      // ZUNMQL

      SRNAMT = 'ZUNMQL';
      INFOT = 1;
      zunmql('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zunmql('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zunmql('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zunmql('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zunmql('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zunmql('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zunmql('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zunmql('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('ZUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zunmql('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      zunmql('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      zunmql('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMQL', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      zunmql('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('ZUNMQL', INFOT, NOUT, LERR, OK );

      // ZUNM2L

      SRNAMT = 'ZUNM2L';
      INFOT = 1;
      zunm2l('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      zunm2l('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      zunm2l('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      zunm2l('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zunm2l('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zunm2l('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      zunm2l('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zunm2l('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO );
      chkxer('ZUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      zunm2l('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNM2L', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      zunm2l('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('ZUNM2L', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of ZERRQL

      }
