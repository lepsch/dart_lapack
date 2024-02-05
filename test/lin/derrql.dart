import 'common.dart';

      void derrql(PATH, NUNIT ) {

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
      double             A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGEQL2, DGEQLF, DGEQLS, DORG2L, DORGQL, DORM2L, DORMQL
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC / srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      infoc.NOUT = NUNIT;
      WRITE( infoc.NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I, J] = 1.0 / (I+J).toDouble();
            AF[I, J] = 1.0 / (I+J).toDouble();
         } // 10
         B[J] = 0.0;
         W[J] = 0.0;
         X[J] = 0.0;
      } // 20
      infoc.OK = true;

      // Error exits for QL factorization

      // DGEQLF

      srnamc.SRNAMT = 'DGEQLF';
      infoc.INFOT = 1;
      dgeqlf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('DGEQLF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgeqlf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('DGEQLF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgeqlf(2, 1, A, 1, B, W, 1, INFO );
      chkxer('DGEQLF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dgeqlf(1, 2, A, 1, B, W, 1, INFO );
      chkxer('DGEQLF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGEQL2

      srnamc.SRNAMT = 'DGEQL2';
      infoc.INFOT = 1;
      dgeql2(-1, 0, A, 1, B, W, INFO );
      chkxer('DGEQL2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgeql2(0, -1, A, 1, B, W, INFO );
      chkxer('DGEQL2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgeql2(2, 1, A, 1, B, W, INFO );
      chkxer('DGEQL2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGEQLS

      srnamc.SRNAMT = 'DGEQLS';
      infoc.INFOT = 1;
      dgeqls(-1, 0, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGEQLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgeqls(0, -1, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGEQLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgeqls(1, 2, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGEQLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dgeqls(0, 0, -1, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGEQLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dgeqls(2, 1, 0, A, 1, X, B, 2, W, 1, INFO );
      chkxer('DGEQLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dgeqls(2, 1, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('DGEQLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dgeqls(1, 1, 2, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGEQLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DORGQL

      srnamc.SRNAMT = 'DORGQL';
      infoc.INFOT = 1;
      dorgql(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorgql(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorgql(1, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('DORGQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorgql(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('DORGQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorgql(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('DORGQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dorgql(2, 1, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dorgql(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('DORGQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DORG2L

      srnamc.SRNAMT = 'DORG2L';
      infoc.INFOT = 1;
      dorg2l(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('DORG2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorg2l(0, -1, 0, A, 1, X, W, INFO );
      chkxer('DORG2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorg2l(1, 2, 0, A, 1, X, W, INFO );
      chkxer('DORG2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorg2l(0, 0, -1, A, 1, X, W, INFO );
      chkxer('DORG2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorg2l(2, 1, 2, A, 2, X, W, INFO );
      chkxer('DORG2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dorg2l(2, 1, 0, A, 1, X, W, INFO );
      chkxer('DORG2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DORMQL

      srnamc.SRNAMT = 'DORMQL';
      infoc.INFOT = 1;
      dormql('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dormql('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dormql('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dormql('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dormql('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dormql('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dormql('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dormql('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dormql('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dormql('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      dormql('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      dormql('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMQL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DORM2L

      srnamc.SRNAMT = 'DORM2L';
      infoc.INFOT = 1;
      dorm2l('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorm2l('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorm2l('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dorm2l('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dorm2l('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dorm2l('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dorm2l('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dorm2l('L', 'N', 2, 1, 0, A, 1, X, AF, 2, W, INFO );
      chkxer('DORM2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dorm2l('R', 'N', 1, 2, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORM2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dorm2l('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('DORM2L', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      return;
      }
