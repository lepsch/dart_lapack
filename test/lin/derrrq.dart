import 'common.dart';
      void derrrq(PATH, NUNIT ) {

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
      // EXTERNAL ALAESM, CHKXER, DGERQ2, DGERQF, DGERQS, DORGR2, DORGRQ, DORMR2, DORMRQ
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
            A[I][J] = 1.0 / (I+J).toDouble();
            AF[I][J] = 1.0 / (I+J).toDouble();
         } // 10
         B[J] = 0.0;
         W[J] = 0.0;
         X[J] = 0.0;
      } // 20
      infoc.OK = true;

      // Error exits for RQ factorization

      // DGERQF

      srnamc.SRNAMT = 'DGERQF';
      infoc.INFOT = 1;
      dgerqf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('DGERQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgerqf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('DGERQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgerqf(2, 1, A, 1, B, W, 2, INFO );
      chkxer('DGERQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dgerqf(2, 1, A, 2, B, W, 1, INFO );
      chkxer('DGERQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGERQ2

      srnamc.SRNAMT = 'DGERQ2';
      infoc.INFOT = 1;
      dgerq2(-1, 0, A, 1, B, W, INFO );
      chkxer('DGERQ2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgerq2(0, -1, A, 1, B, W, INFO );
      chkxer('DGERQ2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgerq2(2, 1, A, 1, B, W, INFO );
      chkxer('DGERQ2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGERQS

      srnamc.SRNAMT = 'DGERQS';
      infoc.INFOT = 1;
      dgerqs(-1, 0, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGERQS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgerqs(0, -1, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGERQS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgerqs(2, 1, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('DGERQS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dgerqs(0, 0, -1, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGERQS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dgerqs(2, 2, 0, A, 1, X, B, 2, W, 1, INFO );
      chkxer('DGERQS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dgerqs(2, 2, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('DGERQS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dgerqs(1, 1, 2, A, 1, X, B, 1, W, 1, INFO );
      chkxer('DGERQS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DORGRQ

      srnamc.SRNAMT = 'DORGRQ';
      infoc.INFOT = 1;
      dorgrq(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorgrq(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorgrq(2, 1, 0, A, 2, X, W, 2, INFO );
      chkxer('DORGRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorgrq(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('DORGRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorgrq(1, 2, 2, A, 1, X, W, 1, INFO );
      chkxer('DORGRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dorgrq(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('DORGRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dorgrq(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('DORGRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DORGR2

      srnamc.SRNAMT = 'DORGR2';
      infoc.INFOT = 1;
      dorgr2(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('DORGR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorgr2(0, -1, 0, A, 1, X, W, INFO );
      chkxer('DORGR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorgr2(2, 1, 0, A, 2, X, W, INFO );
      chkxer('DORGR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorgr2(0, 0, -1, A, 1, X, W, INFO );
      chkxer('DORGR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorgr2(1, 2, 2, A, 2, X, W, INFO );
      chkxer('DORGR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dorgr2(2, 2, 0, A, 1, X, W, INFO );
      chkxer('DORGR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DORMRQ

      srnamc.SRNAMT = 'DORMRQ';
      infoc.INFOT = 1;
      dormrq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dormrq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dormrq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dormrq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dormrq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dormrq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dormrq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dormrq('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dormrq('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dormrq('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      dormrq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      dormrq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMRQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DORMR2

      srnamc.SRNAMT = 'DORMR2';
      infoc.INFOT = 1;
      dormr2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dormr2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dormr2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dormr2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dormr2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dormr2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dormr2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dormr2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO );
      chkxer('DORMR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dormr2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dormr2('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORMR2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      return;
      }
