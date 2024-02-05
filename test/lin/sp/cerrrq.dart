      void cerrrq(PATH, NUNIT ) {

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
      Complex            A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CGERQ2, CGERQF, CGERQS, CHKXER, CUNGR2, CUNGRQ, CUNMR2, CUNMRQ
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
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
            A[I, J] = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
            AF[I, J] = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
         } // 10
         B[J] = 0.;
         W[J] = 0.;
         X[J] = 0.;
      } // 20
      OK = true;

      // Error exits for RQ factorization

      // CGERQF

     srnamc.SRNAMT = 'CGERQF';
      INFOT = 1;
      cgerqf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('CGERQF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgerqf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('CGERQF', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgerqf(2, 1, A, 1, B, W, 2, INFO );
      chkxer('CGERQF', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cgerqf(2, 1, A, 2, B, W, 1, INFO );
      chkxer('CGERQF', INFOT, NOUT, LERR, OK );

      // CGERQ2

     srnamc.SRNAMT = 'CGERQ2';
      INFOT = 1;
      cgerq2(-1, 0, A, 1, B, W, INFO );
      chkxer('CGERQ2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgerq2(0, -1, A, 1, B, W, INFO );
      chkxer('CGERQ2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgerq2(2, 1, A, 1, B, W, INFO );
      chkxer('CGERQ2', INFOT, NOUT, LERR, OK );

      // CGERQS

     srnamc.SRNAMT = 'CGERQS';
      INFOT = 1;
      cgerqs(-1, 0, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('CGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgerqs(0, -1, 0, A, 1, X, B, 1, W, 1, INFO );
      chkxer('CGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgerqs(2, 1, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('CGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cgerqs(0, 0, -1, A, 1, X, B, 1, W, 1, INFO );
      chkxer('CGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cgerqs(2, 2, 0, A, 1, X, B, 2, W, 1, INFO );
      chkxer('CGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cgerqs(2, 2, 0, A, 2, X, B, 1, W, 1, INFO );
      chkxer('CGERQS', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cgerqs(1, 1, 2, A, 1, X, B, 1, W, 1, INFO );
      chkxer('CGERQS', INFOT, NOUT, LERR, OK );

      // CUNGRQ

     srnamc.SRNAMT = 'CUNGRQ';
      INFOT = 1;
      cungrq(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('CUNGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cungrq(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('CUNGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cungrq(2, 1, 0, A, 2, X, W, 2, INFO );
      chkxer('CUNGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cungrq(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('CUNGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cungrq(1, 2, 2, A, 1, X, W, 1, INFO );
      chkxer('CUNGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cungrq(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('CUNGRQ', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cungrq(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('CUNGRQ', INFOT, NOUT, LERR, OK );

      // CUNGR2

     srnamc.SRNAMT = 'CUNGR2';
      INFOT = 1;
      cungr2(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('CUNGR2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cungr2(0, -1, 0, A, 1, X, W, INFO );
      chkxer('CUNGR2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cungr2(2, 1, 0, A, 2, X, W, INFO );
      chkxer('CUNGR2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cungr2(0, 0, -1, A, 1, X, W, INFO );
      chkxer('CUNGR2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cungr2(1, 2, 2, A, 2, X, W, INFO );
      chkxer('CUNGR2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cungr2(2, 2, 0, A, 1, X, W, INFO );
      chkxer('CUNGR2', INFOT, NOUT, LERR, OK );

      // CUNMRQ

     srnamc.SRNAMT = 'CUNMRQ';
      INFOT = 1;
      cunmrq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cunmrq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cunmrq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cunmrq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmrq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmrq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmrq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunmrq('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('CUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunmrq('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cunmrq('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cunmrq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMRQ', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cunmrq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('CUNMRQ', INFOT, NOUT, LERR, OK );

      // CUNMR2

     srnamc.SRNAMT = 'CUNMR2';
      INFOT = 1;
      cunmr2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cunmr2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cunmr2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cunmr2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmr2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmr2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmr2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunmr2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO );
      chkxer('CUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunmr2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNMR2', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cunmr2('L', 'N', 2, 1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNMR2', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
