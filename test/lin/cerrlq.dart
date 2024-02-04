      void cerrlq(PATH, NUNIT ) {

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
      // EXTERNAL ALAESM, CGELQ2, CGELQF, CHKXER, CUNGL2, CUNGLQ, CUNML2, CUNMLQ
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
            A[I, J] = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
            AF[I, J] = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
         } // 10
         B[J] = 0.;
         W[J] = 0.;
         X[J] = 0.;
      } // 20
      OK = true;

      // Error exits for LQ factorization

      // CGELQF

      SRNAMT = 'CGELQF';
      INFOT = 1;
      cgelqf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('CGELQF', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgelqf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('CGELQF', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgelqf(2, 1, A, 1, B, W, 2, INFO );
      chkxer('CGELQF', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cgelqf(2, 1, A, 2, B, W, 1, INFO );
      chkxer('CGELQF', INFOT, NOUT, LERR, OK );

      // CGELQ2

      SRNAMT = 'CGELQ2';
      INFOT = 1;
      cgelq2(-1, 0, A, 1, B, W, INFO );
      chkxer('CGELQ2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cgelq2(0, -1, A, 1, B, W, INFO );
      chkxer('CGELQ2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cgelq2(2, 1, A, 1, B, W, INFO );
      chkxer('CGELQ2', INFOT, NOUT, LERR, OK );

      // CUNGLQ

      SRNAMT = 'CUNGLQ';
      INFOT = 1;
      cunglq(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('CUNGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cunglq(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('CUNGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cunglq(2, 1, 0, A, 2, X, W, 2, INFO );
      chkxer('CUNGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cunglq(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('CUNGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cunglq(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('CUNGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunglq(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('CUNGLQ', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      cunglq(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('CUNGLQ', INFOT, NOUT, LERR, OK );

      // CUNGL2

      SRNAMT = 'CUNGL2';
      INFOT = 1;
      cungl2(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('CUNGL2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cungl2(0, -1, 0, A, 1, X, W, INFO );
      chkxer('CUNGL2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cungl2(2, 1, 0, A, 2, X, W, INFO );
      chkxer('CUNGL2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cungl2(0, 0, -1, A, 1, X, W, INFO );
      chkxer('CUNGL2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cungl2(1, 1, 2, A, 1, X, W, INFO );
      chkxer('CUNGL2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cungl2(2, 2, 0, A, 1, X, W, INFO );
      chkxer('CUNGL2', INFOT, NOUT, LERR, OK );

      // CUNMLQ

      SRNAMT = 'CUNMLQ';
      INFOT = 1;
      cunmlq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cunmlq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cunmlq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cunmlq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmlq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmlq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunmlq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunmlq('L', 'N', 2, 0, 2, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('CUNMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunmlq('R', 'N', 0, 2, 2, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cunmlq('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('CUNMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cunmlq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('CUNMLQ', INFOT, NOUT, LERR, OK );
      INFOT = 12;
      cunmlq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('CUNMLQ', INFOT, NOUT, LERR, OK );

      // CUNML2

      SRNAMT = 'CUNML2';
      INFOT = 1;
      cunml2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNML2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      cunml2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNML2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      cunml2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNML2', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      cunml2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNML2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunml2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNML2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunml2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNML2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      cunml2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNML2', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunml2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO );
      chkxer('CUNML2', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      cunml2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO );
      chkxer('CUNML2', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      cunml2('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('CUNML2', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
