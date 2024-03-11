      void zerrlq(PATH, infoc.NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                infoc.NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      Complex         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGELQ2, ZGELQF, ZUNGL2, ZUNGLQ, ZUNML2, ZUNMLQ
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX

      NOUT = infoc.NUNIT;
      NOUT.println( * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / (I+J).toDouble() );
         } // 10
         B[J] = 0.0;
         W[J] = 0.0;
         X[J] = 0.0;
      } // 20
      infoc.OK.value = true;

      // Error exits for LQ factorization

      // ZGELQF

     srnamc.SRNAMT = 'ZGELQF';
      infoc.INFOT = 1;
      zgelqf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('ZGELQF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgelqf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('ZGELQF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgelqf(2, 1, A, 1, B, W, 2, INFO );
      chkxer('ZGELQF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zgelqf(2, 1, A, 2, B, W, 1, INFO );
      chkxer('ZGELQF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZGELQ2

     srnamc.SRNAMT = 'ZGELQ2';
      infoc.INFOT = 1;
      zgelq2(-1, 0, A, 1, B, W, INFO );
      chkxer('ZGELQ2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgelq2(0, -1, A, 1, B, W, INFO );
      chkxer('ZGELQ2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgelq2(2, 1, A, 1, B, W, INFO );
      chkxer('ZGELQ2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZUNGLQ

     srnamc.SRNAMT = 'ZUNGLQ';
      infoc.INFOT = 1;
      zunglq(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zunglq(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zunglq(2, 1, 0, A, 2, X, W, 2, INFO );
      chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zunglq(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zunglq(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunglq(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      zunglq(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('ZUNGLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZUNGL2

     srnamc.SRNAMT = 'ZUNGL2';
      infoc.INFOT = 1;
      zungl2(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('ZUNGL2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zungl2(0, -1, 0, A, 1, X, W, INFO );
      chkxer('ZUNGL2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zungl2(2, 1, 0, A, 2, X, W, INFO );
      chkxer('ZUNGL2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zungl2(0, 0, -1, A, 1, X, W, INFO );
      chkxer('ZUNGL2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zungl2(1, 1, 2, A, 1, X, W, INFO );
      chkxer('ZUNGL2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zungl2(2, 2, 0, A, 1, X, W, INFO );
      chkxer('ZUNGL2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZUNMLQ

     srnamc.SRNAMT = 'ZUNMLQ';
      infoc.INFOT = 1;
      zunmlq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zunmlq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zunmlq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zunmlq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunmlq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunmlq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunmlq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zunmlq('L', 'N', 2, 0, 2, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zunmlq('R', 'N', 0, 2, 2, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      zunmlq('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      zunmlq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      zunmlq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('ZUNMLQ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZUNML2

     srnamc.SRNAMT = 'ZUNML2';
      infoc.INFOT = 1;
      zunml2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zunml2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zunml2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zunml2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunml2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunml2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zunml2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zunml2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO );
      chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zunml2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO );
      chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      zunml2('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('ZUNML2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      }
