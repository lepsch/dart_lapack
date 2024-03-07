      void zerrqrt(PATH, infoc.NUNIT ) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                infoc.NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      Complex         A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGEQRT2, ZGEQRT3, ZGEQRT, ZGEMQRT
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
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I][J] = 1.0 / DCMPLX( (I+J).toDouble(), 0.0 );
            C[I][J] = 1.0 / DCMPLX( (I+J).toDouble(), 0.0 );
            T[I][J] = 1.0 / DCMPLX( (I+J).toDouble(), 0.0 );
         }
         W[J] = 0.0;
      }
      infoc.OK.value = true;

      // Error exits for QRT factorization

      // ZGEQRT

     srnamc.SRNAMT = 'ZGEQRT';
      infoc.INFOT = 1;
      zgeqrt(-1, 0, 1, A, 1, T, 1, W, INFO );
      chkxer('ZGEQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgeqrt(0, -1, 1, A, 1, T, 1, W, INFO );
      chkxer('ZGEQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zgeqrt(0, 0, 0, A, 1, T, 1, W, INFO );
      chkxer('ZGEQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zgeqrt(2, 1, 1, A, 1, T, 1, W, INFO );
      chkxer('ZGEQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      zgeqrt(2, 2, 2, A, 2, T, 1, W, INFO );
      chkxer('ZGEQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZGEQRT2

     srnamc.SRNAMT = 'ZGEQRT2';
      infoc.INFOT = 1;
      zgeqrt2(-1, 0, A, 1, T, 1, INFO );
      chkxer('ZGEQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgeqrt2(0, -1, A, 1, T, 1, INFO );
      chkxer('ZGEQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgeqrt2(2, 1, A, 1, T, 1, INFO );
      chkxer('ZGEQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      zgeqrt2(2, 2, A, 2, T, 1, INFO );
      chkxer('ZGEQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZGEQRT3

     srnamc.SRNAMT = 'ZGEQRT3';
      infoc.INFOT = 1;
      zgeqrt3(-1, 0, A, 1, T, 1, INFO );
      chkxer('ZGEQRT3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgeqrt3(0, -1, A, 1, T, 1, INFO );
      chkxer('ZGEQRT3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgeqrt3(2, 1, A, 1, T, 1, INFO );
      chkxer('ZGEQRT3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      zgeqrt3(2, 2, A, 2, T, 1, INFO );
      chkxer('ZGEQRT3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZGEMQRT

     srnamc.SRNAMT = 'ZGEMQRT';
      infoc.INFOT = 1;
      zgemqrt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      zgemqrt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      zgemqrt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      zgemqrt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zgemqrt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      zgemqrt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      zgemqrt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      zgemqrt('R', 'N', 1, 2, 1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      zgemqrt('L', 'N', 2, 1, 1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      zgemqrt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO );
      chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      zgemqrt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO );
      chkxer('ZGEMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      }
