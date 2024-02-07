      void zerrqrtp(PATH, infoc.NUNIT ) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                infoc.NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      Complex         A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZTPQRT2, ZTPQRT, ZTPMQRT
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
            A[I][J] = 1.0 / DCMPLX((I+J).toDouble(),0.0);
            C[I][J] = 1.0 / DCMPLX((I+J).toDouble(),0.0);
            T[I][J] = 1.0 / DCMPLX((I+J).toDouble(),0.0);
         }
         W[J] = DCMPLX(0.0,0.0);
      }
      infoc.OK = true;

      // Error exits for TPQRT factorization

      // ZTPQRT

     srnamc.SRNAMT = 'ZTPQRT';
      infoc.INFOT = 1;
      ztpqrt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztpqrt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztpqrt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztpqrt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztpqrt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztpqrt(0, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      ztpqrt(1, 2, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      ztpqrt(2, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      ztpqrt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO );
      chkxer('ZTPQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZTPQRT2

     srnamc.SRNAMT = 'ZTPQRT2';
      infoc.INFOT = 1;
      ztpqrt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('ZTPQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztpqrt2(0, -1, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('ZTPQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztpqrt2(0, 0, -1, A, 1, B, 1, T, 1, INFO );
      chkxer('ZTPQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      ztpqrt2(2, 2, 0, A, 1, B, 2, T, 2, INFO );
      chkxer('ZTPQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      ztpqrt2(2, 2, 0, A, 2, B, 1, T, 2, INFO );
      chkxer('ZTPQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      ztpqrt2(2, 2, 0, A, 2, B, 2, T, 1, INFO );
      chkxer('ZTPQRT2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // ZTPMQRT

     srnamc.SRNAMT = 'ZTPMQRT';
      infoc.INFOT = 1;
      ztpmqrt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      ztpmqrt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      ztpmqrt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      ztpmqrt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      ztpmqrt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      infoc.INFOT = 6;
      ztpmqrt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      ztpmqrt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      ztpmqrt('R', 'N', 1, 2, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      ztpmqrt('L', 'N', 2, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 11;
      ztpmqrt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO );
      chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 13;
      ztpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO );
      chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 15;
      ztpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO );
      chkxer('ZTPMQRT', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      return;
      }
