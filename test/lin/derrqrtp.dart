import 'common.dart';
      void derrqrtp(PATH, final int NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      double             A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DTPQRT2, DTPQRT, DTPMQRT
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

      infoc.NOUT = NUNIT;
      WRITE( infoc.NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I][J] = 1.0 / (I+J).toDouble();
            C[I][J] = 1.0 / (I+J).toDouble();
            T[I][J] = 1.0 / (I+J).toDouble();
         }
         W[J] = 0.0;
      }
      infoc.OK = true;

      // Error exits for TPQRT factorization

      // DTPQRT

      srnamc.SRNAMT = 'DTPQRT';
      infoc.INFOT = 1;
      dtpqrt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dtpqrt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dtpqrt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dtpqrt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dtpqrt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dtpqrt(0, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      dtpqrt(1, 2, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dtpqrt(2, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('DTPQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dtpqrt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO );
      chkxer('DTPQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DTPQRT2

      srnamc.SRNAMT = 'DTPQRT2';
      infoc.INFOT = 1;
      dtpqrt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('DTPQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dtpqrt2(0, -1, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('DTPQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dtpqrt2(0, 0, -1, A, 1, B, 1, T, 1, INFO );
      chkxer('DTPQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dtpqrt2(2, 2, 0, A, 1, B, 2, T, 2, INFO );
      chkxer('DTPQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dtpqrt2(2, 2, 0, A, 2, B, 1, T, 2, INFO );
      chkxer('DTPQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      dtpqrt2(2, 2, 0, A, 2, B, 2, T, 1, INFO );
      chkxer('DTPQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DTPMQRT

      srnamc.SRNAMT = 'DTPMQRT';
      infoc.INFOT = 1;
      dtpmqrt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dtpmqrt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dtpmqrt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dtpmqrt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dtpmqrt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      infoc.INFOT = 6;
      dtpmqrt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dtpmqrt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      dtpmqrt('R', 'N', 1, 2, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 9;
      dtpmqrt('L', 'N', 2, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 11;
      dtpmqrt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO );
      chkxer('DTPMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 13;
      dtpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO );
      chkxer('DTPMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 15;
      dtpmqrt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO );
      chkxer('DTPMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      }
