import 'common.dart';
      void derrqrt(final int PATH, final int NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      double             A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGEQRT2, DGEQRT3, DGEQRT, DGEMQRT
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

      // Error exits for QRT factorization

      // DGEQRT

      srnamc.SRNAMT = 'DGEQRT';
      infoc.INFOT = 1;
      dgeqrt(-1, 0, 1, A, 1, T, 1, W, INFO );
      chkxer('DGEQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgeqrt(0, -1, 1, A, 1, T, 1, W, INFO );
      chkxer('DGEQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dgeqrt(0, 0, 0, A, 1, T, 1, W, INFO );
      chkxer('DGEQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dgeqrt(2, 1, 1, A, 1, T, 1, W, INFO );
      chkxer('DGEQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dgeqrt(2, 2, 2, A, 2, T, 1, W, INFO );
      chkxer('DGEQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGEQRT2

      srnamc.SRNAMT = 'DGEQRT2';
      infoc.INFOT = 1;
      dgeqrt2(-1, 0, A, 1, T, 1, INFO );
      chkxer('DGEQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgeqrt2(0, -1, A, 1, T, 1, INFO );
      chkxer('DGEQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgeqrt2(2, 1, A, 1, T, 1, INFO );
      chkxer('DGEQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      dgeqrt2(2, 2, A, 2, T, 1, INFO );
      chkxer('DGEQRT2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGEQRT3

      srnamc.SRNAMT = 'DGEQRT3';
      infoc.INFOT = 1;
      dgeqrt3(-1, 0, A, 1, T, 1, INFO );
      chkxer('DGEQRT3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgeqrt3(0, -1, A, 1, T, 1, INFO );
      chkxer('DGEQRT3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgeqrt3(2, 1, A, 1, T, 1, INFO );
      chkxer('DGEQRT3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      dgeqrt3(2, 2, A, 2, T, 1, INFO );
      chkxer('DGEQRT3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGEMQRT

      srnamc.SRNAMT = 'DGEMQRT';
      infoc.INFOT = 1;
      dgemqrt('/', 'N', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgemqrt('L', '/', 0, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dgemqrt('L', 'N', -1, 0, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgemqrt('L', 'N', 0, -1, 0, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dgemqrt('L', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dgemqrt('R', 'N', 0, 0, -1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 6;
      dgemqrt('L', 'N', 0, 0, 0, 0, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dgemqrt('R', 'N', 1, 2, 1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dgemqrt('L', 'N', 2, 1, 1, 1, A, 1, T, 1, C, 1, W, INFO );
      chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dgemqrt('R', 'N', 1, 1, 1, 1, A, 1, T, 0, C, 1, W, INFO );
      chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      dgemqrt('L', 'N', 1, 1, 1, 1, A, 1, T, 1, C, 0, W, INFO );
      chkxer('DGEMQRT', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      }
