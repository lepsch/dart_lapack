import 'common.dart';

      void derrlq(PATH, final int NUNIT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      double             A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGELQ2, DGELQF, DORGL2, DORGLQ, DORML2, DORMLQ
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

      // Error exits for LQ factorization

      // DGELQF

      srnamc.SRNAMT = 'DGELQF';
      infoc.INFOT = 1;
      dgelqf(-1, 0, A, 1, B, W, 1, INFO );
      chkxer('DGELQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgelqf(0, -1, A, 1, B, W, 1, INFO );
      chkxer('DGELQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgelqf(2, 1, A, 1, B, W, 2, INFO );
      chkxer('DGELQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dgelqf(2, 1, A, 2, B, W, 1, INFO );
      chkxer('DGELQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DGELQ2

      srnamc.SRNAMT = 'DGELQ2';
      infoc.INFOT = 1;
      dgelq2(-1, 0, A, 1, B, W, INFO );
      chkxer('DGELQ2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dgelq2(0, -1, A, 1, B, W, INFO );
      chkxer('DGELQ2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dgelq2(2, 1, A, 1, B, W, INFO );
      chkxer('DGELQ2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DORGLQ

      srnamc.SRNAMT = 'DORGLQ';
      infoc.INFOT = 1;
      dorglq(-1, 0, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorglq(0, -1, 0, A, 1, X, W, 1, INFO );
      chkxer('DORGLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorglq(2, 1, 0, A, 2, X, W, 2, INFO );
      chkxer('DORGLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorglq(0, 0, -1, A, 1, X, W, 1, INFO );
      chkxer('DORGLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorglq(1, 1, 2, A, 1, X, W, 1, INFO );
      chkxer('DORGLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dorglq(2, 2, 0, A, 1, X, W, 2, INFO );
      chkxer('DORGLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 8;
      dorglq(2, 2, 0, A, 2, X, W, 1, INFO );
      chkxer('DORGLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DORGL2

      srnamc.SRNAMT = 'DORGL2';
      infoc.INFOT = 1;
      dorgl2(-1, 0, 0, A, 1, X, W, INFO );
      chkxer('DORGL2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorgl2(0, -1, 0, A, 1, X, W, INFO );
      chkxer('DORGL2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorgl2(2, 1, 0, A, 2, X, W, INFO );
      chkxer('DORGL2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorgl2(0, 0, -1, A, 1, X, W, INFO );
      chkxer('DORGL2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorgl2(1, 1, 2, A, 1, X, W, INFO );
      chkxer('DORGL2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dorgl2(2, 2, 0, A, 1, X, W, INFO );
      chkxer('DORGL2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DORMLQ

      srnamc.SRNAMT = 'DORMLQ';
      infoc.INFOT = 1;
      dormlq('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dormlq('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dormlq('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dormlq('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dormlq('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dormlq('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dormlq('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dormlq('L', 'N', 2, 0, 2, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dormlq('R', 'N', 0, 2, 2, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dormlq('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      dormlq('L', 'N', 1, 2, 0, A, 1, X, AF, 1, W, 1, INFO );
      chkxer('DORMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 12;
      dormlq('R', 'N', 2, 1, 0, A, 1, X, AF, 2, W, 1, INFO );
      chkxer('DORMLQ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // DORML2

      srnamc.SRNAMT = 'DORML2';
      infoc.INFOT = 1;
      dorml2('/', 'N', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 2;
      dorml2('L', '/', 0, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 3;
      dorml2('L', 'N', -1, 0, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 4;
      dorml2('L', 'N', 0, -1, 0, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dorml2('L', 'N', 0, 0, -1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dorml2('L', 'N', 0, 1, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 5;
      dorml2('R', 'N', 1, 0, 1, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dorml2('L', 'N', 2, 1, 2, A, 1, X, AF, 2, W, INFO );
      chkxer('DORML2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 7;
      dorml2('R', 'N', 1, 2, 2, A, 1, X, AF, 1, W, INFO );
      chkxer('DORML2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      infoc.INFOT = 10;
      dorml2('L', 'N', 2, 1, 0, A, 2, X, AF, 1, W, INFO );
      chkxer('DORML2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      }
