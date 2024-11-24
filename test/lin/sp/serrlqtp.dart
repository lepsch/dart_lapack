// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void serrlqtp(final int PATH, final int NUNIT,) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      int                I, INFO, J;
      double               A( NMAX, NMAX ), T( NMAX, NMAX ), W( NMAX ), B( NMAX, NMAX ), C( NMAX, NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, STPLQT2, STPLQT, STPMLQT
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
      // INTRINSIC REAL

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) {
         for (I = 1; I <= NMAX; I++) {
            A[I][J] = 1. / REAL( I+J );
            C[I][J] = 1. / REAL( I+J );
            T[I][J] = 1. / REAL( I+J );
         }
         W[J] = 0.;
      }
      OK = true;

      // Error exits for TPLQT factorization

      // STPLQT

     srnamc.SRNAMT = 'STPLQT';
      INFOT = 1;
      stplqt(-1, 1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      stplqt(1, -1, 0, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stplqt(0, 1, -1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stplqt(0, 1, 1, 1, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      stplqt(0, 1, 0, 0, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      stplqt(1, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 6;
      stplqt(2, 1, 0, 2, A, 1, B, 1, T, 1, W, INFO );
      chkxer('STPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 8;
      stplqt(2, 1, 0, 1, A, 2, B, 1, T, 1, W, INFO );
      chkxer('STPLQT', INFOT, NOUT, LERR, OK );
      INFOT = 10;
      stplqt(2, 2, 1, 2, A, 2, B, 2, T, 1, W, INFO );
      chkxer('STPLQT', INFOT, NOUT, LERR, OK );

      // STPLQT2

     srnamc.SRNAMT = 'STPLQT2';
      INFOT = 1;
      stplqt2(-1, 0, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('STPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      stplqt2(0, -1, 0, A, 1, B, 1, T, 1, INFO );
      chkxer('STPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stplqt2(0, 0, -1, A, 1, B, 1, T, 1, INFO );
      chkxer('STPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      stplqt2(2, 2, 0, A, 1, B, 2, T, 2, INFO );
      chkxer('STPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      stplqt2(2, 2, 0, A, 2, B, 1, T, 2, INFO );
      chkxer('STPLQT2', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      stplqt2(2, 2, 0, A, 2, B, 2, T, 1, INFO );
      chkxer('STPLQT2', INFOT, NOUT, LERR, OK );

      // STPMLQT

     srnamc.SRNAMT = 'STPMLQT';
      INFOT = 1;
      stpmlqt('/', 'N', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 2;
      stpmlqt('L', '/', 0, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 3;
      stpmlqt('L', 'N', -1, 0, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 4;
      stpmlqt('L', 'N', 0, -1, 0, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 5;
      stpmlqt('L', 'N', 0, 0, -1, 0, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      INFOT = 6;
      stpmlqt('L', 'N', 0, 0, 0, -1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 7;
      stpmlqt('L', 'N', 0, 0, 0, 0, 0, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 9;
      stpmlqt('R', 'N', 2, 2, 2, 1, 1, A, 1, T, 1, B, 1, C, 1, W, INFO );
      chkxer('STPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 11;
      stpmlqt('R', 'N', 1, 1, 1, 1, 1, A, 1, T, 0, B, 1, C, 1, W, INFO );
      chkxer('STPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 13;
      stpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 0, C, 1, W, INFO );
      chkxer('STPMLQT', INFOT, NOUT, LERR, OK );
      INFOT = 15;
      stpmlqt('L', 'N', 1, 1, 1, 1, 1, A, 1, T, 1, B, 1, C, 0, W, INFO );
      chkxer('STPMLQT', INFOT, NOUT, LERR, OK );

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      }
