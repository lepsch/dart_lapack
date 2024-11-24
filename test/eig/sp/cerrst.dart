// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cerrst(final int PATH, final int NUNIT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX, LIW, LW;
      const              NMAX = 3, LIW = 12*NMAX, LW = 20*NMAX ;
      String             C2;
      int                I, INFO, J, M, N, NT;
      int                I1( NMAX ), I2( NMAX ), I3( NMAX ), IW( LIW );
      double               D( NMAX ), E( NMAX ), R( LW ), RW( LW ), X( NMAX )       Complex            A( NMAX, NMAX ), C( NMAX, NMAX ), Q( NMAX, NMAX ), TAU( NMAX ), W( LW ), Z( NMAX, NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHBEV, CHBEVD, CHBEVX, CHBTRD, CHEEV, CHEEVD, CHEEVR, CHEEVX, CHETRD, CHKXER, CHPEV, CHPEVD, CHPEVX, CHPTRD, CPTEQR, CSTEDC, CSTEIN, CSTEQR, CUNGTR, CUNMTR, CUPGTR, CUPMTR, CHETD2, CHEEVD_2STAGE, CHEEVR_2STAGE, CHEEVX_2STAGE, CHEEV_2STAGE, CHBEV_2STAGE, CHBEVD_2STAGE, CHBEVX_2STAGE, CHETRD_2STAGE, CHETRD_HE2HB, CHETRD_HB2ST
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
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = 1. / REAL( I+J );
         } // 10
      } // 20
      for (J = 1; J <= NMAX; J++) { // 30
         D[J] = double( J );
         E[J] = 0.0;
         I1[J] = J;
         I2[J] = J;
         TAU[J] = 1.;
      } // 30
      OK = true;
      NT = 0;

      // Test error exits for the ST path.

      if ( lsamen( 2, C2, 'ST' ) ) {

         // CHETRD

        srnamc.SRNAMT = 'CHETRD';
         INFOT = 1;
         chetrd('/', 0, A, 1, D, E, TAU, W, 1, INFO );
         chkxer('CHETRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrd('U', -1, A, 1, D, E, TAU, W, 1, INFO );
         chkxer('CHETRD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetrd('U', 2, A, 1, D, E, TAU, W, 1, INFO );
         chkxer('CHETRD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chetrd('U', 0, A, 1, D, E, TAU, W, 0, INFO );
         chkxer('CHETRD', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // CHETD2

        srnamc.SRNAMT = 'CHETD2';
         INFOT = 1;
         chetd2('/', 0, A, 1, D, E, TAU, INFO );
         chkxer('CHETD2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetd2('U', -1, A, 1, D, E, TAU, INFO );
         chkxer('CHETD2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetd2('U', 2, A, 1, D, E, TAU, INFO );
         chkxer('CHETD2', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // CHETRD_2STAGE

        srnamc.SRNAMT = 'CHETRD_2STAGE';
         INFOT = 1;
         chetrd_2stage('/', 'U', 0, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('CHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         chetrd_2stage('H', 'U', 0, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('CHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrd_2stage('N', '/', 0, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('CHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chetrd_2stage('N', 'U', -1, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('CHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chetrd_2stage('N', 'U', 2, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('CHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         chetrd_2stage('N', 'U', 0, A, 1, D, E, TAU,  C, 0, W, 1, INFO );
         chkxer('CHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         chetrd_2stage('N', 'U', 0, A, 1, D, E, TAU,  C, 1, W, 0, INFO );
         chkxer('CHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // CHETRD_HE2HB

        srnamc.SRNAMT = 'CHETRD_HE2HB';
         INFOT = 1;
         chetrd_he2hb('/', 0, 0, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('CHETRD_HE2HB', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrd_he2hb('U', -1, 0, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('CHETRD_HE2HB', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chetrd_he2hb('U', 0, -1, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('CHETRD_HE2HB', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chetrd_he2hb('U', 2, 0, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('CHETRD_HE2HB', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chetrd_he2hb('U', 0, 2, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('CHETRD_HE2HB', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         chetrd_he2hb('U', 0, 0, A, 1, C, 1, TAU, W, 0, INFO );
         chkxer('CHETRD_HE2HB', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // CHETRD_HB2ST

        srnamc.SRNAMT = 'CHETRD_HB2ST';
         INFOT = 1;
         chetrd_hb2st('/', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrd_hb2st('Y', '/', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrd_hb2st('Y', 'H', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chetrd_hb2st('Y', 'N', '/', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetrd_hb2st('Y', 'N', 'U', -1, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chetrd_hb2st('Y', 'N', 'U', 0, -1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chetrd_hb2st('Y', 'N', 'U', 0, 1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chetrd_hb2st('Y', 'N', 'U', 0, 0, A, 1, D, E,  C, 0, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chetrd_hb2st('Y', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 0, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // CUNGTR

        srnamc.SRNAMT = 'CUNGTR';
         INFOT = 1;
         cungtr('/', 0, A, 1, TAU, W, 1, INFO );
         chkxer('CUNGTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cungtr('U', -1, A, 1, TAU, W, 1, INFO );
         chkxer('CUNGTR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cungtr('U', 2, A, 1, TAU, W, 1, INFO );
         chkxer('CUNGTR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cungtr('U', 3, A, 3, TAU, W, 1, INFO );
         chkxer('CUNGTR', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // CUNMTR

        srnamc.SRNAMT = 'CUNMTR';
         INFOT = 1;
         cunmtr('/', 'U', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cunmtr('L', '/', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cunmtr('L', 'U', '/', 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cunmtr('L', 'U', 'N', -1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cunmtr('L', 'U', 'N', 0, -1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cunmtr('L', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('CUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cunmtr('R', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cunmtr('L', 'U', 'N', 2, 0, A, 2, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cunmtr('L', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('CUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cunmtr('R', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('CUNMTR', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // CHPTRD

        srnamc.SRNAMT = 'CHPTRD';
         INFOT = 1;
         chptrd('/', 0, A, D, E, TAU, INFO );
         chkxer('CHPTRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chptrd('U', -1, A, D, E, TAU, INFO );
         chkxer('CHPTRD', INFOT, NOUT, LERR, OK );
         NT = NT + 2;

         // CUPGTR

        srnamc.SRNAMT = 'CUPGTR';
         INFOT = 1;
         cupgtr('/', 0, A, TAU, Z, 1, W, INFO );
         chkxer('CUPGTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cupgtr('U', -1, A, TAU, Z, 1, W, INFO );
         chkxer('CUPGTR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cupgtr('U', 2, A, TAU, Z, 1, W, INFO );
         chkxer('CUPGTR', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // CUPMTR

        srnamc.SRNAMT = 'CUPMTR';
         INFOT = 1;
         cupmtr('/', 'U', 'N', 0, 0, A, TAU, C, 1, W, INFO );
         chkxer('CUPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cupmtr('L', '/', 'N', 0, 0, A, TAU, C, 1, W, INFO );
         chkxer('CUPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cupmtr('L', 'U', '/', 0, 0, A, TAU, C, 1, W, INFO );
         chkxer('CUPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cupmtr('L', 'U', 'N', -1, 0, A, TAU, C, 1, W, INFO );
         chkxer('CUPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cupmtr('L', 'U', 'N', 0, -1, A, TAU, C, 1, W, INFO );
         chkxer('CUPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cupmtr('L', 'U', 'N', 2, 0, A, TAU, C, 1, W, INFO );
         chkxer('CUPMTR', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // CPTEQR

        srnamc.SRNAMT = 'CPTEQR';
         INFOT = 1;
         cpteqr('/', 0, D, E, Z, 1, RW, INFO );
         chkxer('CPTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cpteqr('N', -1, D, E, Z, 1, RW, INFO );
         chkxer('CPTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cpteqr('V', 2, D, E, Z, 1, RW, INFO );
         chkxer('CPTEQR', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // CSTEIN

        srnamc.SRNAMT = 'CSTEIN';
         INFOT = 1;
         cstein(-1, D, E, 0, X, I1, I2, Z, 1, RW, IW, I3, INFO );
         chkxer('CSTEIN', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cstein(0, D, E, -1, X, I1, I2, Z, 1, RW, IW, I3, INFO );
         chkxer('CSTEIN', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cstein(0, D, E, 1, X, I1, I2, Z, 1, RW, IW, I3, INFO );
         chkxer('CSTEIN', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cstein(2, D, E, 0, X, I1, I2, Z, 1, RW, IW, I3, INFO );
         chkxer('CSTEIN', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // CSTEQR

        srnamc.SRNAMT = 'CSTEQR';
         INFOT = 1;
         csteqr('/', 0, D, E, Z, 1, RW, INFO );
         chkxer('CSTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csteqr('N', -1, D, E, Z, 1, RW, INFO );
         chkxer('CSTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         csteqr('V', 2, D, E, Z, 1, RW, INFO );
         chkxer('CSTEQR', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // CSTEDC

        srnamc.SRNAMT = 'CSTEDC';
         INFOT = 1;
         cstedc('/', 0, D, E, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cstedc('N', -1, D, E, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cstedc('V', 2, D, E, Z, 1, W, 4, RW, 23, IW, 28, INFO );
         chkxer('CSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cstedc('N', 2, D, E, Z, 1, W, 0, RW, 1, IW, 1, INFO );
         chkxer('CSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cstedc('V', 2, D, E, Z, 2, W, 0, RW, 23, IW, 28, INFO );
         chkxer('CSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cstedc('N', 2, D, E, Z, 1, W, 1, RW, 0, IW, 1, INFO );
         chkxer('CSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cstedc('I', 2, D, E, Z, 2, W, 1, RW, 1, IW, 12, INFO );
         chkxer('CSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cstedc('V', 2, D, E, Z, 2, W, 4, RW, 1, IW, 28, INFO );
         chkxer('CSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cstedc('N', 2, D, E, Z, 1, W, 1, RW, 1, IW, 0, INFO );
         chkxer('CSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cstedc('I', 2, D, E, Z, 2, W, 1, RW, 23, IW, 0, INFO );
         chkxer('CSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cstedc('V', 2, D, E, Z, 2, W, 4, RW, 23, IW, 0, INFO );
         chkxer('CSTEDC', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // CHEEVD

        srnamc.SRNAMT = 'CHEEVD';
         INFOT = 1;
         cheevd('/', 'U', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cheevd('N', '/', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cheevd('N', 'U', -1, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cheevd('N', 'U', 2, A, 1, X, W, 3, RW, 2, IW, 1, INFO );
         chkxer('CHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cheevd('N', 'U', 1, A, 1, X, W, 0, RW, 1, IW, 1, INFO );
         chkxer('CHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cheevd('N', 'U', 2, A, 2, X, W, 2, RW, 2, IW, 1, INFO );
         chkxer('CHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cheevd('V', 'U', 2, A, 2, X, W, 3, RW, 25, IW, 12, INFO );
         chkxer('CHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cheevd('N', 'U', 1, A, 1, X, W, 1, RW, 0, IW, 1, INFO );
         chkxer('CHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cheevd('N', 'U', 2, A, 2, X, W, 3, RW, 1, IW, 1, INFO );
         chkxer('CHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cheevd('V', 'U', 2, A, 2, X, W, 8, RW, 18, IW, 12, INFO );
         chkxer('CHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cheevd('N', 'U', 1, A, 1, X, W, 1, RW, 1, IW, 0, INFO );
         chkxer('CHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cheevd('V', 'U', 2, A, 2, X, W, 8, RW, 25, IW, 11, INFO );
         chkxer('CHEEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // CHEEVD_2STAGE

        srnamc.SRNAMT = 'CHEEVD_2STAGE';
         INFOT = 1;
         cheevd_2stage('/', 'U', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         cheevd_2stage('V', 'U', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cheevd_2stage('N', '/', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cheevd_2stage('N', 'U', -1, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cheevd_2stage('N', 'U', 2, A, 1, X, W, 3, RW, 2, IW, 1, INFO );
         chkxer('CHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cheevd_2stage('N', 'U', 1, A, 1, X, W, 0, RW, 1, IW, 1, INFO );
         chkxer('CHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cheevd_2stage('N', 'U', 2, A, 2, X, W, 2, RW, 2, IW, 1, INFO );
         chkxer('CHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 8
          // CALL CHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 3,
      // $                            RW, 25, IW, 12, INFO )
      //     CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 10;
         cheevd_2stage('N', 'U', 1, A, 1, X, W, 1, RW, 0, IW, 1, INFO );
         chkxer('CHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cheevd_2stage('N', 'U', 2, A, 2, X, W, 25, RW, 1, IW, 1, INFO );
         chkxer('CHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 10
          // CALL CHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 8,
      // $                            RW, 18, IW, 12, INFO )
      //     CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 12;
         cheevd_2stage('N', 'U', 1, A, 1, X, W, 1, RW, 1, IW, 0, INFO );
         chkxer('CHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
          // CALL CHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 8,
      // $                            RW, 25, IW, 11, INFO )
      //     CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
         NT = NT + 10;

         // CHEEV

        srnamc.SRNAMT = 'CHEEV ';
         INFOT = 1;
         cheev('/', 'U', 0, A, 1, X, W, 1, RW, INFO );
         chkxer('CHEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cheev('N', '/', 0, A, 1, X, W, 1, RW, INFO );
         chkxer('CHEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cheev('N', 'U', -1, A, 1, X, W, 1, RW, INFO );
         chkxer('CHEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cheev('N', 'U', 2, A, 1, X, W, 3, RW, INFO );
         chkxer('CHEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cheev('N', 'U', 2, A, 2, X, W, 2, RW, INFO );
         chkxer('CHEEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 5;

         // CHEEV_2STAGE

        srnamc.SRNAMT = 'CHEEV_2STAGE ';
         INFOT = 1;
         cheev_2stage('/', 'U', 0, A, 1, X, W, 1, RW, INFO );
         chkxer('CHEEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         cheev_2stage('V', 'U', 0, A, 1, X, W, 1, RW, INFO );
         chkxer('CHEEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cheev_2stage('N', '/', 0, A, 1, X, W, 1, RW, INFO );
         chkxer('CHEEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cheev_2stage('N', 'U', -1, A, 1, X, W, 1, RW, INFO );
         chkxer('CHEEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cheev_2stage('N', 'U', 2, A, 1, X, W, 3, RW, INFO );
         chkxer('CHEEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cheev_2stage('N', 'U', 2, A, 2, X, W, 2, RW, INFO );
         chkxer('CHEEV_2STAGE ', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // CHEEVX

        srnamc.SRNAMT = 'CHEEVX';
         INFOT = 1;
         cheevx('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('CHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cheevx('V', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('CHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cheevx('V', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         INFOT = 4;
         cheevx('V', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('CHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cheevx('V', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, 3, RW, IW, I3, INFO );
         chkxer('CHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cheevx('V', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('CHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cheevx('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('CHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cheevx('V', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 2, W, 3, RW, IW, I3, INFO );
         chkxer('CHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cheevx('V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 3, RW, IW, I3, INFO );
         chkxer('CHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         cheevx('V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, 2, RW, IW, I1, INFO );
         chkxer('CHEEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // CHEEVX_2STAGE

        srnamc.SRNAMT = 'CHEEVX_2STAGE';
         INFOT = 1;
         cheevx_2stage('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('CHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         cheevx_2stage('V', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('CHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cheevx_2stage('N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('CHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cheevx_2stage('N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         INFOT = 4;
         cheevx_2stage('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('CHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cheevx_2stage('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, 3, RW, IW, I3, INFO );
         chkxer('CHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cheevx_2stage('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('CHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cheevx_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('CHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cheevx_2stage('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 2, W, 3, RW, IW, I3, INFO );
         chkxer('CHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cheevx_2stage('N', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 0, W, 3, RW, IW, I3, INFO );
         chkxer('CHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         cheevx_2stage('N', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, 0, RW, IW, I1, INFO );
         chkxer('CHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // CHEEVR

        srnamc.SRNAMT = 'CHEEVR';
         N = 1;
         INFOT = 1;
         cheevr('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cheevr('V', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cheevr('V', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cheevr('V', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cheevr('V', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cheevr('V', 'V', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cheevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 10;

         cheevr('V', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cheevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 0, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         cheevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 0, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         cheevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 0, IW( 2*N-1 ), 10*N, INFO );
         chkxer('CHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         cheevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW, 0, INFO );
         chkxer('CHEEVR', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // CHEEVR_2STAGE

        srnamc.SRNAMT = 'CHEEVR_2STAGE';
         N = 1;
         INFOT = 1;
         cheevr_2stage('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         cheevr_2stage('V', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cheevr_2stage('N', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cheevr_2stage('N', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cheevr_2stage('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cheevr_2stage('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cheevr_2stage('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cheevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cheevr_2stage('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cheevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 0, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         cheevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 0, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         cheevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, RW, 0, IW( 2*N-1 ), 10*N, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         cheevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, RW, 24*N, IW, 0, INFO );
         chkxer('CHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 13;

         // CHPEVD

        srnamc.SRNAMT = 'CHPEVD';
         INFOT = 1;
         chpevd('/', 'U', 0, A, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chpevd('N', '/', 0, A, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chpevd('N', 'U', -1, A, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chpevd('V', 'U', 2, A, X, Z, 1, W, 4, RW, 25, IW, 12, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chpevd('N', 'U', 1, A, X, Z, 1, W, 0, RW, 1, IW, 1, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chpevd('N', 'U', 2, A, X, Z, 2, W, 1, RW, 2, IW, 1, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chpevd('V', 'U', 2, A, X, Z, 2, W, 2, RW, 25, IW, 12, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chpevd('N', 'U', 1, A, X, Z, 1, W, 1, RW, 0, IW, 1, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chpevd('N', 'U', 2, A, X, Z, 2, W, 2, RW, 1, IW, 1, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chpevd('V', 'U', 2, A, X, Z, 2, W, 4, RW, 18, IW, 12, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chpevd('N', 'U', 1, A, X, Z, 1, W, 1, RW, 1, IW, 0, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chpevd('N', 'U', 2, A, X, Z, 2, W, 2, RW, 2, IW, 0, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chpevd('V', 'U', 2, A, X, Z, 2, W, 4, RW, 25, IW, 2, INFO );
         chkxer('CHPEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 13;

         // CHPEV

        srnamc.SRNAMT = 'CHPEV ';
         INFOT = 1;
         chpev('/', 'U', 0, A, X, Z, 1, W, RW, INFO );
         chkxer('CHPEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chpev('N', '/', 0, A, X, Z, 1, W, RW, INFO );
         chkxer('CHPEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chpev('N', 'U', -1, A, X, Z, 1, W, RW, INFO );
         chkxer('CHPEV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chpev('V', 'U', 2, A, X, Z, 1, W, RW, INFO );
         chkxer('CHPEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // CHPEVX

        srnamc.SRNAMT = 'CHPEVX';
         INFOT = 1;
         chpevx('/', 'A', 'U', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chpevx('V', '/', 'U', 0, A, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chpevx('V', 'A', '/', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chpevx('V', 'A', 'U', -1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chpevx('V', 'V', 'U', 1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chpevx('V', 'I', 'U', 1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chpevx('V', 'I', 'U', 2, A, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 2, W, RW, IW, I3, INFO );
         chkxer('CHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         chpevx('V', 'A', 'U', 2, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHPEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the HB path.

      } else if ( lsamen( 2, C2, 'HB' ) ) {

         // CHBTRD

        srnamc.SRNAMT = 'CHBTRD';
         INFOT = 1;
         chbtrd('/', 'U', 0, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('CHBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chbtrd('N', '/', 0, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('CHBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chbtrd('N', 'U', -1, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('CHBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chbtrd('N', 'U', 0, -1, A, 1, D, E, Z, 1, W, INFO );
         chkxer('CHBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         chbtrd('N', 'U', 1, 1, A, 1, D, E, Z, 1, W, INFO );
         chkxer('CHBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         chbtrd('V', 'U', 2, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('CHBTRD', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // CHETRD_HB2ST

        srnamc.SRNAMT = 'CHETRD_HB2ST';
         INFOT = 1;
         chetrd_hb2st('/', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrd_hb2st('N', '/', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrd_hb2st('N', 'H', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chetrd_hb2st('N', 'N', '/', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetrd_hb2st('N', 'N', 'U', -1, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chetrd_hb2st('N', 'N', 'U', 0, -1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chetrd_hb2st('N', 'N', 'U', 0, 1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chetrd_hb2st('N', 'N', 'U', 0, 0, A, 1, D, E,  C, 0, W, 1, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chetrd_hb2st('N', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 0, INFO );
         chkxer('CHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // CHBEVD

        srnamc.SRNAMT = 'CHBEVD';
         INFOT = 1;
         chbevd('/', 'U', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chbevd('N', '/', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chbevd('N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chbevd('N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         chbevd('N', 'U', 2, 1, A, 1, X, Z, 1, W, 2, RW, 2, IW, 1, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chbevd('V', 'U', 2, 1, A, 2, X, Z, 1, W, 8, RW, 25, IW, 12, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chbevd('N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chbevd('N', 'U', 2, 1, A, 2, X, Z, 2, W, 1, RW, 2, IW, 1, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chbevd('V', 'U', 2, 1, A, 2, X, Z, 2, W, 2, RW, 25, IW, 12, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chbevd('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 0, IW, 1, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chbevd('N', 'U', 2, 1, A, 2, X, Z, 2, W, 2, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chbevd('V', 'U', 2, 1, A, 2, X, Z, 2, W, 8, RW, 2, IW, 12, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         chbevd('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 0, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         chbevd('N', 'U', 2, 1, A, 2, X, Z, 2, W, 2, RW, 2, IW, 0, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         chbevd('V', 'U', 2, 1, A, 2, X, Z, 2, W, 8, RW, 25, IW, 2, INFO );
         chkxer('CHBEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 15;

         // CHBEVD_2STAGE

        srnamc.SRNAMT = 'CHBEVD_2STAGE';
         INFOT = 1;
         chbevd_2stage('/', 'U', 0, 0, A, 1, X, Z, 1,  W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         chbevd_2stage('V', 'U', 0, 0, A, 1, X, Z, 1,  W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chbevd_2stage('N', '/', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chbevd_2stage('N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chbevd_2stage('N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         chbevd_2stage('N', 'U', 2, 1, A, 1, X, Z, 1, W, 2, RW, 2, IW, 1, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chbevd_2stage('N', 'U', 2, 1, A, 2, X, Z, 0, W, 8, RW, 25, IW, 12, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chbevd_2stage('N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chbevd_2stage('N', 'U', 2, 1, A, 2, X, Z, 2, W, 1, RW, 2, IW, 1, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 11
          // CALL CHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
      // $                         W, 2, RW, 25, IW, 12, INFO )
      //     CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 13;
         chbevd_2stage('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 0, IW, 1, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chbevd_2stage('N', 'U', 2, 1, A, 2, X, Z, 2, W, 25, RW, 1, IW, 1, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 13
          // CALL CHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
      // $                          W, 25, RW, 2, IW, 12, INFO )
      //     CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 15;
         chbevd_2stage('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 0, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         chbevd_2stage('N', 'U', 2, 1, A, 2, X, Z, 2, W, 25, RW, 2, IW, 0, INFO );
         chkxer('CHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 15
          // CALL CHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
      // $                          W, 25, RW, 25, IW, 2, INFO )
      //     CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         NT = NT + 13;

         // CHBEV

        srnamc.SRNAMT = 'CHBEV ';
         INFOT = 1;
         chbev('/', 'U', 0, 0, A, 1, X, Z, 1, W, RW, INFO );
         chkxer('CHBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chbev('N', '/', 0, 0, A, 1, X, Z, 1, W, RW, INFO );
         chkxer('CHBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chbev('N', 'U', -1, 0, A, 1, X, Z, 1, W, RW, INFO );
         chkxer('CHBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chbev('N', 'U', 0, -1, A, 1, X, Z, 1, W, RW, INFO );
         chkxer('CHBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         chbev('N', 'U', 2, 1, A, 1, X, Z, 1, W, RW, INFO );
         chkxer('CHBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chbev('V', 'U', 2, 0, A, 1, X, Z, 1, W, RW, INFO );
         chkxer('CHBEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // CHBEV_2STAGE

        srnamc.SRNAMT = 'CHBEV_2STAGE ';
         INFOT = 1;
         chbev_2stage('/', 'U', 0, 0, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('CHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         chbev_2stage('V', 'U', 0, 0, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('CHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chbev_2stage('N', '/', 0, 0, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('CHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chbev_2stage('N', 'U', -1, 0, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('CHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chbev_2stage('N', 'U', 0, -1, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('CHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         chbev_2stage('N', 'U', 2, 1, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('CHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chbev_2stage('N', 'U', 2, 0, A, 1, X, Z, 0, W, 0, RW, INFO );
         chkxer('CHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chbev_2stage('N', 'U', 2, 0, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('CHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

         // CHBEVX

        srnamc.SRNAMT = 'CHBEVX';
         INFOT = 1;
         chbevx('/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chbevx('V', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chbevx('V', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         INFOT = 4;
         chbevx('V', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chbevx('V', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chbevx('V', 'A', 'U', 2, 1, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, RW, IW, I3, INFO );
         chkxer('CHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chbevx('V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, RW, IW, I3, INFO );
         chkxer('CHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         chbevx('V', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         chbevx('V', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chbevx('V', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         chbevx('V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('CHBEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // CHBEVX_2STAGE

        srnamc.SRNAMT = 'CHBEVX_2STAGE';
         INFOT = 1;
         chbevx_2stage('/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         INFOT = 1;
         chbevx_2stage('V', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('CHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chbevx_2stage('N', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('CHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chbevx_2stage('N', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         INFOT = 4;
         chbevx_2stage('N', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('CHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chbevx_2stage('N', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('CHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chbevx_2stage('N', 'A', 'U', 2, 1, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, 0, RW, IW, I3, INFO );
         chkxer('CHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 9
          // CALL CHBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 1,
      // $                       0.0, 0.0, 0, 0, 0.0,
      // $                       M, X, Z, 2, W, 0, RW, IW, I3, INFO )
      //     CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 11;
         chbevx_2stage('N', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('CHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         chbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('CHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         chbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('CHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         chbevx_2stage('N', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 0, W, 0, RW, IW, I3, INFO );
         chkxer('CHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         chbevx_2stage('N', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('CHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 12;
      }

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH, NT;
      } else {
         WRITE( NOUT, FMT = 9998 )PATH;
      }

 9999 FORMAT(' ${.a3} routines passed the tests of the error exits (${.i3} tests done)' );
 9998 FORMAT( ' *** ${.a3} routines failed the tests of the error exits ***' );

      }
