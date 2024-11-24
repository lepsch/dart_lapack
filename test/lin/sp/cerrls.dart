// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cerrls(final int PATH, final int NUNIT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      String             C2;
      int                INFO, IRNK;
      double               RCOND;
      int                IP( NMAX );
      double               RW( NMAX ), S( NMAX );
      Complex            A( NMAX, NMAX ), B( NMAX, NMAX ), W( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CGELS, CGELSD, CGELSS, CGELST, CGELSY, CGETSLS
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT

      NOUT = NUNIT;
      C2 = PATH( 2: 3 );
      A[1][1] = ( 1.0, 0.0 );
      A[1][2] = ( 2.0, 0.0 );
      A[2][2] = ( 3.0, 0.0 );
      A[2][1] = ( 4.0, 0.0 );
      OK = true;
      WRITE( NOUT, FMT = * );

      // Test error exits for the least squares driver routines.

      if ( lsamen( 2, C2, 'LS' ) ) {

         // CGELS

        srnamc.SRNAMT = 'CGELS ';
         INFOT = 1;
         cgels('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('CGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgels('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('CGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgels('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('CGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgels('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('CGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgels('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('CGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgels('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('CGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgels('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('CGELS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cgels('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('CGELS ', INFOT, NOUT, LERR, OK );

         // CGELST

        srnamc.SRNAMT = 'CGELST';
         INFOT = 1;
         cgelst('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('CGELST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgelst('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('CGELST', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgelst('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('CGELST', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgelst('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('CGELST', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgelst('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('CGELST', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgelst('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('CGELST', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgelst('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('CGELST', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cgelst('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('CGELST', INFOT, NOUT, LERR, OK );

         // CGETSLS

        srnamc.SRNAMT = 'CGETSLS';
         INFOT = 1;
         cgetsls('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('CGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgetsls('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('CGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgetsls('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('CGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgetsls('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('CGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgetsls('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('CGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgetsls('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('CGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgetsls('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('CGETSLS', INFOT, NOUT, LERR, OK );

         // CGELSS

        srnamc.SRNAMT = 'CGELSS';
         INFOT = 1;
         cgelss(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO );
         chkxer('CGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgelss(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO );
         chkxer('CGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgelss(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO );
         chkxer('CGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgelss(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, RW, INFO );
         chkxer('CGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgelss(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, RW, INFO );
         chkxer('CGELSS', INFOT, NOUT, LERR, OK );

         // CGELSY

        srnamc.SRNAMT = 'CGELSY';
         INFOT = 1;
         cgelsy(-1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('CGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgelsy(0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('CGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgelsy(0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('CGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgelsy(2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('CGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgelsy(2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('CGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cgelsy(0, 3, 0, A, 1, B, 3, IP, RCOND, IRNK, W, 1, RW, INFO );
         chkxer('CGELSY', INFOT, NOUT, LERR, OK );

         // CGELSD

        srnamc.SRNAMT = 'CGELSD';
         INFOT = 1;
         cgelsd(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('CGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgelsd(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('CGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgelsd(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('CGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgelsd(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('CGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgelsd(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('CGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cgelsd(2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, RW, IP, INFO );
         chkxer('CGELSD', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      }
