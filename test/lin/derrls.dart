import 'common.dart';

      void derrls(PATH, NUNIT ) {

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
      String             C2;
      int                INFO, IRNK;
      double             RCOND;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double             A( NMAX, NMAX ), B( NMAX, NMAX ), S( NMAX ), W( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGELS, DGELSD, DGELSS, DGELST, DGELSY, DGETSLS
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
      // .. Executable Statements ..

      infoc.NOUT = NUNIT;
      WRITE( infoc.NOUT, FMT = * );
      C2 = PATH( 2: 3 );
      A[1, 1] = 1.0;
      A[1, 2] = 2.0;
      A[2, 2] = 3.0;
      A[2, 1] = 4.0;
      infoc.OK = true;

      if ( lsamen( 2, C2, 'LS' ) ) {

         // Test error exits for the least squares driver routines.

         // DGELS

         srnamc.SRNAMT = 'DGELS ';
         infoc.INFOT = 1;
         dgels('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELS ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgels('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELS ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgels('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELS ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dgels('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELS ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dgels('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('DGELS ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dgels('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('DGELS ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dgels('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('DGELS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dgels('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELS ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DGELST

         srnamc.SRNAMT = 'DGELST';
         infoc.INFOT = 1;
         dgelst('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELST', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgelst('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELST', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgelst('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELST', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dgelst('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELST', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dgelst('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('DGELST', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dgelst('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('DGELST', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dgelst('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('DGELST', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dgelst('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELST', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DGETSLS

         srnamc.SRNAMT = 'DGETSLS';
         infoc.INFOT = 1;
         dgetsls('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGETSLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgetsls('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGETSLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgetsls('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGETSLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dgetsls('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('DGETSLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dgetsls('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('DGETSLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dgetsls('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('DGETSLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dgetsls('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('DGETSLS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DGELSS

         srnamc.SRNAMT = 'DGELSS';
         infoc.INFOT = 1;
         dgelss(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO );
         chkxer('DGELSS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgelss(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO );
         chkxer('DGELSS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgelss(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO );
         chkxer('DGELSS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dgelss(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, INFO );
         chkxer('DGELSS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dgelss(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, INFO );
         chkxer('DGELSS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DGELSY

         srnamc.SRNAMT = 'DGELSY';
         infoc.INFOT = 1;
         dgelsy(-1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('DGELSY', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgelsy(0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('DGELSY', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgelsy(0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('DGELSY', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dgelsy(2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('DGELSY', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dgelsy(2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('DGELSY', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         dgelsy(2, 2, 1, A, 2, B, 2, IP, RCOND, IRNK, W, 1, INFO );
         chkxer('DGELSY', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DGELSD

         srnamc.SRNAMT = 'DGELSD';
         infoc.INFOT = 1;
         dgelsd(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('DGELSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgelsd(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('DGELSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgelsd(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('DGELSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         dgelsd(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('DGELSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         dgelsd(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('DGELSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         dgelsd(2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, IP, INFO );
         chkxer('DGELSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      return;
      }
