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
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGELS, DGELSD, DGELSS, DGELST, DGELSY, DGETSLS
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );
      A( 1, 1 ) = 1.0;
      A( 1, 2 ) = 2.0;
      A( 2, 2 ) = 3.0;
      A( 2, 1 ) = 4.0;
      OK = true;

      if ( LSAMEN( 2, C2, 'LS' ) ) {

         // Test error exits for the least squares driver routines.

         // DGELS

         SRNAMT = 'DGELS ';
         INFOT = 1;
         dgels('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgels('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgels('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgels('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dgels('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('DGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dgels('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('DGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dgels('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('DGELS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dgels('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELS ', INFOT, NOUT, LERR, OK );

         // DGELST

         SRNAMT = 'DGELST';
         INFOT = 1;
         dgelst('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgelst('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELST', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgelst('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELST', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgelst('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELST', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dgelst('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('DGELST', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dgelst('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('DGELST', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dgelst('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('DGELST', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dgelst('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGELST', INFOT, NOUT, LERR, OK );

         // DGETSLS

         SRNAMT = 'DGETSLS';
         INFOT = 1;
         dgetsls('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgetsls('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgetsls('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('DGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgetsls('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('DGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dgetsls('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('DGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dgetsls('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('DGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dgetsls('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('DGETSLS', INFOT, NOUT, LERR, OK );

         // DGELSS

         SRNAMT = 'DGELSS';
         INFOT = 1;
         dgelss(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO );
         chkxer('DGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgelss(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO );
         chkxer('DGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgelss(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO );
         chkxer('DGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgelss(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, INFO );
         chkxer('DGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dgelss(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, INFO );
         chkxer('DGELSS', INFOT, NOUT, LERR, OK );

         // DGELSY

         SRNAMT = 'DGELSY';
         INFOT = 1;
         dgelsy(-1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('DGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgelsy(0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('DGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgelsy(0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('DGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgelsy(2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('DGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dgelsy(2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('DGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dgelsy(2, 2, 1, A, 2, B, 2, IP, RCOND, IRNK, W, 1, INFO );
         chkxer('DGELSY', INFOT, NOUT, LERR, OK );

         // DGELSD

         SRNAMT = 'DGELSD';
         INFOT = 1;
         dgelsd(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('DGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgelsd(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('DGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgelsd(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('DGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgelsd(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('DGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dgelsd(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('DGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dgelsd(2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, IP, INFO );
         chkxer('DGELSD', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of DERRLS

      }
