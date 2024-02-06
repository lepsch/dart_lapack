      void zerrls(PATH, infoc.NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                infoc.NUNIT;
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
      double             RW( NMAX ), S( NMAX );
      Complex         A( NMAX, NMAX ), B( NMAX, NMAX ), W( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGELS, ZGELSD, ZGELSS, ZGELST, ZGELSY, ZGETSLS
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
      // .. Executable Statements ..

      NOUT = infoc.NUNIT;
      C2 = PATH( 2: 3 );
      A[1][1] = ( 1.0, 0.0 );
      A[1][2] = ( 2.0, 0.0 );
      A[2][2] = ( 3.0, 0.0 );
      A[2][1] = ( 4.0, 0.0 );
      infoc.OK = true;
      WRITE( NOUT, FMT = * );

      // Test error exits for the least squares driver routines.

      if ( lsamen( 2, C2, 'LS' ) ) {

         // ZGELS

        srnamc.SRNAMT = 'ZGELS ';
         infoc.INFOT = 1;
         zgels('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELS ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgels('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELS ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgels('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELS ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgels('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELS ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zgels('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('ZGELS ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zgels('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('ZGELS ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zgels('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('ZGELS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zgels('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELS ', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGELST

        srnamc.SRNAMT = 'ZGELST';
         infoc.INFOT = 1;
         zgelst('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgelst('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgelst('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgelst('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zgelst('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zgelst('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zgelst('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zgelst('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGELST', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGETSLS

        srnamc.SRNAMT = 'ZGETSLS';
         infoc.INFOT = 1;
         zgetsls('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgetsls('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgetsls('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zgetsls('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zgetsls('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zgetsls('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zgetsls('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('ZGETSLS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGELSS

        srnamc.SRNAMT = 'ZGELSS';
         infoc.INFOT = 1;
         zgelss(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO );
         chkxer('ZGELSS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgelss(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO );
         chkxer('ZGELSS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgelss(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, RW, INFO );
         chkxer('ZGELSS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zgelss(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, RW, INFO );
         chkxer('ZGELSS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zgelss(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, RW, INFO );
         chkxer('ZGELSS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGELSY

        srnamc.SRNAMT = 'ZGELSY';
         infoc.INFOT = 1;
         zgelsy(-1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('ZGELSY', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgelsy(0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('ZGELSY', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgelsy(0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('ZGELSY', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zgelsy(2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('ZGELSY', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zgelsy(2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10, RW, INFO );
         chkxer('ZGELSY', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         zgelsy(0, 3, 0, A, 1, B, 3, IP, RCOND, IRNK, W, 1, RW, INFO );
         chkxer('ZGELSY', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGELSD

        srnamc.SRNAMT = 'ZGELSD';
         infoc.INFOT = 1;
         zgelsd(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('ZGELSD', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgelsd(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('ZGELSD', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgelsd(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('ZGELSD', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zgelsd(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('ZGELSD', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zgelsd(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, RW, IP, INFO );
         chkxer('ZGELSD', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         zgelsd(2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, RW, IP, INFO );
         chkxer('ZGELSD', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      return;
      }
