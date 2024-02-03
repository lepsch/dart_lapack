      SUBROUTINE SERRLS( PATH, NUNIT );

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
      REAL               RCOND;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      REAL               A( NMAX, NMAX ), B( NMAX, NMAX ), S( NMAX ), W( NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SGELS, SGELSD, SGELSS, SGELST, SGELSY, SGETSLS
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

         // SGELS

         SRNAMT = 'SGELS ';
         INFOT = 1;
         sgels('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('SGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgels('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('SGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgels('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('SGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgels('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('SGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgels('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('SGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgels('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('SGELS ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgels('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('DGELS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sgels('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('SGELS ', INFOT, NOUT, LERR, OK );

         // SGELST

         SRNAMT = 'SGELST';
         INFOT = 1;
         sgelst('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('SGELST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgelst('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('SGELST', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgelst('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('SGELST', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgelst('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('SGELST', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgelst('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('SGELST', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgelst('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('SGELST', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgelst('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('SGELST', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sgelst('N', 1, 1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('SGELST', INFOT, NOUT, LERR, OK );

         // SGETSLS

         SRNAMT = 'SGETSLS';
         INFOT = 1;
         sgetsls('/', 0, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('SGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgetsls('N', -1, 0, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('SGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgetsls('N', 0, -1, 0, A, 1, B, 1, W, 1, INFO );
         chkxer('SGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgetsls('N', 0, 0, -1, A, 1, B, 1, W, 1, INFO );
         chkxer('SGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgetsls('N', 2, 0, 0, A, 1, B, 2, W, 2, INFO );
         chkxer('SGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgetsls('N', 2, 0, 0, A, 2, B, 1, W, 2, INFO );
         chkxer('SGETSLS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgetsls('N', 0, 2, 0, A, 1, B, 1, W, 2, INFO );
         chkxer('SGETSLS', INFOT, NOUT, LERR, OK );

         // SGELSS

         SRNAMT = 'SGELSS';
         INFOT = 1;
         sgelss(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO );
         chkxer('SGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgelss(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO );
         chkxer('SGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgelss(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 1, INFO );
         chkxer('SGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgelss(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 2, INFO );
         chkxer('SGELSS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgelss(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 2, INFO );
         chkxer('SGELSS', INFOT, NOUT, LERR, OK );

         // SGELSY

         SRNAMT = 'SGELSY';
         INFOT = 1;
         sgelsy(-1, 0, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('SGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgelsy(0, -1, 0, A, 1, B, 1, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('SGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgelsy(0, 0, -1, A, 1, B, 1, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('SGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgelsy(2, 0, 0, A, 1, B, 2, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('SGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgelsy(2, 0, 0, A, 2, B, 1, IP, RCOND, IRNK, W, 10, INFO );
         chkxer('SGELSY', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sgelsy(2, 2, 1, A, 2, B, 2, IP, RCOND, IRNK, W, 1, INFO );
         chkxer('SGELSY', INFOT, NOUT, LERR, OK );

         // SGELSD

         SRNAMT = 'SGELSD';
         INFOT = 1;
         sgelsd(-1, 0, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('SGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgelsd(0, -1, 0, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('SGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgelsd(0, 0, -1, A, 1, B, 1, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('SGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgelsd(2, 0, 0, A, 1, B, 2, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('SGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgelsd(2, 0, 0, A, 2, B, 1, S, RCOND, IRNK, W, 10, IP, INFO );
         chkxer('SGELSD', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sgelsd(2, 2, 1, A, 2, B, 2, S, RCOND, IRNK, W, 1, IP, INFO );
         chkxer('SGELSD', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of SERRLS

      }
