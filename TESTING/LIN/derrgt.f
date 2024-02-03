      void derrgt(PATH, NUNIT ) {

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
      int                INFO;
      double             ANORM, RCOND;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX ), IW( NMAX );
      double             B( NMAX ), C( NMAX ), CF( NMAX ), D( NMAX ), DF( NMAX ), E( NMAX ), EF( NMAX ), F( NMAX ), R1( NMAX ), R2( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGTCON, DGTRFS, DGTTRF, DGTTRS, DPTCON, DPTRFS, DPTTRF, DPTTRS
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
      D( 1 ) = 1.0;
      D( 2 ) = 2.0;
      DF( 1 ) = 1.0;
      DF( 2 ) = 2.0;
      E( 1 ) = 3.0;
      E( 2 ) = 4.0;
      EF( 1 ) = 3.0;
      EF( 2 ) = 4.0;
      ANORM = 1.0;
      OK = true;

      if ( LSAMEN( 2, C2, 'GT' ) ) {

         // Test error exits for the general tridiagonal routines.

         // DGTTRF

         SRNAMT = 'DGTTRF';
         INFOT = 1;
         dgttrf(-1, C, D, E, F, IP, INFO );
         chkxer('DGTTRF', INFOT, NOUT, LERR, OK );

         // DGTTRS

         SRNAMT = 'DGTTRS';
         INFOT = 1;
         dgttrs('/', 0, 0, C, D, E, F, IP, X, 1, INFO );
         chkxer('DGTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgttrs('N', -1, 0, C, D, E, F, IP, X, 1, INFO );
         chkxer('DGTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgttrs('N', 0, -1, C, D, E, F, IP, X, 1, INFO );
         chkxer('DGTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dgttrs('N', 2, 1, C, D, E, F, IP, X, 1, INFO );
         chkxer('DGTTRS', INFOT, NOUT, LERR, OK );

         // DGTRFS

         SRNAMT = 'DGTRFS';
         INFOT = 1;
         dgtrfs('/', 0, 0, C, D, E, CF, DF, EF, F, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgtrfs('N', -1, 0, C, D, E, CF, DF, EF, F, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgtrfs('N', 0, -1, C, D, E, CF, DF, EF, F, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dgtrfs('N', 2, 1, C, D, E, CF, DF, EF, F, IP, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('DGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dgtrfs('N', 2, 1, C, D, E, CF, DF, EF, F, IP, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('DGTRFS', INFOT, NOUT, LERR, OK );

         // DGTCON

         SRNAMT = 'DGTCON';
         INFOT = 1;
         dgtcon('/', 0, C, D, E, F, IP, ANORM, RCOND, W, IW, INFO );
         chkxer('DGTCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgtcon('I', -1, C, D, E, F, IP, ANORM, RCOND, W, IW, INFO );
         chkxer('DGTCON', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dgtcon('I', 0, C, D, E, F, IP, -ANORM, RCOND, W, IW, INFO );
         chkxer('DGTCON', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PT' ) ) {

         // Test error exits for the positive definite tridiagonal
         // routines.

         // DPTTRF

         SRNAMT = 'DPTTRF';
         INFOT = 1;
         dpttrf(-1, D, E, INFO );
         chkxer('DPTTRF', INFOT, NOUT, LERR, OK );

         // DPTTRS

         SRNAMT = 'DPTTRS';
         INFOT = 1;
         dpttrs(-1, 0, D, E, X, 1, INFO );
         chkxer('DPTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpttrs(0, -1, D, E, X, 1, INFO );
         chkxer('DPTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dpttrs(2, 1, D, E, X, 1, INFO );
         chkxer('DPTTRS', INFOT, NOUT, LERR, OK );

         // DPTRFS

         SRNAMT = 'DPTRFS';
         INFOT = 1;
         dptrfs(-1, 0, D, E, DF, EF, B, 1, X, 1, R1, R2, W, INFO );
         chkxer('DPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dptrfs(0, -1, D, E, DF, EF, B, 1, X, 1, R1, R2, W, INFO );
         chkxer('DPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dptrfs(2, 1, D, E, DF, EF, B, 1, X, 2, R1, R2, W, INFO );
         chkxer('DPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dptrfs(2, 1, D, E, DF, EF, B, 2, X, 1, R1, R2, W, INFO );
         chkxer('DPTRFS', INFOT, NOUT, LERR, OK );

         // DPTCON

         SRNAMT = 'DPTCON';
         INFOT = 1;
         dptcon(-1, D, E, ANORM, RCOND, W, INFO );
         chkxer('DPTCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dptcon(0, D, E, -ANORM, RCOND, W, INFO );
         chkxer('DPTCON', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of DERRGT

      }
