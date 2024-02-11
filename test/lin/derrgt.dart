import 'common.dart';

      void derrgt(final int PATH, final int NUNIT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      String             C2;
      int                INFO;
      double             ANORM, RCOND;
      int                IP( NMAX ), IW( NMAX );
      double             B( NMAX ), C( NMAX ), CF( NMAX ), D( NMAX ), DF( NMAX ), E( NMAX ), EF( NMAX ), F( NMAX ), R1( NMAX ), R2( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DGTCON, DGTRFS, DGTTRF, DGTTRS, DPTCON, DPTRFS, DPTTRF, DPTTRS
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC / srnamc.SRNAMT

      infoc.NOUT = NUNIT;
      WRITE( infoc.NOUT, FMT = * );
      C2 = PATH( 2: 3 );
      D[1] = 1.0;
      D[2] = 2.0;
      DF[1] = 1.0;
      DF[2] = 2.0;
      E[1] = 3.0;
      E[2] = 4.0;
      EF[1] = 3.0;
      EF[2] = 4.0;
      ANORM = 1.0;
      infoc.OK = true;

      if ( lsamen( 2, C2, 'GT' ) ) {

         // Test error exits for the general tridiagonal routines.

         // DGTTRF

         srnamc.SRNAMT = 'DGTTRF';
         infoc.INFOT = 1;
         dgttrf(-1, C, D, E, F, IP, INFO );
         chkxer('DGTTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DGTTRS

         srnamc.SRNAMT = 'DGTTRS';
         infoc.INFOT = 1;
         dgttrs('/', 0, 0, C, D, E, F, IP, X, 1, INFO );
         chkxer('DGTTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgttrs('N', -1, 0, C, D, E, F, IP, X, 1, INFO );
         chkxer('DGTTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgttrs('N', 0, -1, C, D, E, F, IP, X, 1, INFO );
         chkxer('DGTTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dgttrs('N', 2, 1, C, D, E, F, IP, X, 1, INFO );
         chkxer('DGTTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DGTRFS

         srnamc.SRNAMT = 'DGTRFS';
         infoc.INFOT = 1;
         dgtrfs('/', 0, 0, C, D, E, CF, DF, EF, F, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DGTRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgtrfs('N', -1, 0, C, D, E, CF, DF, EF, F, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DGTRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         dgtrfs('N', 0, -1, C, D, E, CF, DF, EF, F, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DGTRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         dgtrfs('N', 2, 1, C, D, E, CF, DF, EF, F, IP, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('DGTRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 15;
         dgtrfs('N', 2, 1, C, D, E, CF, DF, EF, F, IP, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('DGTRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DGTCON

         srnamc.SRNAMT = 'DGTCON';
         infoc.INFOT = 1;
         dgtcon('/', 0, C, D, E, F, IP, ANORM, RCOND, W, IW, INFO );
         chkxer('DGTCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dgtcon('I', -1, C, D, E, F, IP, ANORM, RCOND, W, IW, INFO );
         chkxer('DGTCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dgtcon('I', 0, C, D, E, F, IP, -ANORM, RCOND, W, IW, INFO );
         chkxer('DGTCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'PT' ) ) {

         // Test error exits for the positive definite tridiagonal
         // routines.

         // DPTTRF

         srnamc.SRNAMT = 'DPTTRF';
         infoc.INFOT = 1;
         dpttrf(-1, D, E, INFO );
         chkxer('DPTTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DPTTRS

         srnamc.SRNAMT = 'DPTTRS';
         infoc.INFOT = 1;
         dpttrs(-1, 0, D, E, X, 1, INFO );
         chkxer('DPTTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dpttrs(0, -1, D, E, X, 1, INFO );
         chkxer('DPTTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         dpttrs(2, 1, D, E, X, 1, INFO );
         chkxer('DPTTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DPTRFS

         srnamc.SRNAMT = 'DPTRFS';
         infoc.INFOT = 1;
         dptrfs(-1, 0, D, E, DF, EF, B, 1, X, 1, R1, R2, W, INFO );
         chkxer('DPTRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         dptrfs(0, -1, D, E, DF, EF, B, 1, X, 1, R1, R2, W, INFO );
         chkxer('DPTRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         dptrfs(2, 1, D, E, DF, EF, B, 1, X, 2, R1, R2, W, INFO );
         chkxer('DPTRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         dptrfs(2, 1, D, E, DF, EF, B, 2, X, 1, R1, R2, W, INFO );
         chkxer('DPTRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );

         // DPTCON

         srnamc.SRNAMT = 'DPTCON';
         infoc.INFOT = 1;
         dptcon(-1, D, E, ANORM, RCOND, W, INFO );
         chkxer('DPTCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         dptcon(0, D, E, -ANORM, RCOND, W, INFO );
         chkxer('DPTCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, infoc.NOUT );

      }
