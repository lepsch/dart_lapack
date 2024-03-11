      void zerrgt(PATH, infoc.NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                infoc.NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 2 ;
      String             C2;
      int                I, INFO;
      double             ANORM, RCOND;
      int                IP( NMAX );
      double             D( NMAX ), DF( NMAX ), R1( NMAX ), R2( NMAX ), RW( NMAX );
      Complex         B( NMAX ), DL( NMAX ), DLF( NMAX ), DU( NMAX ), DU2( NMAX ), DUF( NMAX ), E( NMAX ), EF( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGTCON, ZGTRFS, ZGTTRF, ZGTTRS, ZPTCON, ZPTRFS, ZPTTRF, ZPTTRS
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT

      NOUT = infoc.NUNIT;
      NOUT.println( * );
      C2 = PATH.substring( 1, 3 );
      for (I = 1; I <= NMAX; I++) { // 10
         D[I] = 1.0;
         E[I] = 2.0;
         DL[I] = 3.0;
         DU[I] = 4.0;
      } // 10
      ANORM = 1.0;
      infoc.OK.value = true;

      if ( lsamen( 2, C2, 'GT' ) ) {

         // Test error exits for the general tridiagonal routines.

         // ZGTTRF

        srnamc.SRNAMT = 'ZGTTRF';
         infoc.INFOT = 1;
         zgttrf(-1, DL, E, DU, DU2, IP, INFO );
         chkxer('ZGTTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGTTRS

        srnamc.SRNAMT = 'ZGTTRS';
         infoc.INFOT = 1;
         zgttrs('/', 0, 0, DL, E, DU, DU2, IP, X, 1, INFO );
         chkxer('ZGTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgttrs('N', -1, 0, DL, E, DU, DU2, IP, X, 1, INFO );
         chkxer('ZGTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgttrs('N', 0, -1, DL, E, DU, DU2, IP, X, 1, INFO );
         chkxer('ZGTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zgttrs('N', 2, 1, DL, E, DU, DU2, IP, X, 1, INFO );
         chkxer('ZGTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGTRFS

        srnamc.SRNAMT = 'ZGTRFS';
         infoc.INFOT = 1;
         zgtrfs('/', 0, 0, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZGTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgtrfs('N', -1, 0, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZGTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zgtrfs('N', 0, -1, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZGTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 13;
         zgtrfs('N', 2, 1, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 1, X, 2, R1, R2, W, RW, INFO );
         chkxer('ZGTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 15;
         zgtrfs('N', 2, 1, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 2, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZGTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZGTCON

        srnamc.SRNAMT = 'ZGTCON';
         infoc.INFOT = 1;
         zgtcon('/', 0, DL, E, DU, DU2, IP, ANORM, RCOND, W, INFO );
         chkxer('ZGTCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zgtcon('I', -1, DL, E, DU, DU2, IP, ANORM, RCOND, W, INFO );
         chkxer('ZGTCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zgtcon('I', 0, DL, E, DU, DU2, IP, -ANORM, RCOND, W, INFO );
         chkxer('ZGTCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'PT' ) ) {

         // Test error exits for the positive definite tridiagonal
         // routines.

         // ZPTTRF

        srnamc.SRNAMT = 'ZPTTRF';
         infoc.INFOT = 1;
         zpttrf(-1, D, E, INFO );
         chkxer('ZPTTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPTTRS

        srnamc.SRNAMT = 'ZPTTRS';
         infoc.INFOT = 1;
         zpttrs('/', 1, 0, D, E, X, 1, INFO );
         chkxer('ZPTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zpttrs('U', -1, 0, D, E, X, 1, INFO );
         chkxer('ZPTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zpttrs('U', 0, -1, D, E, X, 1, INFO );
         chkxer('ZPTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zpttrs('U', 2, 1, D, E, X, 1, INFO );
         chkxer('ZPTTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPTRFS

        srnamc.SRNAMT = 'ZPTRFS';
         infoc.INFOT = 1;
         zptrfs('/', 1, 0, D, E, DF, EF, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZPTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zptrfs('U', -1, 0, D, E, DF, EF, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZPTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zptrfs('U', 0, -1, D, E, DF, EF, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZPTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         zptrfs('U', 2, 1, D, E, DF, EF, B, 1, X, 2, R1, R2, W, RW, INFO );
         chkxer('ZPTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         zptrfs('U', 2, 1, D, E, DF, EF, B, 2, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZPTRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZPTCON

        srnamc.SRNAMT = 'ZPTCON';
         infoc.INFOT = 1;
         zptcon(-1, D, E, ANORM, RCOND, RW, INFO );
         chkxer('ZPTCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zptcon(0, D, E, -ANORM, RCOND, RW, INFO );
         chkxer('ZPTCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      }
