      void serrgt(final int PATH, final int NUNIT) {

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
      double               ANORM, RCOND;
      int                IP( NMAX ), IW( NMAX );
      double               B( NMAX ), C( NMAX ), CF( NMAX ), D( NMAX ), DF( NMAX ), E( NMAX ), EF( NMAX ), F( NMAX ), R1( NMAX ), R2( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SGTCON, SGTRFS, SGTTRF, SGTTRS, SPTCON, SPTRFS, SPTTRF, SPTTRS
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
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );
      D[1] = 1.;
      D[2] = 2.;
      DF[1] = 1.;
      DF[2] = 2.;
      E[1] = 3.;
      E[2] = 4.;
      EF[1] = 3.;
      EF[2] = 4.;
      ANORM = 1.0;
      OK = true;

      if ( lsamen( 2, C2, 'GT' ) ) {

         // Test error exits for the general tridiagonal routines.

         // SGTTRF

        srnamc.SRNAMT = 'SGTTRF';
         INFOT = 1;
         sgttrf(-1, C, D, E, F, IP, INFO );
         chkxer('SGTTRF', INFOT, NOUT, LERR, OK );

         // SGTTRS

        srnamc.SRNAMT = 'SGTTRS';
         INFOT = 1;
         sgttrs('/', 0, 0, C, D, E, F, IP, X, 1, INFO );
         chkxer('SGTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgttrs('N', -1, 0, C, D, E, F, IP, X, 1, INFO );
         chkxer('SGTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgttrs('N', 0, -1, C, D, E, F, IP, X, 1, INFO );
         chkxer('SGTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sgttrs('N', 2, 1, C, D, E, F, IP, X, 1, INFO );
         chkxer('SGTTRS', INFOT, NOUT, LERR, OK );

         // SGTRFS

        srnamc.SRNAMT = 'SGTRFS';
         INFOT = 1;
         sgtrfs('/', 0, 0, C, D, E, CF, DF, EF, F, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgtrfs('N', -1, 0, C, D, E, CF, DF, EF, F, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgtrfs('N', 0, -1, C, D, E, CF, DF, EF, F, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         sgtrfs('N', 2, 1, C, D, E, CF, DF, EF, F, IP, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('SGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         sgtrfs('N', 2, 1, C, D, E, CF, DF, EF, F, IP, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('SGTRFS', INFOT, NOUT, LERR, OK );

         // SGTCON

        srnamc.SRNAMT = 'SGTCON';
         INFOT = 1;
         sgtcon('/', 0, C, D, E, F, IP, ANORM, RCOND, W, IW, INFO );
         chkxer('SGTCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgtcon('I', -1, C, D, E, F, IP, ANORM, RCOND, W, IW, INFO );
         chkxer('SGTCON', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgtcon('I', 0, C, D, E, F, IP, -ANORM, RCOND, W, IW, INFO );
         chkxer('SGTCON', INFOT, NOUT, LERR, OK );

      } else if ( lsamen( 2, C2, 'PT' ) ) {

         // Test error exits for the positive definite tridiagonal
         // routines.

         // SPTTRF

        srnamc.SRNAMT = 'SPTTRF';
         INFOT = 1;
         spttrf(-1, D, E, INFO );
         chkxer('SPTTRF', INFOT, NOUT, LERR, OK );

         // SPTTRS

        srnamc.SRNAMT = 'SPTTRS';
         INFOT = 1;
         spttrs(-1, 0, D, E, X, 1, INFO );
         chkxer('SPTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spttrs(0, -1, D, E, X, 1, INFO );
         chkxer('SPTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         spttrs(2, 1, D, E, X, 1, INFO );
         chkxer('SPTTRS', INFOT, NOUT, LERR, OK );

         // SPTRFS

        srnamc.SRNAMT = 'SPTRFS';
         INFOT = 1;
         sptrfs(-1, 0, D, E, DF, EF, B, 1, X, 1, R1, R2, W, INFO );
         chkxer('SPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sptrfs(0, -1, D, E, DF, EF, B, 1, X, 1, R1, R2, W, INFO );
         chkxer('SPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sptrfs(2, 1, D, E, DF, EF, B, 1, X, 2, R1, R2, W, INFO );
         chkxer('SPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sptrfs(2, 1, D, E, DF, EF, B, 2, X, 1, R1, R2, W, INFO );
         chkxer('SPTRFS', INFOT, NOUT, LERR, OK );

         // SPTCON

        srnamc.SRNAMT = 'SPTCON';
         INFOT = 1;
         sptcon(-1, D, E, ANORM, RCOND, W, INFO );
         chkxer('SPTCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sptcon(0, D, E, -ANORM, RCOND, W, INFO );
         chkxer('SPTCON', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      }
