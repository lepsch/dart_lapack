      SUBROUTINE ZERRGT( PATH, NUNIT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX;
      const              NMAX = 2 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, INFO;
      double             ANORM, RCOND;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double             D( NMAX ), DF( NMAX ), R1( NMAX ), R2( NMAX ), RW( NMAX );
      COMPLEX*16         B( NMAX ), DL( NMAX ), DLF( NMAX ), DU( NMAX ), DU2( NMAX ), DUF( NMAX ), E( NMAX ), EF( NMAX ), W( NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZGTCON, ZGTRFS, ZGTTRF, ZGTTRS, ZPTCON, ZPTRFS, ZPTTRF, ZPTTRS
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

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
      for (I = 1; I <= NMAX; I++) { // 10
         D( I ) = 1.D0
         E( I ) = 2.D0
         DL( I ) = 3.D0
         DU( I ) = 4.D0
      } // 10
      ANORM = 1.0D0
      OK = .TRUE.

      if ( LSAMEN( 2, C2, 'GT' ) ) {

         // Test error exits for the general tridiagonal routines.

         // ZGTTRF

         SRNAMT = 'ZGTTRF'
         INFOT = 1
         zgttrf(-1, DL, E, DU, DU2, IP, INFO );
         chkxer('ZGTTRF', INFOT, NOUT, LERR, OK );

         // ZGTTRS

         SRNAMT = 'ZGTTRS'
         INFOT = 1
         zgttrs('/', 0, 0, DL, E, DU, DU2, IP, X, 1, INFO );
         chkxer('ZGTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgttrs('N', -1, 0, DL, E, DU, DU2, IP, X, 1, INFO );
         chkxer('ZGTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgttrs('N', 0, -1, DL, E, DU, DU2, IP, X, 1, INFO );
         chkxer('ZGTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zgttrs('N', 2, 1, DL, E, DU, DU2, IP, X, 1, INFO );
         chkxer('ZGTTRS', INFOT, NOUT, LERR, OK );

         // ZGTRFS

         SRNAMT = 'ZGTRFS'
         INFOT = 1
         zgtrfs('/', 0, 0, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgtrfs('N', -1, 0, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgtrfs('N', 0, -1, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 13
         zgtrfs('N', 2, 1, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 1, X, 2, R1, R2, W, RW, INFO );
         chkxer('ZGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 15
         zgtrfs('N', 2, 1, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 2, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZGTRFS', INFOT, NOUT, LERR, OK );

         // ZGTCON

         SRNAMT = 'ZGTCON'
         INFOT = 1
         zgtcon('/', 0, DL, E, DU, DU2, IP, ANORM, RCOND, W, INFO );
         chkxer('ZGTCON', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgtcon('I', -1, DL, E, DU, DU2, IP, ANORM, RCOND, W, INFO );
         chkxer('ZGTCON', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zgtcon('I', 0, DL, E, DU, DU2, IP, -ANORM, RCOND, W, INFO );
         chkxer('ZGTCON', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PT' ) ) {

         // Test error exits for the positive definite tridiagonal
         // routines.

         // ZPTTRF

         SRNAMT = 'ZPTTRF'
         INFOT = 1
         zpttrf(-1, D, E, INFO );
         chkxer('ZPTTRF', INFOT, NOUT, LERR, OK );

         // ZPTTRS

         SRNAMT = 'ZPTTRS'
         INFOT = 1
         zpttrs('/', 1, 0, D, E, X, 1, INFO );
         chkxer('ZPTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zpttrs('U', -1, 0, D, E, X, 1, INFO );
         chkxer('ZPTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zpttrs('U', 0, -1, D, E, X, 1, INFO );
         chkxer('ZPTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zpttrs('U', 2, 1, D, E, X, 1, INFO );
         chkxer('ZPTTRS', INFOT, NOUT, LERR, OK );

         // ZPTRFS

         SRNAMT = 'ZPTRFS'
         INFOT = 1
         zptrfs('/', 1, 0, D, E, DF, EF, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zptrfs('U', -1, 0, D, E, DF, EF, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zptrfs('U', 0, -1, D, E, DF, EF, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zptrfs('U', 2, 1, D, E, DF, EF, B, 1, X, 2, R1, R2, W, RW, INFO );
         chkxer('ZPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zptrfs('U', 2, 1, D, E, DF, EF, B, 2, X, 1, R1, R2, W, RW, INFO );
         chkxer('ZPTRFS', INFOT, NOUT, LERR, OK );

         // ZPTCON

         SRNAMT = 'ZPTCON'
         INFOT = 1
         zptcon(-1, D, E, ANORM, RCOND, RW, INFO );
         chkxer('ZPTCON', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zptcon(0, D, E, -ANORM, RCOND, RW, INFO );
         chkxer('ZPTCON', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of ZERRGT

      }
