      SUBROUTINE CERRGT( PATH, NUNIT )

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
      REAL               ANORM, RCOND
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      REAL               D( NMAX ), DF( NMAX ), R1( NMAX ), R2( NMAX ), RW( NMAX )       COMPLEX            B( NMAX ), DL( NMAX ), DLF( NMAX ), DU( NMAX ), DU2( NMAX ), DUF( NMAX ), E( NMAX ), EF( NMAX ), W( NMAX ), X( NMAX )
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CGTCON, CGTRFS, CGTTRF, CGTTRS, CHKXER, CPTCON, CPTRFS, CPTTRF, CPTTRS
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )
      DO 10 I = 1, NMAX
         D( I ) = 1.
         E( I ) = 2.
         DL( I ) = 3.
         DU( I ) = 4.
   10 CONTINUE
      ANORM = 1.0
      OK = .TRUE.

      if ( LSAMEN( 2, C2, 'GT' ) ) {

         // Test error exits for the general tridiagonal routines.

         // CGTTRF

         SRNAMT = 'CGTTRF'
         INFOT = 1
         cgttrf(-1, DL, E, DU, DU2, IP, INFO );
         chkxer('CGTTRF', INFOT, NOUT, LERR, OK );

         // CGTTRS

         SRNAMT = 'CGTTRS'
         INFOT = 1
         cgttrs('/', 0, 0, DL, E, DU, DU2, IP, X, 1, INFO );
         chkxer('CGTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cgttrs('N', -1, 0, DL, E, DU, DU2, IP, X, 1, INFO );
         chkxer('CGTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cgttrs('N', 0, -1, DL, E, DU, DU2, IP, X, 1, INFO );
         chkxer('CGTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 10
         cgttrs('N', 2, 1, DL, E, DU, DU2, IP, X, 1, INFO );
         chkxer('CGTTRS', INFOT, NOUT, LERR, OK );

         // CGTRFS

         SRNAMT = 'CGTRFS'
         INFOT = 1
         cgtrfs('/', 0, 0, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cgtrfs('N', -1, 0, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cgtrfs('N', 0, -1, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 13
         cgtrfs('N', 2, 1, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 1, X, 2, R1, R2, W, RW, INFO );
         chkxer('CGTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 15
         cgtrfs('N', 2, 1, DL, E, DU, DLF, EF, DUF, DU2, IP, B, 2, X, 1, R1, R2, W, RW, INFO );
         chkxer('CGTRFS', INFOT, NOUT, LERR, OK );

         // CGTCON

         SRNAMT = 'CGTCON'
         INFOT = 1
         cgtcon('/', 0, DL, E, DU, DU2, IP, ANORM, RCOND, W, INFO );
         chkxer('CGTCON', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cgtcon('I', -1, DL, E, DU, DU2, IP, ANORM, RCOND, W, INFO );
         chkxer('CGTCON', INFOT, NOUT, LERR, OK );
         INFOT = 8
         cgtcon('I', 0, DL, E, DU, DU2, IP, -ANORM, RCOND, W, INFO );
         chkxer('CGTCON', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'PT' ) ) {

         // Test error exits for the positive definite tridiagonal
         // routines.

         // CPTTRF

         SRNAMT = 'CPTTRF'
         INFOT = 1
         cpttrf(-1, D, E, INFO );
         chkxer('CPTTRF', INFOT, NOUT, LERR, OK );

         // CPTTRS

         SRNAMT = 'CPTTRS'
         INFOT = 1
         cpttrs('/', 1, 0, D, E, X, 1, INFO );
         chkxer('CPTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cpttrs('U', -1, 0, D, E, X, 1, INFO );
         chkxer('CPTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cpttrs('U', 0, -1, D, E, X, 1, INFO );
         chkxer('CPTTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7
         cpttrs('U', 2, 1, D, E, X, 1, INFO );
         chkxer('CPTTRS', INFOT, NOUT, LERR, OK );

         // CPTRFS

         SRNAMT = 'CPTRFS'
         INFOT = 1
         cptrfs('/', 1, 0, D, E, DF, EF, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         cptrfs('U', -1, 0, D, E, DF, EF, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         cptrfs('U', 0, -1, D, E, DF, EF, B, 1, X, 1, R1, R2, W, RW, INFO );
         chkxer('CPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 9
         cptrfs('U', 2, 1, D, E, DF, EF, B, 1, X, 2, R1, R2, W, RW, INFO );
         chkxer('CPTRFS', INFOT, NOUT, LERR, OK );
         INFOT = 11
         cptrfs('U', 2, 1, D, E, DF, EF, B, 2, X, 1, R1, R2, W, RW, INFO );
         chkxer('CPTRFS', INFOT, NOUT, LERR, OK );

         // CPTCON

         SRNAMT = 'CPTCON'
         INFOT = 1
         cptcon(-1, D, E, ANORM, RCOND, RW, INFO );
         chkxer('CPTCON', INFOT, NOUT, LERR, OK );
         INFOT = 4
         cptcon(0, D, E, -ANORM, RCOND, RW, INFO );
         chkxer('CPTCON', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of CERRGT

      }
