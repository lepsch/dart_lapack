      void cerred(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX, LW;
      const              NMAX = 4, LW = 5*NMAX ;
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, IHI, ILO, INFO, J, NS, NT, SDIM;
      REAL               ABNRM;
      // ..
      // .. Local Arrays ..
      bool               B( NMAX );
      int                IW( 4*NMAX );
      REAL               R1( NMAX ), R2( NMAX ), RW( LW ), S( NMAX );
      COMPLEX            A( NMAX, NMAX ), U( NMAX, NMAX ), VL( NMAX, NMAX ), VR( NMAX, NMAX ), VT( NMAX, NMAX ), W( 10*NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, CGEES, CGEESX, CGEEV, CGEEVX, CGEJSV, CGESDD, CGESVD, CGESVDX, CGESVDQ
      // ..
      // .. External Functions ..
      bool               LSAMEN, CSLECT;
      // EXTERNAL LSAMEN, CSLECT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LEN_TRIM
      // ..
      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      REAL               SELWI( 20 ), SELWR( 20 );
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT, SELDIM, SELOPT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // COMMON / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Initialize A

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = ZERO;
         } // 10
      } // 20
      for (I = 1; I <= NMAX; I++) { // 30
         A( I, I ) = ONE;
      } // 30
      OK = true;
      NT = 0;

      if ( LSAMEN( 2, C2, 'EV' ) ) {

         // Test CGEEV

         SRNAMT = 'CGEEV ';
         INFOT = 1;
         cgeev('X', 'N', 0, A, 1, X, VL, 1, VR, 1, W, 1, RW, INFO );
         chkxer('CGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgeev('N', 'X', 0, A, 1, X, VL, 1, VR, 1, W, 1, RW, INFO );
         chkxer('CGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgeev('N', 'N', -1, A, 1, X, VL, 1, VR, 1, W, 1, RW, INFO );
         chkxer('CGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgeev('N', 'N', 2, A, 1, X, VL, 1, VR, 1, W, 4, RW, INFO );
         chkxer('CGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgeev('V', 'N', 2, A, 2, X, VL, 1, VR, 1, W, 4, RW, INFO );
         chkxer('CGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cgeev('N', 'V', 2, A, 2, X, VL, 1, VR, 1, W, 4, RW, INFO );
         chkxer('CGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cgeev('V', 'V', 1, A, 1, X, VL, 1, VR, 1, W, 1, RW, INFO );
         chkxer('CGEEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

      } else if ( LSAMEN( 2, C2, 'ES' ) ) {

         // Test CGEES

         SRNAMT = 'CGEES ';
         INFOT = 1;
         cgees('X', 'N', CSLECT, 0, A, 1, SDIM, X, VL, 1, W, 1, RW, B, INFO );
         chkxer('CGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgees('N', 'X', CSLECT, 0, A, 1, SDIM, X, VL, 1, W, 1, RW, B, INFO );
         chkxer('CGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgees('N', 'S', CSLECT, -1, A, 1, SDIM, X, VL, 1, W, 1, RW, B, INFO );
         chkxer('CGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgees('N', 'S', CSLECT, 2, A, 1, SDIM, X, VL, 1, W, 4, RW, B, INFO );
         chkxer('CGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cgees('V', 'S', CSLECT, 2, A, 2, SDIM, X, VL, 1, W, 4, RW, B, INFO );
         chkxer('CGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cgees('N', 'S', CSLECT, 1, A, 1, SDIM, X, VL, 1, W, 1, RW, B, INFO );
         chkxer('CGEES ', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

      } else if ( LSAMEN( 2, C2, 'VX' ) ) {

         // Test CGEEVX

         SRNAMT = 'CGEEVX';
         INFOT = 1;
         cgeevx('X', 'N', 'N', 'N', 0, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO );
         chkxer('CGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgeevx('N', 'X', 'N', 'N', 0, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO );
         chkxer('CGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgeevx('N', 'N', 'X', 'N', 0, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO );
         chkxer('CGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgeevx('N', 'N', 'N', 'X', 0, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO );
         chkxer('CGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgeevx('N', 'N', 'N', 'N', -1, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO );
         chkxer('CGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgeevx('N', 'N', 'N', 'N', 2, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 4, RW, INFO );
         chkxer('CGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cgeevx('N', 'V', 'N', 'N', 2, A, 2, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 4, RW, INFO );
         chkxer('CGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cgeevx('N', 'N', 'V', 'N', 2, A, 2, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 4, RW, INFO );
         chkxer('CGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         cgeevx('N', 'N', 'N', 'N', 1, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO );
         chkxer('CGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         cgeevx('N', 'N', 'V', 'V', 1, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 2, RW, INFO );
         chkxer('CGEEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

      } else if ( LSAMEN( 2, C2, 'SX' ) ) {

         // Test CGEESX

         SRNAMT = 'CGEESX';
         INFOT = 1;
         cgeesx('X', 'N', CSLECT, 'N', 0, A, 1, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 1, RW, B, INFO );
         chkxer('CGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgeesx('N', 'X', CSLECT, 'N', 0, A, 1, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 1, RW, B, INFO );
         chkxer('CGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgeesx('N', 'N', CSLECT, 'X', 0, A, 1, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 1, RW, B, INFO );
         chkxer('CGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgeesx('N', 'N', CSLECT, 'N', -1, A, 1, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 1, RW, B, INFO );
         chkxer('CGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgeesx('N', 'N', CSLECT, 'N', 2, A, 1, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 4, RW, B, INFO );
         chkxer('CGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cgeesx('V', 'N', CSLECT, 'N', 2, A, 2, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 4, RW, B, INFO );
         chkxer('CGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cgeesx('N', 'N', CSLECT, 'N', 1, A, 1, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 1, RW, B, INFO );
         chkxer('CGEESX', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

      } else if ( LSAMEN( 2, C2, 'BD' ) ) {

         // Test CGESVD

         SRNAMT = 'CGESVD';
         INFOT = 1;
         cgesvd('X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO );
         chkxer('CGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgesvd('N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO );
         chkxer('CGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgesvd('O', 'O', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO );
         chkxer('CGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgesvd('N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO );
         chkxer('CGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgesvd('N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO );
         chkxer('CGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgesvd('N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, RW, INFO );
         chkxer('CGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cgesvd('A', 'N', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, RW, INFO );
         chkxer('CGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cgesvd('N', 'A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, RW, INFO );
         chkxer('CGESVD', INFOT, NOUT, LERR, OK );
         NT = NT + 8;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }

         // Test CGESDD

         SRNAMT = 'CGESDD';
         INFOT = 1;
         cgesdd('X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgesdd('N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgesdd('N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgesdd('N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, RW, IW, INFO );
         chkxer('CGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgesdd('A', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, RW, IW, INFO );
         chkxer('CGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cgesdd('A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, RW, IW, INFO );
         chkxer('CGESDD', INFOT, NOUT, LERR, OK );
         NT = NT - 2;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }

         // Test CGEJSV

         SRNAMT = 'CGEJSV';
         INFOT = 1;
         cgejsv('X', 'U', 'V', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('CGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgejsv('G', 'X', 'V', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('CGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgejsv('G', 'U', 'X', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('CGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgejsv('G', 'U', 'V', 'X', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('CGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgejsv('G', 'U', 'V', 'R', 'X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('CGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgejsv('G', 'U', 'V', 'R', 'N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('CGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgejsv('G', 'U', 'V', 'R', 'N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('CGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgejsv('G', 'U', 'V', 'R', 'N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('CGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('CGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         cgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 2, A, 2, S, U, 1, VT, 2, W, 1, RW, 1, IW, INFO);
         chkxer('CGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 2, A, 2, S, U, 2, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('CGEJSV', INFOT, NOUT, LERR, OK );
         NT = 11;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }

         // Test CGESVDX

         SRNAMT = 'CGESVDX';
         INFOT = 1;
         cgesvdx('X', 'N', 'A', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgesvdx('N', 'X', 'A', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgesvdx('N', 'N', 'X', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgesvdx('N', 'N', 'A', -1, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgesvdx('N', 'N', 'A', 0, -1, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgesvdx('N', 'N', 'A', 2, 1, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cgesvdx('N', 'N', 'V', 2, 1, A, 2, -ONE, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cgesvdx('N', 'N', 'V', 2, 1, A, 2, ONE, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cgesvdx('N', 'N', 'I', 2, 2, A, 2, ZERO, ZERO, 0, 1, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cgesvdx('V', 'N', 'I', 2, 2, A, 2, ZERO, ZERO, 1, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cgesvdx('V', 'N', 'A', 2, 2, A, 2, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         cgesvdx('N', 'V', 'A', 2, 2, A, 2, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('CGESVDX', INFOT, NOUT, LERR, OK );
         NT = 12;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }

         // Test CGESVDQ

         SRNAMT = 'CGESVDQ';
         INFOT = 1;
         cgesvdq('X', 'P', 'T', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('CGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgesvdq('A', 'X', 'T', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('CGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgesvdq('A', 'P', 'X', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('CGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgesvdq('A', 'P', 'T', 'X', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('CGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgesvdq('A', 'P', 'T', 'A', 'X', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('CGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cgesvdq('A', 'P', 'T', 'A', 'A', -1, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('CGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgesvdq('A', 'P', 'T', 'A', 'A', 0, 1, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('CGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 0, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('CGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, -1, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('CGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         cgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, 1, VT, -1, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('CGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         cgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, 1, VT, 1, NS, IW, -5, W, 1, RW, 1, INFO );
         chkxer('CGESVDQ', INFOT, NOUT, LERR, OK );
         NT = 11;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }
      }

      // Print a summary line.

      if ( !LSAMEN( 2, C2, 'BD' ) ) {
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }
      }

 9999 FORMAT( 1X, A, ' passed the tests of the error exits (', I3, ' tests done)' );
 9998 FORMAT( ' *** ', A, ' failed the tests of the error exits ***' );
      return;

      // End of CERRED

      }
