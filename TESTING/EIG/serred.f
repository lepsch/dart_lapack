      void serred(PATH, NUNIT ) {

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
      REAL               ONE, ZERO;
      const              NMAX = 4, ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, IHI, ILO, INFO, J, NS, NT, SDIM;
      REAL               ABNRM;
      // ..
      // .. Local Arrays ..
      bool               B( NMAX );
      int                IW( 2*NMAX );
      REAL               A( NMAX, NMAX ), R1( NMAX ), R2( NMAX ), S( NMAX ), U( NMAX, NMAX ), VL( NMAX, NMAX ), VR( NMAX, NMAX ), VT( NMAX, NMAX ), W( 10*NMAX ), WI( NMAX ), WR( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, SGEES, SGEESX, SGEEV, SGEEVX, SGEJSV, SGESDD, SGESVD, SGESVDX, SGESVDQ
      // ..
      // .. External Functions ..
      bool               SSLECT, LSAMEN;
      // EXTERNAL SSLECT, LSAMEN
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

         // Test SGEEV

         SRNAMT = 'SGEEV ';
         INFOT = 1;
         sgeev('X', 'N', 0, A, 1, WR, WI, VL, 1, VR, 1, W, 1, INFO );
         chkxer('SGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgeev('N', 'X', 0, A, 1, WR, WI, VL, 1, VR, 1, W, 1, INFO );
         chkxer('SGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgeev('N', 'N', -1, A, 1, WR, WI, VL, 1, VR, 1, W, 1, INFO );
         chkxer('SGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgeev('N', 'N', 2, A, 1, WR, WI, VL, 1, VR, 1, W, 6, INFO );
         chkxer('SGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sgeev('V', 'N', 2, A, 2, WR, WI, VL, 1, VR, 1, W, 8, INFO );
         chkxer('SGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sgeev('N', 'V', 2, A, 2, WR, WI, VL, 1, VR, 1, W, 8, INFO );
         chkxer('SGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         sgeev('V', 'V', 1, A, 1, WR, WI, VL, 1, VR, 1, W, 3, INFO );
         chkxer('SGEEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

      } else if ( LSAMEN( 2, C2, 'ES' ) ) {

         // Test SGEES

         SRNAMT = 'SGEES ';
         INFOT = 1;
         sgees('X', 'N', SSLECT, 0, A, 1, SDIM, WR, WI, VL, 1, W, 1, B, INFO );
         chkxer('SGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgees('N', 'X', SSLECT, 0, A, 1, SDIM, WR, WI, VL, 1, W, 1, B, INFO );
         chkxer('SGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgees('N', 'S', SSLECT, -1, A, 1, SDIM, WR, WI, VL, 1, W, 1, B, INFO );
         chkxer('SGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgees('N', 'S', SSLECT, 2, A, 1, SDIM, WR, WI, VL, 1, W, 6, B, INFO );
         chkxer('SGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sgees('V', 'S', SSLECT, 2, A, 2, SDIM, WR, WI, VL, 1, W, 6, B, INFO );
         chkxer('SGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         sgees('N', 'S', SSLECT, 1, A, 1, SDIM, WR, WI, VL, 1, W, 2, B, INFO );
         chkxer('SGEES ', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

      } else if ( LSAMEN( 2, C2, 'VX' ) ) {

         // Test SGEEVX

         SRNAMT = 'SGEEVX';
         INFOT = 1;
         sgeevx('X', 'N', 'N', 'N', 0, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('SGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgeevx('N', 'X', 'N', 'N', 0, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('SGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgeevx('N', 'N', 'X', 'N', 0, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('SGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgeevx('N', 'N', 'N', 'X', 0, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('SGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgeevx('N', 'N', 'N', 'N', -1, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('SGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgeevx('N', 'N', 'N', 'N', 2, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('SGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sgeevx('N', 'V', 'N', 'N', 2, A, 2, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 6, IW, INFO );
         chkxer('SGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         sgeevx('N', 'N', 'V', 'N', 2, A, 2, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 6, IW, INFO );
         chkxer('SGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 21;
         sgeevx('N', 'N', 'N', 'N', 1, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('SGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 21;
         sgeevx('N', 'V', 'N', 'N', 1, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 2, IW, INFO );
         chkxer('SGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 21;
         sgeevx('N', 'N', 'V', 'V', 1, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 3, IW, INFO );
         chkxer('SGEEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

      } else if ( LSAMEN( 2, C2, 'SX' ) ) {

         // Test SGEESX

         SRNAMT = 'SGEESX';
         INFOT = 1;
         sgeesx('X', 'N', SSLECT, 'N', 0, A, 1, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 1, IW, 1, B, INFO );
         chkxer('SGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgeesx('N', 'X', SSLECT, 'N', 0, A, 1, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 1, IW, 1, B, INFO );
         chkxer('SGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgeesx('N', 'N', SSLECT, 'X', 0, A, 1, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 1, IW, 1, B, INFO );
         chkxer('SGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgeesx('N', 'N', SSLECT, 'N', -1, A, 1, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 1, IW, 1, B, INFO );
         chkxer('SGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgeesx('N', 'N', SSLECT, 'N', 2, A, 1, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 6, IW, 1, B, INFO );
         chkxer('SGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sgeesx('V', 'N', SSLECT, 'N', 2, A, 2, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 6, IW, 1, B, INFO );
         chkxer('SGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         sgeesx('N', 'N', SSLECT, 'N', 1, A, 1, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 2, IW, 1, B, INFO );
         chkxer('SGEESX', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

      } else if ( LSAMEN( 2, C2, 'BD' ) ) {

         // Test SGESVD

         SRNAMT = 'SGESVD';
         INFOT = 1;
         sgesvd('X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, INFO );
         chkxer('SGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgesvd('N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, INFO );
         chkxer('SGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgesvd('O', 'O', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, INFO );
         chkxer('SGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgesvd('N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, INFO );
         chkxer('SGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgesvd('N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, INFO );
         chkxer('SGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgesvd('N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, INFO );
         chkxer('SGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sgesvd('A', 'N', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, INFO );
         chkxer('SGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sgesvd('N', 'A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, INFO );
         chkxer('SGESVD', INFOT, NOUT, LERR, OK );
         NT = 8;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }

         // Test SGESDD

         SRNAMT = 'SGESDD';
         INFOT = 1;
         sgesdd('X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgesdd('N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgesdd('N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgesdd('N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, IW, INFO );
         chkxer('SGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgesdd('A', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, IW, INFO );
         chkxer('SGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sgesdd('A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, IW, INFO );
         chkxer('SGESDD', INFOT, NOUT, LERR, OK );
         NT = 6;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }

         // Test SGEJSV

         SRNAMT = 'SGEJSV';
         INFOT = 1;
         sgejsv('X', 'U', 'V', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('SGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgejsv('G', 'X', 'V', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('SGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgejsv('G', 'U', 'X', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('SGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgejsv('G', 'U', 'V', 'X', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('SGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgejsv('G', 'U', 'V', 'R', 'X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('SGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgejsv('G', 'U', 'V', 'R', 'N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('SGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgejsv('G', 'U', 'V', 'R', 'N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('SGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgejsv('G', 'U', 'V', 'R', 'N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('SGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('SGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         sgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 2, A, 2, S, U, 1, VT, 2, W, 1, IW, INFO);
         chkxer('SGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         sgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 2, A, 2, S, U, 2, VT, 1, W, 1, IW, INFO);
         chkxer('SGEJSV', INFOT, NOUT, LERR, OK );
         NT = 11;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }

         // Test SGESVDX

         SRNAMT = 'SGESVDX';
         INFOT = 1;
         sgesvdx('X', 'N', 'A', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgesvdx('N', 'X', 'A', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgesvdx('N', 'N', 'X', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgesvdx('N', 'N', 'A', -1, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgesvdx('N', 'N', 'A', 0, -1, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgesvdx('N', 'N', 'A', 2, 1, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sgesvdx('N', 'N', 'V', 2, 1, A, 2, -ONE, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sgesvdx('N', 'N', 'V', 2, 1, A, 2, ONE, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sgesvdx('N', 'N', 'I', 2, 2, A, 2, ZERO, ZERO, 0, 1, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sgesvdx('V', 'N', 'I', 2, 2, A, 2, ZERO, ZERO, 1, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         sgesvdx('V', 'N', 'A', 2, 2, A, 2, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         sgesvdx('N', 'V', 'A', 2, 2, A, 2, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('SGESVDX', INFOT, NOUT, LERR, OK );
         NT = 12;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }

         // Test SGESVDQ

         SRNAMT = 'SGESVDQ';
         INFOT = 1;
         sgesvdq('X', 'P', 'T', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('SGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgesvdq('A', 'X', 'T', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('SGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgesvdq('A', 'P', 'X', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('SGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgesvdq('A', 'P', 'T', 'X', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('SGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgesvdq('A', 'P', 'T', 'A', 'X', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('SGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sgesvdq('A', 'P', 'T', 'A', 'A', -1, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('SGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgesvdq('A', 'P', 'T', 'A', 'A', 0, 1, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('SGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 0, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('SGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, -1, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('SGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         sgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, 1, VT, -1, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('SGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         sgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, 1, VT, 1, NS, IW, -5, W, 1, W, 1, INFO );
         chkxer('SGESVDQ', INFOT, NOUT, LERR, OK );
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
      }
