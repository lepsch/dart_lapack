      void zerred(PATH, NUNIT ) {

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
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, IHI, ILO, INFO, J, NS, NT, SDIM;
      double             ABNRM;
      // ..
      // .. Local Arrays ..
      bool               B( NMAX );
      int                IW( 4*NMAX );
      double             R1( NMAX ), R2( NMAX ), RW( LW ), S( NMAX );
      Complex         A( NMAX, NMAX ), U( NMAX, NMAX ), VL( NMAX, NMAX ), VR( NMAX, NMAX ), VT( NMAX, NMAX ), W( 10*NMAX ), X( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZGEES, ZGEESX, ZGEEV, ZGEEVX, ZGESVJ, ZGESDD, ZGESVD, ZGESVDX, ZGESVDQ
      // ..
      // .. External Functions ..
      //- bool               LSAMEN, ZSLECT;
      // EXTERNAL LSAMEN, ZSLECT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LEN_TRIM
      // ..
      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      double             SELWI( 20 ), SELWR( 20 );
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
            A[I, J] = ZERO;
         } // 10
      } // 20
      for (I = 1; I <= NMAX; I++) { // 30
         A[I, I] = ONE;
      } // 30
      OK = true;
      NT = 0;

      if ( lsamen( 2, C2, 'EV' ) ) {

         // Test ZGEEV

         SRNAMT = 'ZGEEV ';
         INFOT = 1;
         zgeev('X', 'N', 0, A, 1, X, VL, 1, VR, 1, W, 1, RW, INFO );
         chkxer('ZGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgeev('N', 'X', 0, A, 1, X, VL, 1, VR, 1, W, 1, RW, INFO );
         chkxer('ZGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgeev('N', 'N', -1, A, 1, X, VL, 1, VR, 1, W, 1, RW, INFO );
         chkxer('ZGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgeev('N', 'N', 2, A, 1, X, VL, 1, VR, 1, W, 4, RW, INFO );
         chkxer('ZGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zgeev('V', 'N', 2, A, 2, X, VL, 1, VR, 1, W, 4, RW, INFO );
         chkxer('ZGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zgeev('N', 'V', 2, A, 2, X, VL, 1, VR, 1, W, 4, RW, INFO );
         chkxer('ZGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zgeev('V', 'V', 1, A, 1, X, VL, 1, VR, 1, W, 1, RW, INFO );
         chkxer('ZGEEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

      } else if ( lsamen( 2, C2, 'ES' ) ) {

         // Test ZGEES

         SRNAMT = 'ZGEES ';
         INFOT = 1;
         zgees('X', 'N', ZSLECT, 0, A, 1, SDIM, X, VL, 1, W, 1, RW, B, INFO );
         chkxer('ZGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgees('N', 'X', ZSLECT, 0, A, 1, SDIM, X, VL, 1, W, 1, RW, B, INFO );
         chkxer('ZGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgees('N', 'S', ZSLECT, -1, A, 1, SDIM, X, VL, 1, W, 1, RW, B, INFO );
         chkxer('ZGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zgees('N', 'S', ZSLECT, 2, A, 1, SDIM, X, VL, 1, W, 4, RW, B, INFO );
         chkxer('ZGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zgees('V', 'S', ZSLECT, 2, A, 2, SDIM, X, VL, 1, W, 4, RW, B, INFO );
         chkxer('ZGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zgees('N', 'S', ZSLECT, 1, A, 1, SDIM, X, VL, 1, W, 1, RW, B, INFO );
         chkxer('ZGEES ', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

      } else if ( lsamen( 2, C2, 'VX' ) ) {

         // Test ZGEEVX

         SRNAMT = 'ZGEEVX';
         INFOT = 1;
         zgeevx('X', 'N', 'N', 'N', 0, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO );
         chkxer('ZGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgeevx('N', 'X', 'N', 'N', 0, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO );
         chkxer('ZGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgeevx('N', 'N', 'X', 'N', 0, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO );
         chkxer('ZGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgeevx('N', 'N', 'N', 'X', 0, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO );
         chkxer('ZGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgeevx('N', 'N', 'N', 'N', -1, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO );
         chkxer('ZGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zgeevx('N', 'N', 'N', 'N', 2, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 4, RW, INFO );
         chkxer('ZGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zgeevx('N', 'V', 'N', 'N', 2, A, 2, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 4, RW, INFO );
         chkxer('ZGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zgeevx('N', 'N', 'V', 'N', 2, A, 2, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 4, RW, INFO );
         chkxer('ZGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         zgeevx('N', 'N', 'N', 'N', 1, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, RW, INFO );
         chkxer('ZGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         zgeevx('N', 'N', 'V', 'V', 1, A, 1, X, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 2, RW, INFO );
         chkxer('ZGEEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

      } else if ( lsamen( 2, C2, 'SX' ) ) {

         // Test ZGEESX

         SRNAMT = 'ZGEESX';
         INFOT = 1;
         zgeesx('X', 'N', ZSLECT, 'N', 0, A, 1, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 1, RW, B, INFO );
         chkxer('ZGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgeesx('N', 'X', ZSLECT, 'N', 0, A, 1, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 1, RW, B, INFO );
         chkxer('ZGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgeesx('N', 'N', ZSLECT, 'X', 0, A, 1, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 1, RW, B, INFO );
         chkxer('ZGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgeesx('N', 'N', ZSLECT, 'N', -1, A, 1, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 1, RW, B, INFO );
         chkxer('ZGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zgeesx('N', 'N', ZSLECT, 'N', 2, A, 1, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 4, RW, B, INFO );
         chkxer('ZGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zgeesx('V', 'N', ZSLECT, 'N', 2, A, 2, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 4, RW, B, INFO );
         chkxer('ZGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zgeesx('N', 'N', ZSLECT, 'N', 1, A, 1, SDIM, X, VL, 1, R1( 1 ), R2( 1 ), W, 1, RW, B, INFO );
         chkxer('ZGEESX', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

      } else if ( lsamen( 2, C2, 'BD' ) ) {

         // Test ZGESVD

         SRNAMT = 'ZGESVD';
         INFOT = 1;
         zgesvd('X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO );
         chkxer('ZGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgesvd('N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO );
         chkxer('ZGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgesvd('O', 'O', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO );
         chkxer('ZGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgesvd('N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO );
         chkxer('ZGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgesvd('N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, RW, INFO );
         chkxer('ZGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zgesvd('N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, RW, INFO );
         chkxer('ZGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zgesvd('A', 'N', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, RW, INFO );
         chkxer('ZGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zgesvd('N', 'A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, RW, INFO );
         chkxer('ZGESVD', INFOT, NOUT, LERR, OK );
         NT = NT + 8;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }

         // Test ZGESDD

         SRNAMT = 'ZGESDD';
         INFOT = 1;
         zgesdd('X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgesdd('N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgesdd('N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgesdd('N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, RW, IW, INFO );
         chkxer('ZGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zgesdd('A', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, RW, IW, INFO );
         chkxer('ZGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zgesdd('A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, RW, IW, INFO );
         chkxer('ZGESDD', INFOT, NOUT, LERR, OK );
         NT = NT - 2;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }

         // Test ZGEJSV

         SRNAMT = 'ZGEJSV';
         INFOT = 1;
         zgejsv('X', 'U', 'V', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('ZGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgejsv('G', 'X', 'V', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('ZGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgejsv('G', 'U', 'X', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('ZGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgejsv('G', 'U', 'V', 'X', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('ZGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgejsv('G', 'U', 'V', 'R', 'X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('ZGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zgejsv('G', 'U', 'V', 'R', 'N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('ZGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zgejsv('G', 'U', 'V', 'R', 'N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('ZGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zgejsv('G', 'U', 'V', 'R', 'N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('ZGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('ZGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 2, A, 2, S, U, 1, VT, 2, W, 1, RW, 1, IW, INFO);
         chkxer('ZGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 2, A, 2, S, U, 2, VT, 1, W, 1, RW, 1, IW, INFO);
         chkxer('ZGEJSV', INFOT, NOUT, LERR, OK );
         NT = 11;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }

         // Test ZGESVDX

         SRNAMT = 'ZGESVDX';
         INFOT = 1;
         zgesvdx('X', 'N', 'A', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgesvdx('N', 'X', 'A', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgesvdx('N', 'N', 'X', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgesvdx('N', 'N', 'A', -1, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgesvdx('N', 'N', 'A', 0, -1, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zgesvdx('N', 'N', 'A', 2, 1, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zgesvdx('N', 'N', 'V', 2, 1, A, 2, -ONE, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zgesvdx('N', 'N', 'V', 2, 1, A, 2, ONE, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zgesvdx('N', 'N', 'I', 2, 2, A, 2, ZERO, ZERO, 0, 1, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zgesvdx('V', 'N', 'I', 2, 2, A, 2, ZERO, ZERO, 1, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zgesvdx('V', 'N', 'A', 2, 2, A, 2, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         zgesvdx('N', 'V', 'A', 2, 2, A, 2, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, RW, IW, INFO );
         chkxer('ZGESVDX', INFOT, NOUT, LERR, OK );
         NT = 12;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }

         // Test ZGESVDQ

         SRNAMT = 'ZGESVDQ';
         INFOT = 1;
         zgesvdq('X', 'P', 'T', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('ZGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zgesvdq('A', 'X', 'T', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('ZGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zgesvdq('A', 'P', 'X', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('ZGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zgesvdq('A', 'P', 'T', 'X', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('ZGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zgesvdq('A', 'P', 'T', 'A', 'X', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('ZGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zgesvdq('A', 'P', 'T', 'A', 'A', -1, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('ZGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zgesvdq('A', 'P', 'T', 'A', 'A', 0, 1, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('ZGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 0, S, U, 0, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('ZGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, -1, VT, 0, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('ZGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         zgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, 1, VT, -1, NS, IW, 1, W, 1, RW, 1, INFO );
         chkxer('ZGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         zgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, 1, VT, 1, NS, IW, -5, W, 1, RW, 1, INFO );
         chkxer('ZGESVDQ', INFOT, NOUT, LERR, OK );
         NT = 11;
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT;
         } else {
            WRITE( NOUT, FMT = 9998 );
         }
      }

      // Print a summary line.

      if ( !lsamen( 2, C2, 'BD' ) ) {
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
