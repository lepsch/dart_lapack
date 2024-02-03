      SUBROUTINE DERRED( PATH, NUNIT )

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
      double             ONE, ZERO;
      const              NMAX = 4, ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, IHI, ILO, INFO, J, NS, NT, SDIM;
      double             ABNRM;
      // ..
      // .. Local Arrays ..
      bool               B( NMAX );
      int                IW( 2*NMAX );
      double             A( NMAX, NMAX ), R1( NMAX ), R2( NMAX ), S( NMAX ), U( NMAX, NMAX ), VL( NMAX, NMAX ), VR( NMAX, NMAX ), VT( NMAX, NMAX ), W( 10*NMAX ), WI( NMAX ), WR( NMAX );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, DGEES, DGEESX, DGEEV, DGEEVX, DGEJSV, DGESDD, DGESVD, DGESVDX, DGESVDQ
      // ..
      // .. External Functions ..
      bool               DSLECT, LSAMEN;
      // EXTERNAL DSLECT, LSAMEN
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

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )

      // Initialize A

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = ZERO
         } // 10
      } // 20
      for (I = 1; I <= NMAX; I++) { // 30
         A( I, I ) = ONE
      } // 30
      OK = true;
      NT = 0

      if ( LSAMEN( 2, C2, 'EV' ) ) {

         // Test DGEEV

         SRNAMT = 'DGEEV '
         INFOT = 1
         dgeev('X', 'N', 0, A, 1, WR, WI, VL, 1, VR, 1, W, 1, INFO );
         chkxer('DGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         dgeev('N', 'X', 0, A, 1, WR, WI, VL, 1, VR, 1, W, 1, INFO );
         chkxer('DGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         dgeev('N', 'N', -1, A, 1, WR, WI, VL, 1, VR, 1, W, 1, INFO );
         chkxer('DGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 5
         dgeev('N', 'N', 2, A, 1, WR, WI, VL, 1, VR, 1, W, 6, INFO );
         chkxer('DGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 9
         dgeev('V', 'N', 2, A, 2, WR, WI, VL, 1, VR, 1, W, 8, INFO );
         chkxer('DGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 11
         dgeev('N', 'V', 2, A, 2, WR, WI, VL, 1, VR, 1, W, 8, INFO );
         chkxer('DGEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 13
         dgeev('V', 'V', 1, A, 1, WR, WI, VL, 1, VR, 1, W, 3, INFO );
         chkxer('DGEEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 7

      } else if ( LSAMEN( 2, C2, 'ES' ) ) {

         // Test DGEES

         SRNAMT = 'DGEES '
         INFOT = 1
         dgees('X', 'N', DSLECT, 0, A, 1, SDIM, WR, WI, VL, 1, W, 1, B, INFO );
         chkxer('DGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         dgees('N', 'X', DSLECT, 0, A, 1, SDIM, WR, WI, VL, 1, W, 1, B, INFO );
         chkxer('DGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 4
         dgees('N', 'S', DSLECT, -1, A, 1, SDIM, WR, WI, VL, 1, W, 1, B, INFO );
         chkxer('DGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 6
         dgees('N', 'S', DSLECT, 2, A, 1, SDIM, WR, WI, VL, 1, W, 6, B, INFO );
         chkxer('DGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 11
         dgees('V', 'S', DSLECT, 2, A, 2, SDIM, WR, WI, VL, 1, W, 6, B, INFO );
         chkxer('DGEES ', INFOT, NOUT, LERR, OK );
         INFOT = 13
         dgees('N', 'S', DSLECT, 1, A, 1, SDIM, WR, WI, VL, 1, W, 2, B, INFO );
         chkxer('DGEES ', INFOT, NOUT, LERR, OK );
         NT = NT + 6

      } else if ( LSAMEN( 2, C2, 'VX' ) ) {

         // Test DGEEVX

         SRNAMT = 'DGEEVX'
         INFOT = 1
         dgeevx('X', 'N', 'N', 'N', 0, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('DGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         dgeevx('N', 'X', 'N', 'N', 0, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('DGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         dgeevx('N', 'N', 'X', 'N', 0, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('DGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         dgeevx('N', 'N', 'N', 'X', 0, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('DGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 5
         dgeevx('N', 'N', 'N', 'N', -1, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('DGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7
         dgeevx('N', 'N', 'N', 'N', 2, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('DGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         dgeevx('N', 'V', 'N', 'N', 2, A, 2, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 6, IW, INFO );
         chkxer('DGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 13
         dgeevx('N', 'N', 'V', 'N', 2, A, 2, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 6, IW, INFO );
         chkxer('DGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 21
         dgeevx('N', 'N', 'N', 'N', 1, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 1, IW, INFO );
         chkxer('DGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 21
         dgeevx('N', 'V', 'N', 'N', 1, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 2, IW, INFO );
         chkxer('DGEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 21
         dgeevx('N', 'N', 'V', 'V', 1, A, 1, WR, WI, VL, 1, VR, 1, ILO, IHI, S, ABNRM, R1, R2, W, 3, IW, INFO );
         chkxer('DGEEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 11

      } else if ( LSAMEN( 2, C2, 'SX' ) ) {

         // Test DGEESX

         SRNAMT = 'DGEESX'
         INFOT = 1
         dgeesx('X', 'N', DSLECT, 'N', 0, A, 1, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 1, IW, 1, B, INFO );
         chkxer('DGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         dgeesx('N', 'X', DSLECT, 'N', 0, A, 1, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 1, IW, 1, B, INFO );
         chkxer('DGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         dgeesx('N', 'N', DSLECT, 'X', 0, A, 1, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 1, IW, 1, B, INFO );
         chkxer('DGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 5
         dgeesx('N', 'N', DSLECT, 'N', -1, A, 1, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 1, IW, 1, B, INFO );
         chkxer('DGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 7
         dgeesx('N', 'N', DSLECT, 'N', 2, A, 1, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 6, IW, 1, B, INFO );
         chkxer('DGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 12
         dgeesx('V', 'N', DSLECT, 'N', 2, A, 2, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 6, IW, 1, B, INFO );
         chkxer('DGEESX', INFOT, NOUT, LERR, OK );
         INFOT = 16
         dgeesx('N', 'N', DSLECT, 'N', 1, A, 1, SDIM, WR, WI, VL, 1, R1( 1 ), R2( 1 ), W, 2, IW, 1, B, INFO );
         chkxer('DGEESX', INFOT, NOUT, LERR, OK );
         NT = NT + 7

      } else if ( LSAMEN( 2, C2, 'BD' ) ) {

         // Test DGESVD

         SRNAMT = 'DGESVD'
         INFOT = 1
         dgesvd('X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, INFO );
         chkxer('DGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 2
         dgesvd('N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, INFO );
         chkxer('DGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 2
         dgesvd('O', 'O', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, INFO );
         chkxer('DGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 3
         dgesvd('N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, INFO );
         chkxer('DGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 4
         dgesvd('N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, INFO );
         chkxer('DGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 6
         dgesvd('N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, INFO );
         chkxer('DGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 9
         dgesvd('A', 'N', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, INFO );
         chkxer('DGESVD', INFOT, NOUT, LERR, OK );
         INFOT = 11
         dgesvd('N', 'A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, INFO );
         chkxer('DGESVD', INFOT, NOUT, LERR, OK );
         NT = 8
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT
         } else {
            WRITE( NOUT, FMT = 9998 )
         }

         // Test DGESDD

         SRNAMT = 'DGESDD'
         INFOT = 1
         dgesdd('X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 2
         dgesdd('N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 3
         dgesdd('N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 5
         dgesdd('N', 2, 1, A, 1, S, U, 1, VT, 1, W, 5, IW, INFO );
         chkxer('DGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 8
         dgesdd('A', 2, 1, A, 2, S, U, 1, VT, 1, W, 5, IW, INFO );
         chkxer('DGESDD', INFOT, NOUT, LERR, OK );
         INFOT = 10
         dgesdd('A', 1, 2, A, 1, S, U, 1, VT, 1, W, 5, IW, INFO );
         chkxer('DGESDD', INFOT, NOUT, LERR, OK );
         NT = 6
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT
         } else {
            WRITE( NOUT, FMT = 9998 )
         }

         // Test DGEJSV

         SRNAMT = 'DGEJSV'
         INFOT = 1
         dgejsv('X', 'U', 'V', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('DGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 2
         dgejsv('G', 'X', 'V', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('DGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 3
         dgejsv('G', 'U', 'X', 'R', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('DGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 4
         dgejsv('G', 'U', 'V', 'X', 'N', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('DGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 5
         dgejsv('G', 'U', 'V', 'R', 'X', 'N', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('DGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 6
         dgejsv('G', 'U', 'V', 'R', 'N', 'X', 0, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('DGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 7
         dgejsv('G', 'U', 'V', 'R', 'N', 'N', -1, 0, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('DGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 8
         dgejsv('G', 'U', 'V', 'R', 'N', 'N', 0, -1, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('DGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 10
         dgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 1, A, 1, S, U, 1, VT, 1, W, 1, IW, INFO);
         chkxer('DGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 13
         dgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 2, A, 2, S, U, 1, VT, 2, W, 1, IW, INFO);
         chkxer('DGEJSV', INFOT, NOUT, LERR, OK );
         INFOT = 15
         dgejsv('G', 'U', 'V', 'R', 'N', 'N', 2, 2, A, 2, S, U, 2, VT, 1, W, 1, IW, INFO);
         chkxer('DGEJSV', INFOT, NOUT, LERR, OK );
         NT = 11
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT
         } else {
            WRITE( NOUT, FMT = 9998 )
         }

         // Test DGESVDX

         SRNAMT = 'DGESVDX'
         INFOT = 1
         dgesvdx('X', 'N', 'A', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         dgesvdx('N', 'X', 'A', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         dgesvdx('N', 'N', 'X', 0, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         dgesvdx('N', 'N', 'A', -1, 0, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 5
         dgesvdx('N', 'N', 'A', 0, -1, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 7
         dgesvdx('N', 'N', 'A', 2, 1, A, 1, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         dgesvdx('N', 'N', 'V', 2, 1, A, 2, -ONE, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         dgesvdx('N', 'N', 'V', 2, 1, A, 2, ONE, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         dgesvdx('N', 'N', 'I', 2, 2, A, 2, ZERO, ZERO, 0, 1, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 11
         dgesvdx('V', 'N', 'I', 2, 2, A, 2, ZERO, ZERO, 1, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 15
         dgesvdx('V', 'N', 'A', 2, 2, A, 2, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESVDX', INFOT, NOUT, LERR, OK );
         INFOT = 17
         dgesvdx('N', 'V', 'A', 2, 2, A, 2, ZERO, ZERO, 0, 0, NS, S, U, 1, VT, 1, W, 1, IW, INFO );
         chkxer('DGESVDX', INFOT, NOUT, LERR, OK );
         NT = 12
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT
         } else {
            WRITE( NOUT, FMT = 9998 )
         }

         // Test DGESVDQ

         SRNAMT = 'DGESVDQ'
         INFOT = 1
         dgesvdq('X', 'P', 'T', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('DGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         dgesvdq('A', 'X', 'T', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('DGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         dgesvdq('A', 'P', 'X', 'A', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('DGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 4
         dgesvdq('A', 'P', 'T', 'X', 'A', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('DGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 5
         dgesvdq('A', 'P', 'T', 'A', 'X', 0, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('DGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 6
         dgesvdq('A', 'P', 'T', 'A', 'A', -1, 0, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('DGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 7
         dgesvdq('A', 'P', 'T', 'A', 'A', 0, 1, A, 1, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('DGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 9
         dgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 0, S, U, 0, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('DGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 12
         dgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, -1, VT, 0, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('DGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 14
         dgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, 1, VT, -1, NS, IW, 1, W, 1, W, 1, INFO );
         chkxer('DGESVDQ', INFOT, NOUT, LERR, OK );
         INFOT = 17
         dgesvdq('A', 'P', 'T', 'A', 'A', 1, 1, A, 1, S, U, 1, VT, 1, NS, IW, -5, W, 1, W, 1, INFO );
         chkxer('DGESVDQ', INFOT, NOUT, LERR, OK );
         NT = 11
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT
         } else {
            WRITE( NOUT, FMT = 9998 )
         }
      }

      // Print a summary line.

      if ( .NOT.LSAMEN( 2, C2, 'BD' ) ) {
         if ( OK ) {
            WRITE( NOUT, FMT = 9999 )SRNAMT( 1:LEN_TRIM( SRNAMT ) ), NT
         } else {
            WRITE( NOUT, FMT = 9998 )
         }
      }

 9999 FORMAT( 1X, A, ' passed the tests of the error exits (', I3, ' tests done)' )
 9998 FORMAT( ' *** ', A, ' failed the tests of the error exits ***' )
      RETURN

      // End of DERRED
      }
