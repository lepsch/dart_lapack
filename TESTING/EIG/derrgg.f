      SUBROUTINE DERRGG( PATH, NUNIT );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NMAX, LW;
      const              NMAX = 3, LW = 6*NMAX ;
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                DUMMYK, DUMMYL, I, IFST, ILO, IHI, ILST, INFO, J, M, NCYCLE, NT, SDIM, LWORK;
      double             ANRM, BNRM, DIF, SCALE, TOLA, TOLB;
      // ..
      // .. Local Arrays ..
      bool               BW( NMAX ), SEL( NMAX );
      int                IW( NMAX ), IDUM(NMAX);
      double             A( NMAX, NMAX ), B( NMAX, NMAX ), LS( NMAX ), Q( NMAX, NMAX ), R1( NMAX ), R2( NMAX ), R3( NMAX ), RCE( 2 ), RCV( 2 ), RS( NMAX ), TAU( NMAX ), U( NMAX, NMAX ), V( NMAX, NMAX ), W( LW ), Z( NMAX, NMAX );
      // ..
      // .. External Functions ..
      bool               DLCTES, DLCTSX, LSAMEN;
      // EXTERNAL DLCTES, DLCTSX, LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, DGGES, DGGESX, DGGEV, DGGEVX, DGGGLM, DGGHRD, DGGLSE, DGGQRF, DGGRQF, DHGEQZ, DORCSD, DTGEVC, DTGEXC, DTGSEN, DTGSJA, DTGSNA, DTGSYL, DGGHD3, DGGES3, DGGEV3, DGGSVD3, DGGSVP3, XLAENV
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

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         SEL( J ) = true;
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = ZERO;
            B( I, J ) = ZERO;
         } // 10
      } // 20
      for (I = 1; I <= NMAX; I++) { // 30
         A( I, I ) = ONE;
         B( I, I ) = ONE;
      } // 30
      OK = true;
      TOLA = 1.0;
      TOLB = 1.0;
      IFST = 1;
      ILST = 1;
      NT = 0;
      LWORK = 1;

      // Call XLAENV to set the parameters used in CLAQZ0

      xlaenv(12, 10 );
      xlaenv(13, 12 );
      xlaenv(14, 13 );
      xlaenv(15, 2 );
      xlaenv(17, 10 );

      // Test error exits for the GG path.

      if ( LSAMEN( 2, C2, 'GG' ) ) {

         // DGGHRD

         SRNAMT = 'DGGHRD';
         INFOT = 1;
         dgghrd('/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('DGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgghrd('N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('DGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgghrd('N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('DGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgghrd('N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('DGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgghrd('N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('DGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dgghrd('N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, INFO );
         chkxer('DGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dgghrd('N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, INFO );
         chkxer('DGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dgghrd('V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO );
         chkxer('DGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dgghrd('N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO );
         chkxer('DGGHRD', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // DGGHD3

         SRNAMT = 'DGGHD3';
         INFOT = 1;
         dgghd3('/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgghd3('N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgghd3('N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dgghd3('N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgghd3('N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dgghd3('N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dgghd3('N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dgghd3('V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dgghd3('N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DGGHD3', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // DHGEQZ

         SRNAMT = 'DHGEQZ';
         INFOT = 1;
         dhgeqz('/', 'N', 'N', 0, 1, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dhgeqz('E', '/', 'N', 0, 1, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dhgeqz('E', 'N', '/', 0, 1, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dhgeqz('E', 'N', 'N', -1, 0, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dhgeqz('E', 'N', 'N', 0, 0, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dhgeqz('E', 'N', 'N', 0, 1, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dhgeqz('E', 'N', 'N', 2, 1, 1, A, 1, B, 2, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dhgeqz('E', 'N', 'N', 2, 1, 1, A, 2, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dhgeqz('E', 'V', 'N', 2, 1, 1, A, 2, B, 2, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         dhgeqz('E', 'N', 'V', 2, 1, 1, A, 2, B, 2, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('DHGEQZ', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // DTGEVC

         SRNAMT = 'DTGEVC';
         INFOT = 1;
         dtgevc('/', 'A', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, INFO );
         chkxer('DTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dtgevc('R', '/', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, INFO );
         chkxer('DTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dtgevc('R', 'A', SEL, -1, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, INFO );
         chkxer('DTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dtgevc('R', 'A', SEL, 2, A, 1, B, 2, Q, 1, Z, 2, 0, M, W, INFO );
         chkxer('DTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dtgevc('R', 'A', SEL, 2, A, 2, B, 1, Q, 1, Z, 2, 0, M, W, INFO );
         chkxer('DTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dtgevc('L', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W, INFO );
         chkxer('DTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dtgevc('R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W, INFO );
         chkxer('DTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dtgevc('R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 2, 1, M, W, INFO );
         chkxer('DTGEVC', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the GSV path.

      } else if ( LSAMEN( 3, PATH, 'GSV' ) ) {

         // DGGSVD3

         SRNAMT = 'DGGSVD3';
         INFOT = 1;
         dggsvd3('/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('DGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dggsvd3('N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('DGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dggsvd3('N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('DGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dggsvd3('N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('DGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dggsvd3('N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('DGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dggsvd3('N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('DGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dggsvd3('N', 'N', 'N', 2, 1, 1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('DGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dggsvd3('N', 'N', 'N', 1, 1, 2, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('DGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dggsvd3('U', 'N', 'N', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B, 2, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('DGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         dggsvd3('N', 'V', 'N', 1, 1, 2, DUMMYK, DUMMYL, A, 1, B, 2, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('DGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         dggsvd3('N', 'N', 'Q', 1, 2, 1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('DGGSVD3', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // DGGSVP3

         SRNAMT = 'DGGSVP3';
         INFOT = 1;
         dggsvp3('/', 'N', 'N', 0, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('DGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dggsvp3('N', '/', 'N', 0, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('DGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dggsvp3('N', 'N', '/', 0, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('DGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dggsvp3('N', 'N', 'N', -1, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('DGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dggsvp3('N', 'N', 'N', 0, -1, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('DGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dggsvp3('N', 'N', 'N', 0, 0, -1, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('DGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dggsvp3('N', 'N', 'N', 2, 1, 1, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('DGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dggsvp3('N', 'N', 'N', 1, 2, 1, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('DGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dggsvp3('U', 'N', 'N', 2, 2, 2, A, 2, B, 2, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('DGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         dggsvp3('N', 'V', 'N', 1, 2, 1, A, 1, B, 2, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('DGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         dggsvp3('N', 'N', 'Q', 1, 1, 2, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('DGGSVP3', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // DTGSJA

         SRNAMT = 'DTGSJA';
         INFOT = 1;
         dtgsja('/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('DTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dtgsja('N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('DTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dtgsja('N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('DTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dtgsja('N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('DTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dtgsja('N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('DTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dtgsja('N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('DTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dtgsja('N', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 0, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('DTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dtgsja('N', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 0, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('DTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         dtgsja('U', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 0, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('DTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         dtgsja('N', 'V', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 0, Q, 1, W, NCYCLE, INFO );
         chkxer('DTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         dtgsja('N', 'N', 'Q', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 0, W, NCYCLE, INFO );
         chkxer('DTGSJA', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

      // Test error exits for the GLM path.

      } else if ( LSAMEN( 3, PATH, 'GLM' ) ) {

         // DGGGLM

         SRNAMT = 'DGGGLM';
         INFOT = 1;
         dggglm(-1, 0, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('DGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dggglm(0, -1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('DGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dggglm(0, 1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('DGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dggglm(0, 0, -1, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('DGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dggglm(1, 0, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('DGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dggglm(0, 0, 0, A, 0, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('DGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dggglm(0, 0, 0, A, 1, B, 0, R1, R2, R3, W, LW, INFO );
         chkxer('DGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dggglm(1, 1, 1, A, 1, B, 1, R1, R2, R3, W, 1, INFO );
         chkxer('DGGGLM', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the LSE path.

      } else if ( LSAMEN( 3, PATH, 'LSE' ) ) {

         // DGGLSE

         SRNAMT = 'DGGLSE';
         INFOT = 1;
         dgglse(-1, 0, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('DGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgglse(0, -1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('DGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgglse(0, 0, -1, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('DGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgglse(0, 0, 1, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('DGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgglse(0, 1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('DGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgglse(0, 0, 0, A, 0, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('DGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dgglse(0, 0, 0, A, 1, B, 0, R1, R2, R3, W, LW, INFO );
         chkxer('DGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dgglse(1, 1, 1, A, 1, B, 1, R1, R2, R3, W, 1, INFO );
         chkxer('DGGLSE', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the CSD path.

      } else if ( LSAMEN( 3, PATH, 'CSD' ) ) {

         // DORCSD

         SRNAMT = 'DORCSD';
         INFOT = 7;
         dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', -1, 0, 0, A, 1, A, 1, A, 1, A, 1, A, A, 1, A, 1, A, 1, A, 1, W, LW, IW, INFO );
         chkxer('DORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, -1, 0, A, 1, A, 1, A, 1, A, 1, A, A, 1, A, 1, A, 1, A, 1, W, LW, IW, INFO );
         chkxer('DORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, -1, A, 1, A, 1, A, 1, A, 1, A, A, 1, A, 1, A, 1, A, 1, W, LW, IW, INFO );
         chkxer('DORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, -1, A, 1, A, 1, A, 1, A, A, 1, A, 1, A, 1, A, 1, W, LW, IW, INFO );
         chkxer('DORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, A, A, -1, A, 1, A, 1, A, 1, W, LW, IW, INFO );
         chkxer('DORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, A, A, 1, A, -1, A, 1, A, 1, W, LW, IW, INFO );
         chkxer('DORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 24;
         dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, A, A, 1, A, 1, A, -1, A, 1, W, LW, IW, INFO );
         chkxer('DORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 26;
         dorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, A, A, 1, A, 1, A, 1, A, -1, W, LW, IW, INFO );
         chkxer('DORCSD', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the GQR path.

      } else if ( LSAMEN( 3, PATH, 'GQR' ) ) {

         // DGGQRF

         SRNAMT = 'DGGQRF';
         INFOT = 1;
         dggqrf(-1, 0, 0, A, 1, R1, B, 1, R2, W, LW, INFO );
         chkxer('DGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dggqrf(0, -1, 0, A, 1, R1, B, 1, R2, W, LW, INFO );
         chkxer('DGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dggqrf(0, 0, -1, A, 1, R1, B, 1, R2, W, LW, INFO );
         chkxer('DGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dggqrf(0, 0, 0, A, 0, R1, B, 1, R2, W, LW, INFO );
         chkxer('DGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dggqrf(0, 0, 0, A, 1, R1, B, 0, R2, W, LW, INFO );
         chkxer('DGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dggqrf(1, 1, 2, A, 1, R1, B, 1, R2, W, 1, INFO );
         chkxer('DGGQRF', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // DGGRQF

         SRNAMT = 'DGGRQF';
         INFOT = 1;
         dggrqf(-1, 0, 0, A, 1, R1, B, 1, R2, W, LW, INFO );
         chkxer('DGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dggrqf(0, -1, 0, A, 1, R1, B, 1, R2, W, LW, INFO );
         chkxer('DGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dggrqf(0, 0, -1, A, 1, R1, B, 1, R2, W, LW, INFO );
         chkxer('DGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dggrqf(0, 0, 0, A, 0, R1, B, 1, R2, W, LW, INFO );
         chkxer('DGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dggrqf(0, 0, 0, A, 1, R1, B, 0, R2, W, LW, INFO );
         chkxer('DGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dggrqf(1, 1, 2, A, 1, R1, B, 1, R2, W, 1, INFO );
         chkxer('DGGRQF', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

      // Test error exits for the DGS, DGV, DGX, and DXV paths.

      } else if ( LSAMEN( 3, PATH, 'DGS' ) || LSAMEN( 3, PATH, 'DGV' ) || LSAMEN( 3, PATH, 'DGX' ) || LSAMEN( 3, PATH, 'DXV' ) ) {

         // DGGES

         SRNAMT = 'DGGES ';
         INFOT = 1;
         dgges('/', 'N', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('DGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgges('N', '/', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('DGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgges('N', 'V', '/', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('DGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgges('N', 'V', 'S', DLCTES, -1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('DGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dgges('N', 'V', 'S', DLCTES, 1, A, 0, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('DGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dgges('N', 'V', 'S', DLCTES, 1, A, 1, B, 0, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('DGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dgges('N', 'V', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 0, U, 1, W, 1, BW, INFO );
         chkxer('DGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dgges('V', 'V', 'S', DLCTES, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 1, U, 2, W, 1, BW, INFO );
         chkxer('DGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         dgges('N', 'V', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 0, W, 1, BW, INFO );
         chkxer('DGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         dgges('V', 'V', 'S', DLCTES, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 1, W, 1, BW, INFO );
         chkxer('DGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 19;
         dgges('V', 'V', 'S', DLCTES, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 2, W, 1, BW, INFO );
         chkxer('DGGES ', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // DGGES3

         SRNAMT = 'DGGES3 ';
         INFOT = 1;
         dgges3('/', 'N', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('DGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dgges3('N', '/', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('DGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dgges3('N', 'V', '/', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('DGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dgges3('N', 'V', 'S', DLCTES, -1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('DGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dgges3('N', 'V', 'S', DLCTES, 1, A, 0, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('DGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dgges3('N', 'V', 'S', DLCTES, 1, A, 1, B, 0, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('DGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dgges3('N', 'V', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 0, U, 1, W, 1, BW, INFO );
         chkxer('DGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dgges3('V', 'V', 'S', DLCTES, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 1, U, 2, W, 1, BW, INFO );
         chkxer('DGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         dgges3('N', 'V', 'S', DLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 0, W, 1, BW, INFO );
         chkxer('DGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         dgges3('V', 'V', 'S', DLCTES, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 1, W, 1, BW, INFO );
         chkxer('DGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 19;
         dgges3('V', 'V', 'S', DLCTES, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 2, W, 1, BW, INFO );
         chkxer('DGGES3 ', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // DGGESX

         SRNAMT = 'DGGESX';
         INFOT = 1;
         dggesx('/', 'N', 'S', DLCTSX, 'N', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dggesx('N', '/', 'S', DLCTSX, 'N', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dggesx('V', 'V', '/', DLCTSX, 'N', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dggesx('V', 'V', 'S', DLCTSX, '/', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dggesx('V', 'V', 'S', DLCTSX, 'B', -1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dggesx('V', 'V', 'S', DLCTSX, 'B', 1, A, 0, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dggesx('V', 'V', 'S', DLCTSX, 'B', 1, A, 1, B, 0, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dggesx('V', 'V', 'S', DLCTSX, 'B', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 0, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dggesx('V', 'V', 'S', DLCTSX, 'B', 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         dggesx('V', 'V', 'S', DLCTSX, 'B', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 0, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         dggesx('V', 'V', 'S', DLCTSX, 'B', 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         dggesx('V', 'V', 'S', DLCTSX, 'B', 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 2, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 24;
         dggesx('V', 'V', 'S', DLCTSX, 'V', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 32, IW, 0, BW, INFO );
         chkxer('DGGESX', INFOT, NOUT, LERR, OK );
         NT = NT + 13;

         // DGGEV

         SRNAMT = 'DGGEV ';
         INFOT = 1;
         dggev('/', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('DGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dggev('N', '/', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('DGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dggev('V', 'V', -1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('DGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dggev('V', 'V', 1, A, 0, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('DGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dggev('V', 'V', 1, A, 1, B, 0, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('DGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dggev('N', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 0, U, 1, W, 1, INFO );
         chkxer('DGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dggev('V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 1, U, 2, W, 1, INFO );
         chkxer('DGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dggev('V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 0, W, 1, INFO );
         chkxer('DGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dggev('V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 1, W, 1, INFO );
         chkxer('DGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dggev('V', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('DGGEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // DGGEV3

         xlaenv(12, 20 );
         xlaenv(13, 4 );
         xlaenv(14, 13 );
         xlaenv(15, 2 );
         xlaenv(17, 10 );
         SRNAMT = 'DGGEV3 ';
         INFOT = 1;
         dggev3('/', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('DGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dggev3('N', '/', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('DGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dggev3('V', 'V', -1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('DGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dggev3('V', 'V', 1, A, 0, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('DGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dggev3('V', 'V', 1, A, 1, B, 0, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('DGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dggev3('N', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 0, U, 1, W, 1, INFO );
         chkxer('DGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dggev3('V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 1, U, 2, W, 1, INFO );
         chkxer('DGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dggev3('V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 0, W, 1, INFO );
         chkxer('DGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dggev3('V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 1, W, 1, INFO );
         chkxer('DGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dggev3('V', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('DGGEV3 ', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // DGGEVX

         SRNAMT = 'DGGEVX';
         INFOT = 1;
         dggevx('/', 'N', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('DGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dggevx('N', '/', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('DGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dggevx('N', 'N', '/', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('DGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dggevx('N', 'N', 'N', '/', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('DGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dggevx('N', 'N', 'N', 'N', -1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('DGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dggevx('N', 'N', 'N', 'N', 1, A, 0, B, 1, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('DGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dggevx('N', 'N', 'N', 'N', 1, A, 1, B, 0, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('DGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dggevx('N', 'N', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 0, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('DGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dggevx('N', 'V', 'N', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 1, U, 2, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('DGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dggevx('N', 'N', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 0, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('DGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dggevx('N', 'N', 'V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('DGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 26;
         dggevx('N', 'N', 'V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 2, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('DGGEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // DTGEXC

         SRNAMT = 'DTGEXC';
         INFOT = 3;
         dtgexc( true , true , -1, A, 1, B, 1, Q, 1, Z, 1, IFST, ILST, W, 1, INFO );
         chkxer('DTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dtgexc( true , true , 1, A, 0, B, 1, Q, 1, Z, 1, IFST, ILST, W, 1, INFO );
         chkxer('DTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dtgexc( true , true , 1, A, 1, B, 0, Q, 1, Z, 1, IFST, ILST, W, 1, INFO );
         chkxer('DTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dtgexc( false , true , 1, A, 1, B, 1, Q, 0, Z, 1, IFST, ILST, W, 1, INFO );
         chkxer('DTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dtgexc( true , true , 1, A, 1, B, 1, Q, 0, Z, 1, IFST, ILST, W, 1, INFO );
         chkxer('DTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dtgexc( true , false , 1, A, 1, B, 1, Q, 1, Z, 0, IFST, ILST, W, 1, INFO );
         chkxer('DTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dtgexc( true , true , 1, A, 1, B, 1, Q, 1, Z, 0, IFST, ILST, W, 1, INFO );
         chkxer('DTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dtgexc( true , true , 1, A, 1, B, 1, Q, 1, Z, 1, IFST, ILST, W, 0, INFO );
         chkxer('DTGEXC', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

         // DTGSEN

         SRNAMT = 'DTGSEN';
         INFOT = 1;
         dtgsen(-1, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('DTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dtgsen(1, true , true , SEL, -1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('DTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dtgsen(1, true , true , SEL, 1, A, 0, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('DTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dtgsen(1, true , true , SEL, 1, A, 1, B, 0, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('DTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dtgsen(1, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 0, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('DTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dtgsen(1, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 0, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('DTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         dtgsen(0, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('DTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         dtgsen(1, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('DTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         dtgsen(2, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('DTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 24;
         dtgsen(0, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 0, INFO );
         chkxer('DTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 24;
         dtgsen(1, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 0, INFO );
         chkxer('DTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 24;
         dtgsen(2, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 1, INFO );
         chkxer('DTGSEN', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // DTGSNA

         SRNAMT = 'DTGSNA';
         INFOT = 1;
         dtgsna('/', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('DTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dtgsna('B', '/', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('DTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dtgsna('B', 'A', SEL, -1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('DTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dtgsna('B', 'A', SEL, 1, A, 0, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('DTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dtgsna('B', 'A', SEL, 1, A, 1, B, 0, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('DTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dtgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 0, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('DTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dtgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 0, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('DTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dtgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 0, M, W, 1, IW, INFO );
         chkxer('DTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         dtgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 0, IW, INFO );
         chkxer('DTGSNA', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // DTGSYL

         SRNAMT = 'DTGSYL';
         INFOT = 1;
         dtgsyl('/', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('DTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dtgsyl('N', -1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('DTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dtgsyl('N', 0, 0, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('DTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dtgsyl('N', 0, 1, 0, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('DTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dtgsyl('N', 0, 1, 1, A, 0, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('DTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dtgsyl('N', 0, 1, 1, A, 1, B, 0, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('DTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dtgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 0, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('DTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dtgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 0, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('DTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dtgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 0, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('DTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         dtgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 0, SCALE, DIF, W, 1, IW, INFO );
         chkxer('DTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         dtgsyl('N', 1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('DTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         dtgsyl('N', 2, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('DTGSYL', INFOT, NOUT, LERR, OK );
         NT = NT + 12;
      }

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH, NT;
      } else {
         WRITE( NOUT, FMT = 9998 )PATH;
      }

 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits (', I3, ' tests done)' );
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ', 'exits ***' );

      RETURN;

      // End of DERRGG

      }
