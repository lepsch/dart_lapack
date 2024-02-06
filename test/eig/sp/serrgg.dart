      void serrgg(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                NUNIT;
      // ..

      int                NMAX, LW;
      const              NMAX = 3, LW = 6*NMAX ;
      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      String             C2;
      int                DUMMYK, DUMMYL, I, IFST, ILO, IHI, ILST, INFO, J, M, NCYCLE, NT, SDIM, LWORK;
      double               ANRM, BNRM, DIF, SCALE, TOLA, TOLB;
      bool               BW( NMAX ), SEL( NMAX );
      int                IW( NMAX ), IDUM(NMAX);
      double               A( NMAX, NMAX ), B( NMAX, NMAX ), LS( NMAX ), Q( NMAX, NMAX ), R1( NMAX ), R2( NMAX ), R3( NMAX ), RCE( 2 ), RCV( 2 ), RS( NMAX ), TAU( NMAX ), U( NMAX, NMAX ), V( NMAX, NMAX ), W( LW ), Z( NMAX, NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN, SLCTES, SLCTSX;
      // EXTERNAL LSAMEN, SLCTES, SLCTSX
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, SGGES, SGGESX, SGGEV, SGGEVX, SGGGLM, SGGHRD, SGGLSE, SGGQRF, SGGRQF, SHGEQZ, SORCSD, STGEVC, STGEXC, STGSEN, STGSJA, STGSNA, STGSYL, SGGES3, SGGEV3, SGGHD3, SGGSVD3, SGGSVP3, XLAENV
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

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         SEL[J] = true;
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = ZERO;
            B[I][J] = ZERO;
         } // 10
      } // 20
      for (I = 1; I <= NMAX; I++) { // 30
         A[I][I] = ONE;
         B[I][I] = ONE;
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

      if ( lsamen( 2, C2, 'GG' ) ) {

         // SGGHRD

        srnamc.SRNAMT = 'SGGHRD';
         INFOT = 1;
         sgghrd('/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('SGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgghrd('N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('SGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgghrd('N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('SGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgghrd('N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('SGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgghrd('N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('SGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgghrd('N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, INFO );
         chkxer('SGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sgghrd('N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, INFO );
         chkxer('SGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sgghrd('V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO );
         chkxer('SGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         sgghrd('N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO );
         chkxer('SGGHRD', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // SGGHD3

        srnamc.SRNAMT = 'SGGHD3';
         INFOT = 1;
         sgghd3('/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgghd3('N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgghd3('N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sgghd3('N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgghd3('N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgghd3('N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sgghd3('N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sgghd3('V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         sgghd3('N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SGGHD3', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // SHGEQZ

        srnamc.SRNAMT = 'SHGEQZ';
         INFOT = 1;
         shgeqz('/', 'N', 'N', 0, 1, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         shgeqz('E', '/', 'N', 0, 1, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         shgeqz('E', 'N', '/', 0, 1, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         shgeqz('E', 'N', 'N', -1, 0, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         shgeqz('E', 'N', 'N', 0, 0, 0, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         shgeqz('E', 'N', 'N', 0, 1, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         shgeqz('E', 'N', 'N', 2, 1, 1, A, 1, B, 2, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         shgeqz('E', 'N', 'N', 2, 1, 1, A, 2, B, 1, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         shgeqz('E', 'V', 'N', 2, 1, 1, A, 2, B, 2, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         shgeqz('E', 'N', 'V', 2, 1, 1, A, 2, B, 2, R1, R2, R3, Q, 1, Z, 1, W, LW, INFO );
         chkxer('SHGEQZ', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // STGEVC

        srnamc.SRNAMT = 'STGEVC';
         INFOT = 1;
         stgevc('/', 'A', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, INFO );
         chkxer('STGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         stgevc('R', '/', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, INFO );
         chkxer('STGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         stgevc('R', 'A', SEL, -1, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, INFO );
         chkxer('STGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         stgevc('R', 'A', SEL, 2, A, 1, B, 2, Q, 1, Z, 2, 0, M, W, INFO );
         chkxer('STGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         stgevc('R', 'A', SEL, 2, A, 2, B, 1, Q, 1, Z, 2, 0, M, W, INFO );
         chkxer('STGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         stgevc('L', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W, INFO );
         chkxer('STGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         stgevc('R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W, INFO );
         chkxer('STGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         stgevc('R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 2, 1, M, W, INFO );
         chkxer('STGEVC', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the GSV path.

      } else if ( lsamen( 3, PATH, 'GSV' ) ) {

         // SGGSVD3

        srnamc.SRNAMT = 'SGGSVD3';
         INFOT = 1;
         sggsvd3('/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('SGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sggsvd3('N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('SGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sggsvd3('N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('SGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sggsvd3('N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('SGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sggsvd3('N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('SGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sggsvd3('N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('SGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sggsvd3('N', 'N', 'N', 2, 1, 1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('SGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sggsvd3('N', 'N', 'N', 1, 1, 2, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('SGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         sggsvd3('U', 'N', 'N', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B, 2, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('SGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         sggsvd3('N', 'V', 'N', 1, 1, 2, DUMMYK, DUMMYL, A, 1, B, 2, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('SGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         sggsvd3('N', 'N', 'Q', 1, 2, 1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, IDUM, INFO );
         chkxer('SGGSVD3', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // SGGSVP3

        srnamc.SRNAMT = 'SGGSVP3';
         INFOT = 1;
         sggsvp3('/', 'N', 'N', 0, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('SGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sggsvp3('N', '/', 'N', 0, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('SGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sggsvp3('N', 'N', '/', 0, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('SGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sggsvp3('N', 'N', 'N', -1, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('SGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sggsvp3('N', 'N', 'N', 0, -1, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('SGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sggsvp3('N', 'N', 'N', 0, 0, -1, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('SGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sggsvp3('N', 'N', 'N', 2, 1, 1, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('SGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sggsvp3('N', 'N', 'N', 1, 2, 1, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('SGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         sggsvp3('U', 'N', 'N', 2, 2, 2, A, 2, B, 2, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('SGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         sggsvp3('N', 'V', 'N', 1, 2, 1, A, 1, B, 2, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('SGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         sggsvp3('N', 'N', 'Q', 1, 1, 2, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, TAU, W, LWORK, INFO );
         chkxer('SGGSVP3', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // STGSJA

        srnamc.SRNAMT = 'STGSJA';
         INFOT = 1;
         stgsja('/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('STGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         stgsja('N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('STGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         stgsja('N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('STGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         stgsja('N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('STGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         stgsja('N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('STGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         stgsja('N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('STGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         stgsja('N', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 0, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('STGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         stgsja('N', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 0, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('STGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         stgsja('U', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 0, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('STGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         stgsja('N', 'V', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 0, Q, 1, W, NCYCLE, INFO );
         chkxer('STGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         stgsja('N', 'N', 'Q', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 0, W, NCYCLE, INFO );
         chkxer('STGSJA', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

      // Test error exits for the GLM path.

      } else if ( lsamen( 3, PATH, 'GLM' ) ) {

         // SGGGLM

        srnamc.SRNAMT = 'SGGGLM';
         INFOT = 1;
         sggglm(-1, 0, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('SGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sggglm(0, -1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('SGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sggglm(0, 1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('SGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sggglm(0, 0, -1, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('SGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sggglm(1, 0, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('SGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sggglm(0, 0, 0, A, 0, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('SGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sggglm(0, 0, 0, A, 1, B, 0, R1, R2, R3, W, LW, INFO );
         chkxer('SGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sggglm(1, 1, 1, A, 1, B, 1, R1, R2, R3, W, 1, INFO );
         chkxer('SGGGLM', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the LSE path.

      } else if ( lsamen( 3, PATH, 'LSE' ) ) {

         // SGGLSE

        srnamc.SRNAMT = 'SGGLSE';
         INFOT = 1;
         sgglse(-1, 0, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('SGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgglse(0, -1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('SGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgglse(0, 0, -1, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('SGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgglse(0, 0, 1, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('SGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgglse(0, 1, 0, A, 1, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('SGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgglse(0, 0, 0, A, 0, B, 1, R1, R2, R3, W, LW, INFO );
         chkxer('SGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgglse(0, 0, 0, A, 1, B, 0, R1, R2, R3, W, LW, INFO );
         chkxer('SGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sgglse(1, 1, 1, A, 1, B, 1, R1, R2, R3, W, 1, INFO );
         chkxer('SGGLSE', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the CSD path.

      } else if ( lsamen( 3, PATH, 'CSD' ) ) {

         // SORCSD

        srnamc.SRNAMT = 'SORCSD';
         INFOT = 7;
         sorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', -1, 0, 0, A, 1, A, 1, A, 1, A, 1, A, A, 1, A, 1, A, 1, A, 1, W, LW, IW, INFO );
         chkxer('SORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, -1, 0, A, 1, A, 1, A, 1, A, 1, A, A, 1, A, 1, A, 1, A, 1, W, LW, IW, INFO );
         chkxer('SORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, -1, A, 1, A, 1, A, 1, A, 1, A, A, 1, A, 1, A, 1, A, 1, W, LW, IW, INFO );
         chkxer('SORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, -1, A, 1, A, 1, A, 1, A, A, 1, A, 1, A, 1, A, 1, W, LW, IW, INFO );
         chkxer('SORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         sorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, A, A, -1, A, 1, A, 1, A, 1, W, LW, IW, INFO );
         chkxer('SORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         sorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, A, A, 1, A, -1, A, 1, A, 1, W, LW, IW, INFO );
         chkxer('SORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 24;
         sorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, A, A, 1, A, 1, A, -1, A, 1, W, LW, IW, INFO );
         chkxer('SORCSD', INFOT, NOUT, LERR, OK );
         INFOT = 26;
         sorcsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, A, A, 1, A, 1, A, 1, A, -1, W, LW, IW, INFO );
         chkxer('SORCSD', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the GQR path.

      } else if ( lsamen( 3, PATH, 'GQR' ) ) {

         // SGGQRF

        srnamc.SRNAMT = 'SGGQRF';
         INFOT = 1;
         sggqrf(-1, 0, 0, A, 1, R1, B, 1, R2, W, LW, INFO );
         chkxer('SGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sggqrf(0, -1, 0, A, 1, R1, B, 1, R2, W, LW, INFO );
         chkxer('SGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sggqrf(0, 0, -1, A, 1, R1, B, 1, R2, W, LW, INFO );
         chkxer('SGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sggqrf(0, 0, 0, A, 0, R1, B, 1, R2, W, LW, INFO );
         chkxer('SGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sggqrf(0, 0, 0, A, 1, R1, B, 0, R2, W, LW, INFO );
         chkxer('SGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sggqrf(1, 1, 2, A, 1, R1, B, 1, R2, W, 1, INFO );
         chkxer('SGGQRF', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // SGGRQF

        srnamc.SRNAMT = 'SGGRQF';
         INFOT = 1;
         sggrqf(-1, 0, 0, A, 1, R1, B, 1, R2, W, LW, INFO );
         chkxer('SGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sggrqf(0, -1, 0, A, 1, R1, B, 1, R2, W, LW, INFO );
         chkxer('SGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sggrqf(0, 0, -1, A, 1, R1, B, 1, R2, W, LW, INFO );
         chkxer('SGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sggrqf(0, 0, 0, A, 0, R1, B, 1, R2, W, LW, INFO );
         chkxer('SGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sggrqf(0, 0, 0, A, 1, R1, B, 0, R2, W, LW, INFO );
         chkxer('SGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sggrqf(1, 1, 2, A, 1, R1, B, 1, R2, W, 1, INFO );
         chkxer('SGGRQF', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

      // Test error exits for the SGS, SGV, SGX, and SXV paths.

      } else if ( lsamen( 3, PATH, 'SGS' ) || lsamen( 3, PATH, 'SGV' ) || lsamen( 3, PATH, 'SGX' ) || lsamen( 3, PATH, 'SXV' ) ) {

         // SGGES

        srnamc.SRNAMT = 'SGGES ';
         INFOT = 1;
         sgges('/', 'N', 'S', SLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('SGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgges('N', '/', 'S', SLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('SGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgges('N', 'V', '/', SLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('SGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgges('N', 'V', 'S', SLCTES, -1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('SGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgges('N', 'V', 'S', SLCTES, 1, A, 0, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('SGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sgges('N', 'V', 'S', SLCTES, 1, A, 1, B, 0, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('SGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         sgges('N', 'V', 'S', SLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 0, U, 1, W, 1, BW, INFO );
         chkxer('SGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         sgges('V', 'V', 'S', SLCTES, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 1, U, 2, W, 1, BW, INFO );
         chkxer('SGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         sgges('N', 'V', 'S', SLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 0, W, 1, BW, INFO );
         chkxer('SGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         sgges('V', 'V', 'S', SLCTES, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 1, W, 1, BW, INFO );
         chkxer('SGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 19;
         sgges('V', 'V', 'S', SLCTES, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 2, W, 1, BW, INFO );
         chkxer('SGGES ', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // SGGES3

        srnamc.SRNAMT = 'SGGES3';
         INFOT = 1;
         sgges3('/', 'N', 'S', SLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('SGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sgges3('N', '/', 'S', SLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('SGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sgges3('N', 'V', '/', SLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('SGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sgges3('N', 'V', 'S', SLCTES, -1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('SGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sgges3('N', 'V', 'S', SLCTES, 1, A, 0, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('SGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sgges3('N', 'V', 'S', SLCTES, 1, A, 1, B, 0, SDIM, R1, R2, R3, Q, 1, U, 1, W, 1, BW, INFO );
         chkxer('SGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         sgges3('N', 'V', 'S', SLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 0, U, 1, W, 1, BW, INFO );
         chkxer('SGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         sgges3('V', 'V', 'S', SLCTES, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 1, U, 2, W, 1, BW, INFO );
         chkxer('SGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         sgges3('N', 'V', 'S', SLCTES, 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 0, W, 1, BW, INFO );
         chkxer('SGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         sgges3('V', 'V', 'S', SLCTES, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 1, W, 1, BW, INFO );
         chkxer('SGGES3 ', INFOT, NOUT, LERR, OK );
         INFOT = 19;
         sgges3('V', 'V', 'S', SLCTES, 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 2, W, 1, BW, INFO );
         chkxer('SGGES3 ', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // SGGESX

        srnamc.SRNAMT = 'SGGESX';
         INFOT = 1;
         sggesx('/', 'N', 'S', SLCTSX, 'N', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sggesx('N', '/', 'S', SLCTSX, 'N', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sggesx('V', 'V', '/', SLCTSX, 'N', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sggesx('V', 'V', 'S', SLCTSX, '/', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sggesx('V', 'V', 'S', SLCTSX, 'B', -1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sggesx('V', 'V', 'S', SLCTSX, 'B', 1, A, 0, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sggesx('V', 'V', 'S', SLCTSX, 'B', 1, A, 1, B, 0, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         sggesx('V', 'V', 'S', SLCTSX, 'B', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 0, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         sggesx('V', 'V', 'S', SLCTSX, 'B', 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         sggesx('V', 'V', 'S', SLCTSX, 'B', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 0, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         sggesx('V', 'V', 'S', SLCTSX, 'B', 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 1, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         sggesx('V', 'V', 'S', SLCTSX, 'B', 2, A, 2, B, 2, SDIM, R1, R2, R3, Q, 2, U, 2, RCE, RCV, W, 1, IW, 1, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 24;
         sggesx('V', 'V', 'S', SLCTSX, 'V', 1, A, 1, B, 1, SDIM, R1, R2, R3, Q, 1, U, 1, RCE, RCV, W, 32, IW, 0, BW, INFO );
         chkxer('SGGESX', INFOT, NOUT, LERR, OK );
         NT = NT + 13;

         // SGGEV

        srnamc.SRNAMT = 'SGGEV ';
         INFOT = 1;
         sggev('/', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('SGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sggev('N', '/', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('SGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sggev('V', 'V', -1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('SGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sggev('V', 'V', 1, A, 0, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('SGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sggev('V', 'V', 1, A, 1, B, 0, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('SGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sggev('N', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 0, U, 1, W, 1, INFO );
         chkxer('SGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sggev('V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 1, U, 2, W, 1, INFO );
         chkxer('SGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         sggev('V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 0, W, 1, INFO );
         chkxer('SGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         sggev('V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 1, W, 1, INFO );
         chkxer('SGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         sggev('V', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('SGGEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // SGGEV3

        srnamc.SRNAMT = 'SGGEV3 ';
         INFOT = 1;
         sggev3('/', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('SGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sggev3('N', '/', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('SGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sggev3('V', 'V', -1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('SGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sggev3('V', 'V', 1, A, 0, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('SGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sggev3('V', 'V', 1, A, 1, B, 0, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('SGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sggev3('N', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 0, U, 1, W, 1, INFO );
         chkxer('SGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sggev3('V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 1, U, 2, W, 1, INFO );
         chkxer('SGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         sggev3('V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 0, W, 1, INFO );
         chkxer('SGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         sggev3('V', 'V', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 1, W, 1, INFO );
         chkxer('SGGEV3 ', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         sggev3('V', 'V', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, W, 1, INFO );
         chkxer('SGGEV3 ', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // SGGEVX

        srnamc.SRNAMT = 'SGGEVX';
         INFOT = 1;
         sggevx('/', 'N', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('SGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sggevx('N', '/', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('SGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sggevx('N', 'N', '/', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('SGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sggevx('N', 'N', 'N', '/', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('SGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sggevx('N', 'N', 'N', 'N', -1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('SGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sggevx('N', 'N', 'N', 'N', 1, A, 0, B, 1, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('SGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sggevx('N', 'N', 'N', 'N', 1, A, 1, B, 0, R1, R2, R3, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('SGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         sggevx('N', 'N', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 0, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('SGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         sggevx('N', 'V', 'N', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 1, U, 2, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('SGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         sggevx('N', 'N', 'N', 'N', 1, A, 1, B, 1, R1, R2, R3, Q, 1, U, 0, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('SGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         sggevx('N', 'N', 'V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('SGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 26;
         sggevx('N', 'N', 'V', 'N', 2, A, 2, B, 2, R1, R2, R3, Q, 2, U, 2, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, IW, BW, INFO );
         chkxer('SGGEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // STGEXC

        srnamc.SRNAMT = 'STGEXC';
         INFOT = 3;
         stgexc( true , true , -1, A, 1, B, 1, Q, 1, Z, 1, IFST, ILST, W, 1, INFO );
         chkxer('STGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         stgexc( true , true , 1, A, 0, B, 1, Q, 1, Z, 1, IFST, ILST, W, 1, INFO );
         chkxer('STGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         stgexc( true , true , 1, A, 1, B, 0, Q, 1, Z, 1, IFST, ILST, W, 1, INFO );
         chkxer('STGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         stgexc( false , true , 1, A, 1, B, 1, Q, 0, Z, 1, IFST, ILST, W, 1, INFO );
         chkxer('STGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         stgexc( true , true , 1, A, 1, B, 1, Q, 0, Z, 1, IFST, ILST, W, 1, INFO );
         chkxer('STGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         stgexc( true , false , 1, A, 1, B, 1, Q, 1, Z, 0, IFST, ILST, W, 1, INFO );
         chkxer('STGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         stgexc( true , true , 1, A, 1, B, 1, Q, 1, Z, 0, IFST, ILST, W, 1, INFO );
         chkxer('STGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         stgexc( true , true , 1, A, 1, B, 1, Q, 1, Z, 1, IFST, ILST, W, 0, INFO );
         chkxer('STGEXC', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

         // STGSEN

        srnamc.SRNAMT = 'STGSEN';
         INFOT = 1;
         stgsen(-1, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('STGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         stgsen(1, true , true , SEL, -1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('STGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         stgsen(1, true , true , SEL, 1, A, 0, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('STGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         stgsen(1, true , true , SEL, 1, A, 1, B, 0, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('STGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         stgsen(1, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 0, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('STGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         stgsen(1, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 0, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('STGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         stgsen(0, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('STGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         stgsen(1, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('STGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         stgsen(2, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('STGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 24;
         stgsen(0, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 0, INFO );
         chkxer('STGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 24;
         stgsen(1, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 0, INFO );
         chkxer('STGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 24;
         stgsen(2, true , true , SEL, 1, A, 1, B, 1, R1, R2, R3, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 1, INFO );
         chkxer('STGSEN', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // STGSNA

        srnamc.SRNAMT = 'STGSNA';
         INFOT = 1;
         stgsna('/', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('STGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         stgsna('B', '/', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('STGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         stgsna('B', 'A', SEL, -1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('STGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         stgsna('B', 'A', SEL, 1, A, 0, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('STGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         stgsna('B', 'A', SEL, 1, A, 1, B, 0, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('STGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         stgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 0, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('STGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         stgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 0, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('STGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         stgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 0, M, W, 1, IW, INFO );
         chkxer('STGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         stgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 0, IW, INFO );
         chkxer('STGSNA', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // STGSYL

        srnamc.SRNAMT = 'STGSYL';
         INFOT = 1;
         stgsyl('/', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('STGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         stgsyl('N', -1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('STGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         stgsyl('N', 0, 0, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('STGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         stgsyl('N', 0, 1, 0, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('STGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         stgsyl('N', 0, 1, 1, A, 0, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('STGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         stgsyl('N', 0, 1, 1, A, 1, B, 0, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('STGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         stgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 0, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('STGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         stgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 0, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('STGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         stgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 0, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('STGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         stgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 0, SCALE, DIF, W, 1, IW, INFO );
         chkxer('STGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         stgsyl('N', 1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('STGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         stgsyl('N', 2, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('STGSYL', INFOT, NOUT, LERR, OK );
         NT = NT + 12;
      }

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH, NT;
      } else {
         WRITE( NOUT, FMT = 9998 )PATH;
      }

 9999 FORMAT(' ${.a3} routines passed the tests of the error exits (${.i3} tests done)' );
 9998 FORMAT( ' *** ${.a3} routines failed the tests of the error exits ***' );

      return;
      }
