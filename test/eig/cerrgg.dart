      void cerrgg(PATH, NUNIT ) {

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
      const              NMAX = 3, LW = 6*NMAX ;
      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                DUMMYK, DUMMYL, I, IFST, IHI, ILO, ILST, INFO, J, M, NCYCLE, NT, SDIM, LWORK;
      double               ANRM, BNRM, DIF, SCALE, TOLA, TOLB;
      // ..
      // .. Local Arrays ..
      bool               BW( NMAX ), SEL( NMAX );
      int                IW( LW ), IDUM(NMAX);
      double               LS( NMAX ), R1( NMAX ), R2( NMAX ), RCE( NMAX ), RCV( NMAX ), RS( NMAX ), RW( LW )       Complex            A( NMAX, NMAX ), ALPHA( NMAX ), B( NMAX, NMAX ), BETA( NMAX ), Q( NMAX, NMAX ), TAU( NMAX ), U( NMAX, NMAX ), V( NMAX, NMAX ), W( LW ), Z( NMAX, NMAX );
      // ..
      // .. External Functions ..
      //- bool               CLCTES, CLCTSX, LSAMEN;
      // EXTERNAL CLCTES, CLCTSX, LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGGES, CGGESX, CGGEV, CGGEVX, CGGGLM, CGGHRD, CGGLSE, CGGQRF, CGGRQF, CHGEQZ, CHKXER, CTGEVC, CTGEXC, CTGSEN, CTGSJA, CTGSNA, CTGSYL, CUNCSD, CGGES3, CGGEV3, CGGHD3, CGGSVD3, CGGSVP3, XLAENV
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
         SEL[J] = true;
         for (I = 1; I <= NMAX; I++) { // 10
            A[I, J] = ZERO;
            B[I, J] = ZERO;
         } // 10
      } // 20
      for (I = 1; I <= NMAX; I++) { // 30
         A[I, I] = ONE;
         B[I, I] = ONE;
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

         // CGGHRD

         SRNAMT = 'CGGHRD';
         INFOT = 1;
         cgghrd('/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('CGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgghrd('N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('CGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgghrd('N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('CGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgghrd('N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('CGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgghrd('N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('CGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgghrd('N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, INFO );
         chkxer('CGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cgghrd('N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, INFO );
         chkxer('CGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cgghrd('V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO );
         chkxer('CGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         cgghrd('N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO );
         chkxer('CGGHRD', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // CGGHD3

         SRNAMT = 'CGGHD3';
         INFOT = 1;
         cgghd3('/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('CGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgghd3('N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('CGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgghd3('N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('CGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cgghd3('N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('CGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgghd3('N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('CGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgghd3('N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, W, LW, INFO );
         chkxer('CGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cgghd3('N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('CGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cgghd3('V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW, INFO );
         chkxer('CGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         cgghd3('N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW, INFO );
         chkxer('CGGHD3', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // CHGEQZ

         SRNAMT = 'CHGEQZ';
         INFOT = 1;
         chgeqz('/', 'N', 'N', 0, 1, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('CHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chgeqz('E', '/', 'N', 0, 1, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('CHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chgeqz('E', 'N', '/', 0, 1, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('CHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chgeqz('E', 'N', 'N', -1, 0, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('CHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chgeqz('E', 'N', 'N', 0, 0, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('CHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         chgeqz('E', 'N', 'N', 0, 1, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('CHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chgeqz('E', 'N', 'N', 2, 1, 1, A, 1, B, 2, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('CHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         chgeqz('E', 'N', 'N', 2, 1, 1, A, 2, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('CHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         chgeqz('E', 'V', 'N', 2, 1, 1, A, 2, B, 2, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('CHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         chgeqz('E', 'N', 'V', 2, 1, 1, A, 2, B, 2, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('CHGEQZ', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // CTGEVC

         SRNAMT = 'CTGEVC';
         INFOT = 1;
         ctgevc('/', 'A', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, RW, INFO );
         chkxer('CTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctgevc('R', '/', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, RW, INFO );
         chkxer('CTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctgevc('R', 'A', SEL, -1, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, RW, INFO );
         chkxer('CTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ctgevc('R', 'A', SEL, 2, A, 1, B, 2, Q, 1, Z, 2, 0, M, W, RW, INFO );
         chkxer('CTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ctgevc('R', 'A', SEL, 2, A, 2, B, 1, Q, 1, Z, 2, 0, M, W, RW, INFO );
         chkxer('CTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ctgevc('L', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W, RW, INFO );
         chkxer('CTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         ctgevc('R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W, RW, INFO );
         chkxer('CTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         ctgevc('R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 2, 1, M, W, RW, INFO );
         chkxer('CTGEVC', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the GSV path.

      } else if ( lsamen( 3, PATH, 'GSV' ) ) {

         // CGGSVD3

         SRNAMT = 'CGGSVD3';
         INFOT = 1;
         cggsvd3('/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('CGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cggsvd3('N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('CGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cggsvd3('N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('CGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cggsvd3('N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('CGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cggsvd3('N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('CGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cggsvd3('N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('CGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cggsvd3('N', 'N', 'N', 2, 1, 1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('CGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cggsvd3('N', 'N', 'N', 1, 1, 2, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('CGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         cggsvd3('U', 'N', 'N', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B, 2, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('CGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         cggsvd3('N', 'V', 'N', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B, 2, R1, R2, U, 2, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('CGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         cggsvd3('N', 'N', 'Q', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B, 2, R1, R2, U, 2, V, 2, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('CGGSVD3', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // CGGSVP3

         SRNAMT = 'CGGSVP3';
         INFOT = 1;
         cggsvp3('/', 'N', 'N', 0, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('CGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cggsvp3('N', '/', 'N', 0, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('CGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cggsvp3('N', 'N', '/', 0, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('CGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cggsvp3('N', 'N', 'N', -1, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('CGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cggsvp3('N', 'N', 'N', 0, -1, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('CGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cggsvp3('N', 'N', 'N', 0, 0, -1, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('CGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cggsvp3('N', 'N', 'N', 2, 1, 1, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('CGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cggsvp3('N', 'N', 'N', 1, 2, 1, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('CGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         cggsvp3('U', 'N', 'N', 2, 2, 2, A, 2, B, 2, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('CGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         cggsvp3('N', 'V', 'N', 2, 2, 2, A, 2, B, 2, TOLA, TOLB, DUMMYK, DUMMYL, U, 2, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('CGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         cggsvp3('N', 'N', 'Q', 2, 2, 2, A, 2, B, 2, TOLA, TOLB, DUMMYK, DUMMYL, U, 2, V, 2, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('CGGSVP3', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // CTGSJA

         SRNAMT = 'CTGSJA';
         INFOT = 1;
         ctgsja('/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('CTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctgsja('N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('CTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctgsja('N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('CTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctgsja('N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('CTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ctgsja('N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('CTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ctgsja('N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('CTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ctgsja('N', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 0, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('CTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         ctgsja('N', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 0, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('CTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         ctgsja('U', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 0, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('CTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         ctgsja('N', 'V', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 0, Q, 1, W, NCYCLE, INFO );
         chkxer('CTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         ctgsja('N', 'N', 'Q', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 0, W, NCYCLE, INFO );
         chkxer('CTGSJA', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

      // Test error exits for the GLM path.

      } else if ( lsamen( 3, PATH, 'GLM' ) ) {

         // CGGGLM

         SRNAMT = 'CGGGLM';
         INFOT = 1;
         cggglm(-1, 0, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cggglm(0, -1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cggglm(0, 1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cggglm(0, 0, -1, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cggglm(1, 0, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cggglm(0, 0, 0, A, 0, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cggglm(0, 0, 0, A, 1, B, 0, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cggglm(1, 1, 1, A, 1, B, 1, TAU, ALPHA, BETA, W, 1, INFO );
         chkxer('CGGGLM', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the LSE path.

      } else if ( lsamen( 3, PATH, 'LSE' ) ) {

         // CGGLSE

         SRNAMT = 'CGGLSE';
         INFOT = 1;
         cgglse(-1, 0, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgglse(0, -1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgglse(0, 0, -1, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgglse(0, 0, 1, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgglse(0, 1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgglse(0, 0, 0, A, 0, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgglse(0, 0, 0, A, 1, B, 0, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('CGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cgglse(1, 1, 1, A, 1, B, 1, TAU, ALPHA, BETA, W, 1, INFO );
         chkxer('CGGLSE', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the CSD path.

      } else if ( lsamen( 3, PATH, 'CSD' ) ) {

         // CUNCSD

         SRNAMT = 'CUNCSD';
         INFOT = 7;
         cuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', -1, 0, 0, A, 1, A, 1, A, 1, A, 1, RS, A, 1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('CUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, -1, 0, A, 1, A, 1, A, 1, A, 1, RS, A, 1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('CUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, -1, A, 1, A, 1, A, 1, A, 1, RS, A, 1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('CUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, -1, A, 1, A, 1, A, 1, RS, A, 1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('CUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         cuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, RS, A, -1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('CUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         cuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, RS, A, 1, A, -1, A, 1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('CUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 24;
         cuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, RS, A, 1, A, 1, A, -1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('CUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 26;
         cuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, RS, A, 1, A, 1, A, 1, A, -1, W, LW, RW, LW, IW, INFO );
         chkxer('CUNCSD', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the GQR path.

      } else if ( lsamen( 3, PATH, 'GQR' ) ) {

         // CGGQRF

         SRNAMT = 'CGGQRF';
         INFOT = 1;
         cggqrf(-1, 0, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('CGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cggqrf(0, -1, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('CGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cggqrf(0, 0, -1, A, 1, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('CGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cggqrf(0, 0, 0, A, 0, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('CGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cggqrf(0, 0, 0, A, 1, ALPHA, B, 0, BETA, W, LW, INFO );
         chkxer('CGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cggqrf(1, 1, 2, A, 1, ALPHA, B, 1, BETA, W, 1, INFO );
         chkxer('CGGQRF', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // CGGRQF

         SRNAMT = 'CGGRQF';
         INFOT = 1;
         cggrqf(-1, 0, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('CGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cggrqf(0, -1, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('CGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cggrqf(0, 0, -1, A, 1, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('CGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cggrqf(0, 0, 0, A, 0, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('CGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cggrqf(0, 0, 0, A, 1, ALPHA, B, 0, BETA, W, LW, INFO );
         chkxer('CGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cggrqf(1, 1, 2, A, 1, ALPHA, B, 1, BETA, W, 1, INFO );
         chkxer('CGGRQF', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

      // Test error exits for the CGS, CGV, CGX, and CXV paths.

      } else if ( lsamen( 3, PATH, 'CGS' ) || lsamen( 3, PATH, 'CGV' ) || lsamen( 3, PATH, 'CGX' ) || lsamen( 3, PATH, 'CXV' ) ) {

         // CGGES

         SRNAMT = 'CGGES ';
         INFOT = 1;
         cgges('/', 'N', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgges('N', '/', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgges('N', 'V', '/', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgges('N', 'V', 'S', CLCTES, -1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgges('N', 'V', 'S', CLCTES, 1, A, 0, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cgges('N', 'V', 'S', CLCTES, 1, A, 1, B, 0, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         cgges('N', 'V', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 0, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         cgges('V', 'V', 'S', CLCTES, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 1, U, 2, W, 1, RW, BW, INFO );
         chkxer('CGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         cgges('N', 'V', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 0, W, 1, RW, BW, INFO );
         chkxer('CGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         cgges('V', 'V', 'S', CLCTES, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         cgges('V', 'V', 'S', CLCTES, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 2, W, 1, RW, BW, INFO );
         chkxer('CGGES ', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // CGGES3

         SRNAMT = 'CGGES3';
         INFOT = 1;
         cgges3('/', 'N', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cgges3('N', '/', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cgges3('N', 'V', '/', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cgges3('N', 'V', 'S', CLCTES, -1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cgges3('N', 'V', 'S', CLCTES, 1, A, 0, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cgges3('N', 'V', 'S', CLCTES, 1, A, 1, B, 0, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         cgges3('N', 'V', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 0, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         cgges3('V', 'V', 'S', CLCTES, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 1, U, 2, W, 1, RW, BW, INFO );
         chkxer('CGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         cgges3('N', 'V', 'S', CLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 0, W, 1, RW, BW, INFO );
         chkxer('CGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         cgges3('V', 'V', 'S', CLCTES, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 1, W, 1, RW, BW, INFO );
         chkxer('CGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         cgges3('V', 'V', 'S', CLCTES, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 2, W, 1, RW, BW, INFO );
         chkxer('CGGES3', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // CGGESX

         SRNAMT = 'CGGESX';
         INFOT = 1;
         cggesx('/', 'N', 'S', CLCTSX, 'N', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cggesx('N', '/', 'S', CLCTSX, 'N', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cggesx('V', 'V', '/', CLCTSX, 'N', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cggesx('V', 'V', 'S', CLCTSX, '/', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cggesx('V', 'V', 'S', CLCTSX, 'B', -1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cggesx('V', 'V', 'S', CLCTSX, 'B', 1, A, 0, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cggesx('V', 'V', 'S', CLCTSX, 'B', 1, A, 1, B, 0, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cggesx('V', 'V', 'S', CLCTSX, 'B', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 0, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cggesx('V', 'V', 'S', CLCTSX, 'B', 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         cggesx('V', 'V', 'S', CLCTSX, 'B', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 0, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         cggesx('V', 'V', 'S', CLCTSX, 'B', 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 21;
         cggesx('V', 'V', 'S', CLCTSX, 'B', 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 2, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 24;
         cggesx('V', 'V', 'S', CLCTSX, 'V', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 32, RW, IW, 0, BW, INFO );
         chkxer('CGGESX', INFOT, NOUT, LERR, OK );
         NT = NT + 13;

         // CGGEV

         SRNAMT = 'CGGEV ';
         INFOT = 1;
         cggev('/', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cggev('N', '/', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cggev('V', 'V', -1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cggev('V', 'V', 1, A, 0, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cggev('V', 'V', 1, A, 1, B, 0, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cggev('N', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 0, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cggev('V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 1, U, 2, W, 1, RW, INFO );
         chkxer('CGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         cggev('V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 0, W, 1, RW, INFO );
         chkxer('CGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         cggev('V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cggev('V', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // CGGEV3

         SRNAMT = 'CGGEV3';
         INFOT = 1;
         cggev3('/', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cggev3('N', '/', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cggev3('V', 'V', -1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cggev3('V', 'V', 1, A, 0, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cggev3('V', 'V', 1, A, 1, B, 0, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cggev3('N', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 0, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         cggev3('V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 1, U, 2, W, 1, RW, INFO );
         chkxer('CGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         cggev3('V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 0, W, 1, RW, INFO );
         chkxer('CGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         cggev3('V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cggev3('V', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('CGGEV3', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // CGGEVX

         SRNAMT = 'CGGEVX';
         INFOT = 1;
         cggevx('/', 'N', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('CGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cggevx('N', '/', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('CGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cggevx('N', 'N', '/', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('CGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cggevx('N', 'N', 'N', '/', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('CGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cggevx('N', 'N', 'N', 'N', -1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('CGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cggevx('N', 'N', 'N', 'N', 1, A, 0, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('CGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         cggevx('N', 'N', 'N', 'N', 1, A, 1, B, 0, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('CGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         cggevx('N', 'N', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 0, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('CGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         cggevx('N', 'V', 'N', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 1, U, 2, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('CGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cggevx('N', 'N', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 0, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('CGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         cggevx('N', 'N', 'V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('CGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 25;
         cggevx('N', 'N', 'V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 2, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 0, RW, IW, BW, INFO );
         chkxer('CGGEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // CTGEXC

         SRNAMT = 'CTGEXC';
         INFOT = 3;
         ctgexc( true , true , -1, A, 1, B, 1, Q, 1, Z, 1, IFST, ILST, INFO );
         chkxer('CTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ctgexc( true , true , 1, A, 0, B, 1, Q, 1, Z, 1, IFST, ILST, INFO );
         chkxer('CTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ctgexc( true , true , 1, A, 1, B, 0, Q, 1, Z, 1, IFST, ILST, INFO );
         chkxer('CTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ctgexc( false , true , 1, A, 1, B, 1, Q, 0, Z, 1, IFST, ILST, INFO );
         chkxer('CTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ctgexc( true , true , 1, A, 1, B, 1, Q, 0, Z, 1, IFST, ILST, INFO );
         chkxer('CTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ctgexc( true , false , 1, A, 1, B, 1, Q, 1, Z, 0, IFST, ILST, INFO );
         chkxer('CTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ctgexc( true , true , 1, A, 1, B, 1, Q, 1, Z, 0, IFST, ILST, INFO );
         chkxer('CTGEXC', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // CTGSEN

         SRNAMT = 'CTGSEN';
         INFOT = 1;
         ctgsen(-1, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('CTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ctgsen(1, true , true , SEL, -1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('CTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ctgsen(1, true , true , SEL, 1, A, 0, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('CTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ctgsen(1, true , true , SEL, 1, A, 1, B, 0, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('CTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         ctgsen(1, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 0, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('CTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         ctgsen(1, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 0, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('CTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 21;
         ctgsen(3, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, -5, IW, 1, INFO );
         chkxer('CTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 23;
         ctgsen(0, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 0, INFO );
         chkxer('CTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 23;
         ctgsen(1, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 0, INFO );
         chkxer('CTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 23;
         ctgsen(5, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 1, INFO );
         chkxer('CTGSEN', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // CTGSNA

         SRNAMT = 'CTGSNA';
         INFOT = 1;
         ctgsna('/', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('CTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctgsna('B', '/', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('CTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctgsna('B', 'A', SEL, -1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('CTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ctgsna('B', 'A', SEL, 1, A, 0, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('CTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ctgsna('B', 'A', SEL, 1, A, 1, B, 0, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('CTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ctgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 0, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('CTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         ctgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 0, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('CTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         ctgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 0, M, W, 1, IW, INFO );
         chkxer('CTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         ctgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 0, IW, INFO );
         chkxer('CTGSNA', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // CTGSYL

         SRNAMT = 'CTGSYL';
         INFOT = 1;
         ctgsyl('/', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('CTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ctgsyl('N', -1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('CTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ctgsyl('N', 0, 0, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('CTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ctgsyl('N', 0, 1, 0, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('CTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ctgsyl('N', 0, 1, 1, A, 0, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('CTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ctgsyl('N', 0, 1, 1, A, 1, B, 0, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('CTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ctgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 0, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('CTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         ctgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 0, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('CTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         ctgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 0, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('CTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 16;
         ctgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 0, SCALE, DIF, W, 1, IW, INFO );
         chkxer('CTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         ctgsyl('N', 1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('CTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         ctgsyl('N', 2, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('CTGSYL', INFOT, NOUT, LERR, OK );
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

      return;
      }
