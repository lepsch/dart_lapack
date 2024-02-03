      SUBROUTINE ZERRGG( PATH, NUNIT )

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
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                DUMMYK, DUMMYL, I, IFST, IHI, ILO, ILST, INFO, J, M, NCYCLE, NT, SDIM, LWORK;
      double             ANRM, BNRM, DIF, SCALE, TOLA, TOLB;
      // ..
      // .. Local Arrays ..
      bool               BW( NMAX ), SEL( NMAX );
      int                IW( LW ), IDUM(NMAX);
      double             LS( NMAX ), R1( NMAX ), R2( NMAX ), RCE( NMAX ), RCV( NMAX ), RS( NMAX ), RW( LW );
      COMPLEX*16         A( NMAX, NMAX ), ALPHA( NMAX ), B( NMAX, NMAX ), BETA( NMAX ), Q( NMAX, NMAX ), TAU( NMAX ), U( NMAX, NMAX ), V( NMAX, NMAX ), W( LW ), Z( NMAX, NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN, ZLCTES, ZLCTSX;
      // EXTERNAL LSAMEN, ZLCTES, ZLCTSX
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZGGES,  ZGGESX, ZGGEV,  ZGGEVX, ZGGGLM, ZGGHRD, ZGGLSE, ZGGQRF, ZGGRQF, ZHGEQZ, ZTGEVC, ZTGEXC, ZTGSEN, ZTGSJA, ZTGSNA, ZTGSYL, ZUNCSD, ZGGES3, ZGGEV3, ZGGHD3, ZGGSVD3, ZGGSVP3, XLAENV
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

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         SEL( J ) = true;
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = ZERO
            B( I, J ) = ZERO
         } // 10
      } // 20
      for (I = 1; I <= NMAX; I++) { // 30
         A( I, I ) = ONE
         B( I, I ) = ONE
      } // 30
      OK = true;
      TOLA = 1.0D0
      TOLB = 1.0D0
      IFST = 1
      ILST = 1
      NT = 0
      LWORK = 1

      // Call XLAENV to set the parameters used in CLAQZ0

      xlaenv(12, 10 );
      xlaenv(13, 12 );
      xlaenv(14, 13 );
      xlaenv(15, 2 );
      xlaenv(17, 10 );

      // Test error exits for the GG path.

      if ( LSAMEN( 2, C2, 'GG' ) ) {

         // ZGGHRD

         SRNAMT = 'ZGGHRD'
         INFOT = 1
         zgghrd('/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('ZGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgghrd('N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('ZGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgghrd('N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('ZGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zgghrd('N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('ZGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zgghrd('N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, INFO );
         chkxer('ZGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zgghrd('N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, INFO );
         chkxer('ZGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zgghrd('N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, INFO );
         chkxer('ZGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zgghrd('V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO );
         chkxer('ZGGHRD', INFOT, NOUT, LERR, OK );
         INFOT = 13
         zgghrd('N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO );
         chkxer('ZGGHRD', INFOT, NOUT, LERR, OK );
         NT = NT + 9

         // ZGGHD3

         SRNAMT = 'ZGGHD3'
         INFOT = 1
         zgghd3('/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('ZGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgghd3('N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('ZGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgghd3('N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('ZGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zgghd3('N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('ZGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zgghd3('N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('ZGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zgghd3('N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, W, LW, INFO );
         chkxer('ZGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zgghd3('N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, W, LW, INFO );
         chkxer('ZGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zgghd3('V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW, INFO );
         chkxer('ZGGHD3', INFOT, NOUT, LERR, OK );
         INFOT = 13
         zgghd3('N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW, INFO );
         chkxer('ZGGHD3', INFOT, NOUT, LERR, OK );
         NT = NT + 9

         // ZHGEQZ

         SRNAMT = 'ZHGEQZ'
         INFOT = 1
         zhgeqz('/', 'N', 'N', 0, 1, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('ZHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhgeqz('E', '/', 'N', 0, 1, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('ZHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhgeqz('E', 'N', '/', 0, 1, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('ZHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhgeqz('E', 'N', 'N', -1, 0, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('ZHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zhgeqz('E', 'N', 'N', 0, 0, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('ZHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zhgeqz('E', 'N', 'N', 0, 1, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('ZHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhgeqz('E', 'N', 'N', 2, 1, 1, A, 1, B, 2, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('ZHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zhgeqz('E', 'N', 'N', 2, 1, 1, A, 2, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('ZHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 14
         zhgeqz('E', 'V', 'N', 2, 1, 1, A, 2, B, 2, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('ZHGEQZ', INFOT, NOUT, LERR, OK );
         INFOT = 16
         zhgeqz('E', 'N', 'V', 2, 1, 1, A, 2, B, 2, ALPHA, BETA, Q, 1, Z, 1, W, 1, RW, INFO );
         chkxer('ZHGEQZ', INFOT, NOUT, LERR, OK );
         NT = NT + 10

         // ZTGEVC

         SRNAMT = 'ZTGEVC'
         INFOT = 1
         ztgevc('/', 'A', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, RW, INFO );
         chkxer('ZTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ztgevc('R', '/', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, RW, INFO );
         chkxer('ZTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 4
         ztgevc('R', 'A', SEL, -1, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, RW, INFO );
         chkxer('ZTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 6
         ztgevc('R', 'A', SEL, 2, A, 1, B, 2, Q, 1, Z, 2, 0, M, W, RW, INFO );
         chkxer('ZTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 8
         ztgevc('R', 'A', SEL, 2, A, 2, B, 1, Q, 1, Z, 2, 0, M, W, RW, INFO );
         chkxer('ZTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 10
         ztgevc('L', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W, RW, INFO );
         chkxer('ZTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 12
         ztgevc('R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W, RW, INFO );
         chkxer('ZTGEVC', INFOT, NOUT, LERR, OK );
         INFOT = 13
         ztgevc('R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 2, 1, M, W, RW, INFO );
         chkxer('ZTGEVC', INFOT, NOUT, LERR, OK );
         NT = NT + 8

      // Test error exits for the GSV path.

      } else if ( LSAMEN( 3, PATH, 'GSV' ) ) {

         // ZGGSVD3

         SRNAMT = 'ZGGSVD3'
         INFOT = 1
         zggsvd3('/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('ZGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zggsvd3('N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('ZGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zggsvd3('N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('ZGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zggsvd3('N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('ZGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zggsvd3('N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('ZGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zggsvd3('N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('ZGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zggsvd3('N', 'N', 'N', 2, 1, 1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('ZGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 12
         zggsvd3('N', 'N', 'N', 1, 1, 2, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('ZGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 16
         zggsvd3('U', 'N', 'N', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B, 2, R1, R2, U, 1, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('ZGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 18
         zggsvd3('N', 'V', 'N', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B, 2, R1, R2, U, 2, V, 1, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('ZGGSVD3', INFOT, NOUT, LERR, OK );
         INFOT = 20
         zggsvd3('N', 'N', 'Q', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B, 2, R1, R2, U, 2, V, 2, Q, 1, W, LWORK, RW, IDUM, INFO );
         chkxer('ZGGSVD3', INFOT, NOUT, LERR, OK );
         NT = NT + 11

         // ZGGSVP3

         SRNAMT = 'ZGGSVP3'
         INFOT = 1
         zggsvp3('/', 'N', 'N', 0, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('ZGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zggsvp3('N', '/', 'N', 0, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('ZGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zggsvp3('N', 'N', '/', 0, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('ZGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zggsvp3('N', 'N', 'N', -1, 0, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('ZGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zggsvp3('N', 'N', 'N', 0, -1, 0, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('ZGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zggsvp3('N', 'N', 'N', 0, 0, -1, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('ZGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zggsvp3('N', 'N', 'N', 2, 1, 1, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('ZGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zggsvp3('N', 'N', 'N', 1, 2, 1, A, 1, B, 1, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('ZGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 16
         zggsvp3('U', 'N', 'N', 2, 2, 2, A, 2, B, 2, TOLA, TOLB, DUMMYK, DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('ZGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 18
         zggsvp3('N', 'V', 'N', 2, 2, 2, A, 2, B, 2, TOLA, TOLB, DUMMYK, DUMMYL, U, 2, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('ZGGSVP3', INFOT, NOUT, LERR, OK );
         INFOT = 20
         zggsvp3('N', 'N', 'Q', 2, 2, 2, A, 2, B, 2, TOLA, TOLB, DUMMYK, DUMMYL, U, 2, V, 2, Q, 1, IW, RW, TAU, W, LWORK, INFO );
         chkxer('ZGGSVP3', INFOT, NOUT, LERR, OK );
         NT = NT + 11

         // ZTGSJA

         SRNAMT = 'ZTGSJA'
         INFOT = 1
         ztgsja('/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('ZTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ztgsja('N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('ZTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 3
         ztgsja('N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('ZTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 4
         ztgsja('N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('ZTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 5
         ztgsja('N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('ZTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 6
         ztgsja('N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('ZTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 10
         ztgsja('N', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 0, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('ZTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 12
         ztgsja('N', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 0, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('ZTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 18
         ztgsja('U', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 0, V, 1, Q, 1, W, NCYCLE, INFO );
         chkxer('ZTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 20
         ztgsja('N', 'V', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 0, Q, 1, W, NCYCLE, INFO );
         chkxer('ZTGSJA', INFOT, NOUT, LERR, OK );
         INFOT = 22
         ztgsja('N', 'N', 'Q', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, TOLA, TOLB, R1, R2, U, 1, V, 1, Q, 0, W, NCYCLE, INFO );
         chkxer('ZTGSJA', INFOT, NOUT, LERR, OK );
         NT = NT + 11

      // Test error exits for the GLM path.

      } else if ( LSAMEN( 3, PATH, 'GLM' ) ) {

         // ZGGGLM

         SRNAMT = 'ZGGGLM'
         INFOT = 1
         zggglm(-1, 0, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zggglm(0, -1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zggglm(0, 1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zggglm(0, 0, -1, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zggglm(1, 0, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zggglm(0, 0, 0, A, 0, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zggglm(0, 0, 0, A, 1, B, 0, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGGLM', INFOT, NOUT, LERR, OK );
         INFOT = 12
         zggglm(1, 1, 1, A, 1, B, 1, TAU, ALPHA, BETA, W, 1, INFO );
         chkxer('ZGGGLM', INFOT, NOUT, LERR, OK );
         NT = NT + 8

      // Test error exits for the LSE path.

      } else if ( LSAMEN( 3, PATH, 'LSE' ) ) {

         // ZGGLSE

         SRNAMT = 'ZGGLSE'
         INFOT = 1
         zgglse(-1, 0, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgglse(0, -1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgglse(0, 0, -1, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgglse(0, 0, 1, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgglse(0, 1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zgglse(0, 0, 0, A, 0, B, 1, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zgglse(0, 0, 0, A, 1, B, 0, TAU, ALPHA, BETA, W, LW, INFO );
         chkxer('ZGGLSE', INFOT, NOUT, LERR, OK );
         INFOT = 12
         zgglse(1, 1, 1, A, 1, B, 1, TAU, ALPHA, BETA, W, 1, INFO );
         chkxer('ZGGLSE', INFOT, NOUT, LERR, OK );
         NT = NT + 8

      // Test error exits for the CSD path.

      } else if ( LSAMEN( 3, PATH, 'CSD' ) ) {

         // ZUNCSD

         SRNAMT = 'ZUNCSD'
         INFOT = 7
         zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', -1, 0, 0, A, 1, A, 1, A, 1, A, 1, RS, A, 1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('ZUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, -1, 0, A, 1, A, 1, A, 1, A, 1, RS, A, 1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('ZUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, -1, A, 1, A, 1, A, 1, A, 1, RS, A, 1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('ZUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, -1, A, 1, A, 1, A, 1, RS, A, 1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('ZUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 20
         zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, RS, A, -1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('ZUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 22
         zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, RS, A, 1, A, -1, A, 1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('ZUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 24
         zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, RS, A, 1, A, 1, A, -1, A, 1, W, LW, RW, LW, IW, INFO );
         chkxer('ZUNCSD', INFOT, NOUT, LERR, OK );
         INFOT = 26
         zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, RS, A, 1, A, 1, A, 1, A, -1, W, LW, RW, LW, IW, INFO );
         chkxer('ZUNCSD', INFOT, NOUT, LERR, OK );
         NT = NT + 8

      // Test error exits for the GQR path.

      } else if ( LSAMEN( 3, PATH, 'GQR' ) ) {

         // ZGGQRF

         SRNAMT = 'ZGGQRF'
         INFOT = 1
         zggqrf(-1, 0, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('ZGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zggqrf(0, -1, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('ZGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zggqrf(0, 0, -1, A, 1, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('ZGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zggqrf(0, 0, 0, A, 0, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('ZGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zggqrf(0, 0, 0, A, 1, ALPHA, B, 0, BETA, W, LW, INFO );
         chkxer('ZGGQRF', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zggqrf(1, 1, 2, A, 1, ALPHA, B, 1, BETA, W, 1, INFO );
         chkxer('ZGGQRF', INFOT, NOUT, LERR, OK );
         NT = NT + 6

         // ZGGRQF

         SRNAMT = 'ZGGRQF'
         INFOT = 1
         zggrqf(-1, 0, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('ZGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zggrqf(0, -1, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('ZGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zggrqf(0, 0, -1, A, 1, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('ZGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zggrqf(0, 0, 0, A, 0, ALPHA, B, 1, BETA, W, LW, INFO );
         chkxer('ZGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zggrqf(0, 0, 0, A, 1, ALPHA, B, 0, BETA, W, LW, INFO );
         chkxer('ZGGRQF', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zggrqf(1, 1, 2, A, 1, ALPHA, B, 1, BETA, W, 1, INFO );
         chkxer('ZGGRQF', INFOT, NOUT, LERR, OK );
         NT = NT + 6

      // Test error exits for the ZGS, ZGV, ZGX, and ZXV paths.

      } else if ( LSAMEN( 3, PATH, 'ZGS' ) .OR. LSAMEN( 3, PATH, 'ZGV' ) .OR. LSAMEN( 3, PATH, 'ZGX' ) .OR. LSAMEN( 3, PATH, 'ZXV' ) ) {

         // ZGGES

         SRNAMT = 'ZGGES '
         INFOT = 1
         zgges('/', 'N', 'S', ZLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgges('N', '/', 'S', ZLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgges('N', 'V', '/', ZLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zgges('N', 'V', 'S', ZLCTES, -1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zgges('N', 'V', 'S', ZLCTES, 1, A, 0, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zgges('N', 'V', 'S', ZLCTES, 1, A, 1, B, 0, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 14
         zgges('N', 'V', 'S', ZLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 0, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 14
         zgges('V', 'V', 'S', ZLCTES, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 1, U, 2, W, 1, RW, BW, INFO );
         chkxer('ZGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 16
         zgges('N', 'V', 'S', ZLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 0, W, 1, RW, BW, INFO );
         chkxer('ZGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 16
         zgges('V', 'V', 'S', ZLCTES, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES ', INFOT, NOUT, LERR, OK );
         INFOT = 18
         zgges('V', 'V', 'S', ZLCTES, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 2, W, 1, RW, BW, INFO );
         chkxer('ZGGES ', INFOT, NOUT, LERR, OK );
         NT = NT + 11

         // ZGGES3

         SRNAMT = 'ZGGES3'
         INFOT = 1
         zgges3('/', 'N', 'S', ZLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zgges3('N', '/', 'S', ZLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zgges3('N', 'V', '/', ZLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zgges3('N', 'V', 'S', ZLCTES, -1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zgges3('N', 'V', 'S', ZLCTES, 1, A, 0, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zgges3('N', 'V', 'S', ZLCTES, 1, A, 1, B, 0, SDIM, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 14
         zgges3('N', 'V', 'S', ZLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 0, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 14
         zgges3('V', 'V', 'S', ZLCTES, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 1, U, 2, W, 1, RW, BW, INFO );
         chkxer('ZGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 16
         zgges3('N', 'V', 'S', ZLCTES, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 0, W, 1, RW, BW, INFO );
         chkxer('ZGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 16
         zgges3('V', 'V', 'S', ZLCTES, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 1, W, 1, RW, BW, INFO );
         chkxer('ZGGES3', INFOT, NOUT, LERR, OK );
         INFOT = 18
         zgges3('V', 'V', 'S', ZLCTES, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 2, W, 1, RW, BW, INFO );
         chkxer('ZGGES3', INFOT, NOUT, LERR, OK );
         NT = NT + 11

         // ZGGESX

         SRNAMT = 'ZGGESX'
         INFOT = 1
         zggesx('/', 'N', 'S', ZLCTSX, 'N', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zggesx('N', '/', 'S', ZLCTSX, 'N', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zggesx('V', 'V', '/', ZLCTSX, 'N', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zggesx('V', 'V', 'S', ZLCTSX, '/', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zggesx('V', 'V', 'S', ZLCTSX, 'B', -1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zggesx('V', 'V', 'S', ZLCTSX, 'B', 1, A, 0, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zggesx('V', 'V', 'S', ZLCTSX, 'B', 1, A, 1, B, 0, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 15
         zggesx('V', 'V', 'S', ZLCTSX, 'B', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 0, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 15
         zggesx('V', 'V', 'S', ZLCTSX, 'B', 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 17
         zggesx('V', 'V', 'S', ZLCTSX, 'B', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 0, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 17
         zggesx('V', 'V', 'S', ZLCTSX, 'B', 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 21
         zggesx('V', 'V', 'S', ZLCTSX, 'B', 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 2, RCE, RCV, W, 1, RW, IW, 1, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         INFOT = 24
         zggesx('V', 'V', 'S', ZLCTSX, 'V', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1, RCE, RCV, W, 32, RW, IW, 0, BW, INFO );
         chkxer('ZGGESX', INFOT, NOUT, LERR, OK );
         NT = NT + 13

         // ZGGEV

         SRNAMT = 'ZGGEV '
         INFOT = 1
         zggev('/', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zggev('N', '/', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zggev('V', 'V', -1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zggev('V', 'V', 1, A, 0, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zggev('V', 'V', 1, A, 1, B, 0, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zggev('N', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 0, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zggev('V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 1, U, 2, W, 1, RW, INFO );
         chkxer('ZGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 13
         zggev('V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 0, W, 1, RW, INFO );
         chkxer('ZGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 13
         zggev('V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV ', INFOT, NOUT, LERR, OK );
         INFOT = 15
         zggev('V', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 10

         // ZGGEV3

         SRNAMT = 'ZGGEV3'
         INFOT = 1
         zggev3('/', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zggev3('N', '/', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zggev3('V', 'V', -1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zggev3('V', 'V', 1, A, 0, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zggev3('V', 'V', 1, A, 1, B, 0, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zggev3('N', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 0, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zggev3('V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 1, U, 2, W, 1, RW, INFO );
         chkxer('ZGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 13
         zggev3('V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 0, W, 1, RW, INFO );
         chkxer('ZGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 13
         zggev3('V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV3', INFOT, NOUT, LERR, OK );
         INFOT = 15
         zggev3('V', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO );
         chkxer('ZGGEV3', INFOT, NOUT, LERR, OK );
         NT = NT + 10

         // ZGGEVX

         SRNAMT = 'ZGGEVX'
         INFOT = 1
         zggevx('/', 'N', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('ZGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zggevx('N', '/', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('ZGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zggevx('N', 'N', '/', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('ZGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zggevx('N', 'N', 'N', '/', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('ZGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zggevx('N', 'N', 'N', 'N', -1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('ZGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zggevx('N', 'N', 'N', 'N', 1, A, 0, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('ZGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zggevx('N', 'N', 'N', 'N', 1, A, 1, B, 0, ALPHA, BETA, Q, 1, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('ZGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 13
         zggevx('N', 'N', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 0, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('ZGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 13
         zggevx('N', 'V', 'N', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 1, U, 2, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('ZGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 15
         zggevx('N', 'N', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 0, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('ZGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 15
         zggevx('N', 'N', 'V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 1, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO );
         chkxer('ZGGEVX', INFOT, NOUT, LERR, OK );
         INFOT = 25
         zggevx('N', 'N', 'V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 2, ILO, IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 0, RW, IW, BW, INFO );
         chkxer('ZGGEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 12

         // ZTGEXC

         SRNAMT = 'ZTGEXC'
         INFOT = 3
         ztgexc( true , true , -1, A, 1, B, 1, Q, 1, Z, 1, IFST, ILST, INFO );
         chkxer('ZTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 5
         ztgexc( true , true , 1, A, 0, B, 1, Q, 1, Z, 1, IFST, ILST, INFO );
         chkxer('ZTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 7
         ztgexc( true , true , 1, A, 1, B, 0, Q, 1, Z, 1, IFST, ILST, INFO );
         chkxer('ZTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 9
         ztgexc( false , true , 1, A, 1, B, 1, Q, 0, Z, 1, IFST, ILST, INFO );
         chkxer('ZTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 9
         ztgexc( true , true , 1, A, 1, B, 1, Q, 0, Z, 1, IFST, ILST, INFO );
         chkxer('ZTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 11
         ztgexc( true , false , 1, A, 1, B, 1, Q, 1, Z, 0, IFST, ILST, INFO );
         chkxer('ZTGEXC', INFOT, NOUT, LERR, OK );
         INFOT = 11
         ztgexc( true , true , 1, A, 1, B, 1, Q, 1, Z, 0, IFST, ILST, INFO );
         chkxer('ZTGEXC', INFOT, NOUT, LERR, OK );
         NT = NT + 7

         // ZTGSEN

         SRNAMT = 'ZTGSEN'
         INFOT = 1
         ztgsen(-1, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('ZTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 5
         ztgsen(1, true , true , SEL, -1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('ZTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 7
         ztgsen(1, true , true , SEL, 1, A, 0, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('ZTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 9
         ztgsen(1, true , true , SEL, 1, A, 1, B, 0, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('ZTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 13
         ztgsen(1, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 0, Z, 1, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('ZTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 15
         ztgsen(1, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 0, M, TOLA, TOLB, RCV, W, 1, IW, 1, INFO );
         chkxer('ZTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 21
         ztgsen(3, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, -5, IW, 1, INFO );
         chkxer('ZTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 23
         ztgsen(0, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 0, INFO );
         chkxer('ZTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 23
         ztgsen(1, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 0, INFO );
         chkxer('ZTGSEN', INFOT, NOUT, LERR, OK );
         INFOT = 23
         ztgsen(5, true , true , SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA, TOLB, RCV, W, 20, IW, 1, INFO );
         chkxer('ZTGSEN', INFOT, NOUT, LERR, OK );
         NT = NT + 11

         // ZTGSNA

         SRNAMT = 'ZTGSNA'
         INFOT = 1
         ztgsna('/', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('ZTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ztgsna('B', '/', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('ZTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 4
         ztgsna('B', 'A', SEL, -1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('ZTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 6
         ztgsna('B', 'A', SEL, 1, A, 0, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('ZTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 8
         ztgsna('B', 'A', SEL, 1, A, 1, B, 0, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('ZTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 10
         ztgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 0, U, 1, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('ZTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 12
         ztgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 0, R1, R2, 1, M, W, 1, IW, INFO );
         chkxer('ZTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 15
         ztgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 0, M, W, 1, IW, INFO );
         chkxer('ZTGSNA', INFOT, NOUT, LERR, OK );
         INFOT = 18
         ztgsna('E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 0, IW, INFO );
         chkxer('ZTGSNA', INFOT, NOUT, LERR, OK );
         NT = NT + 9

         // ZTGSYL

         SRNAMT = 'ZTGSYL'
         INFOT = 1
         ztgsyl('/', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('ZTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 2
         ztgsyl('N', -1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('ZTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 3
         ztgsyl('N', 0, 0, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('ZTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 4
         ztgsyl('N', 0, 1, 0, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('ZTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 6
         ztgsyl('N', 0, 1, 1, A, 0, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('ZTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 8
         ztgsyl('N', 0, 1, 1, A, 1, B, 0, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('ZTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 10
         ztgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 0, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('ZTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 12
         ztgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 0, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('ZTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 14
         ztgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 0, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('ZTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 16
         ztgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 0, SCALE, DIF, W, 1, IW, INFO );
         chkxer('ZTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 20
         ztgsyl('N', 1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('ZTGSYL', INFOT, NOUT, LERR, OK );
         INFOT = 20
         ztgsyl('N', 2, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1, IW, INFO );
         chkxer('ZTGSYL', INFOT, NOUT, LERR, OK );
         NT = NT + 12
      }

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH, NT
      } else {
         WRITE( NOUT, FMT = 9998 )PATH
      }

 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits (', I3, ' tests done)' )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ', 'exits ***' )

      RETURN

      // End of ZERRGG

      }
