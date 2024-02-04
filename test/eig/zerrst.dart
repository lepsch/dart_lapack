      void zerrst(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

// =====================================================================

      // .. Parameters ..
      int                NMAX, LIW, LW;
      const              NMAX = 3, LIW = 12*NMAX, LW = 20*NMAX ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, INFO, J, M, N, NT;
      // ..
      // .. Local Arrays ..
      int                I1( NMAX ), I2( NMAX ), I3( NMAX ), IW( LIW );
      double             D( NMAX ), E( NMAX ), R( LW ), RW( LW ), X( NMAX );
      Complex         A( NMAX, NMAX ), C( NMAX, NMAX ), Q( NMAX, NMAX ), TAU( NMAX ), W( LW ), Z( NMAX, NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, ZHBEV, ZHBEVD, ZHBEVX, ZHBTRD, ZHEEV, ZHEEVD, ZHEEVR, ZHEEVX, ZHETRD, ZHPEV, ZHPEVD, ZHPEVX, ZHPTRD, ZPTEQR, ZSTEDC, ZSTEIN, ZSTEQR, ZUNGTR, ZUNMTR, ZUPGTR, ZUPMTR, ZHETD2, ZHEEVD_2STAGE, ZHEEVR_2STAGE, ZHEEVX_2STAGE, ZHEEV_2STAGE, ZHBEV_2STAGE, ZHBEVD_2STAGE, ZHBEVX_2STAGE, ZHETRD_2STAGE
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
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I, J] = 1.0 / (I+J).toDouble();
         } // 10
      } // 20
      for (J = 1; J <= NMAX; J++) { // 30
         D[J] = J.toDouble();
         E[J] = 0.0;
         I1[J] = J;
         I2[J] = J;
         TAU[J] = 1.0;
      } // 30
      OK = true;
      NT = 0;

      // Test error exits for the ST path.

      if ( LSAMEN( 2, C2, 'ST' ) ) {

         // ZHETRD

         SRNAMT = 'ZHETRD';
         INFOT = 1;
         zhetrd('/', 0, A, 1, D, E, TAU, W, 1, INFO );
         chkxer('ZHETRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhetrd('U', -1, A, 1, D, E, TAU, W, 1, INFO );
         chkxer('ZHETRD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zhetrd('U', 2, A, 1, D, E, TAU, W, 1, INFO );
         chkxer('ZHETRD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zhetrd('U', 0, A, 1, D, E, TAU, W, 0, INFO );
         chkxer('ZHETRD', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // ZHETD2

         SRNAMT = 'ZHETD2';
         INFOT = 1;
         zhetd2('/', 0, A, 1, D, E, TAU, INFO );
         chkxer('ZHETD2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhetd2('U', -1, A, 1, D, E, TAU, INFO );
         chkxer('ZHETD2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zhetd2('U', 2, A, 1, D, E, TAU, INFO );
         chkxer('ZHETD2', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // ZHETRD_2STAGE

         SRNAMT = 'ZHETRD_2STAGE';
         INFOT = 1;
         zhetrd_2stage('/', 'U', 0, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         zhetrd_2stage('H', 'U', 0, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhetrd_2stage('N', '/', 0, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhetrd_2stage('N', 'U', -1, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zhetrd_2stage('N', 'U', 2, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zhetrd_2stage('N', 'U', 0, A, 1, D, E, TAU,  C, 0, W, 1, INFO );
         chkxer('ZHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zhetrd_2stage('N', 'U', 0, A, 1, D, E, TAU,  C, 1, W, 0, INFO );
         chkxer('ZHETRD_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // ZHETRD_HE2HB

         SRNAMT = 'ZHETRD_HE2HB';
         INFOT = 1;
         zhetrd_he2hb('/', 0, 0, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('ZHETRD_HE2HB', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhetrd_he2hb('U', -1, 0, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('ZHETRD_HE2HB', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhetrd_he2hb('U', 0, -1, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('ZHETRD_HE2HB', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zhetrd_he2hb('U', 2, 0, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('ZHETRD_HE2HB', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zhetrd_he2hb('U', 0, 2, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('ZHETRD_HE2HB', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zhetrd_he2hb('U', 0, 0, A, 1, C, 1, TAU, W, 0, INFO );
         chkxer('ZHETRD_HE2HB', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // ZHETRD_HB2ST

         SRNAMT = 'ZHETRD_HB2ST';
         INFOT = 1;
         zhetrd_hb2st('/', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhetrd_hb2st('Y', '/', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhetrd_hb2st('Y', 'H', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhetrd_hb2st('Y', 'N', '/', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zhetrd_hb2st('Y', 'N', 'U', -1, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zhetrd_hb2st('Y', 'N', 'U', 0, -1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zhetrd_hb2st('Y', 'N', 'U', 0, 1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zhetrd_hb2st('Y', 'N', 'U', 0, 0, A, 1, D, E,  C, 0, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zhetrd_hb2st('Y', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 0, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // ZUNGTR

         SRNAMT = 'ZUNGTR';
         INFOT = 1;
         zungtr('/', 0, A, 1, TAU, W, 1, INFO );
         chkxer('ZUNGTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zungtr('U', -1, A, 1, TAU, W, 1, INFO );
         chkxer('ZUNGTR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zungtr('U', 2, A, 1, TAU, W, 1, INFO );
         chkxer('ZUNGTR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zungtr('U', 3, A, 3, TAU, W, 1, INFO );
         chkxer('ZUNGTR', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // ZUNMTR

         SRNAMT = 'ZUNMTR';
         INFOT = 1;
         zunmtr('/', 'U', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zunmtr('L', '/', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zunmtr('L', 'U', '/', 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zunmtr('L', 'U', 'N', -1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zunmtr('L', 'U', 'N', 0, -1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zunmtr('L', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('ZUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zunmtr('R', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zunmtr('L', 'U', 'N', 2, 0, A, 2, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zunmtr('L', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('ZUNMTR', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zunmtr('R', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('ZUNMTR', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // ZHPTRD

         SRNAMT = 'ZHPTRD';
         INFOT = 1;
         zhptrd('/', 0, A, D, E, TAU, INFO );
         chkxer('ZHPTRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhptrd('U', -1, A, D, E, TAU, INFO );
         chkxer('ZHPTRD', INFOT, NOUT, LERR, OK );
         NT = NT + 2;

         // ZUPGTR

         SRNAMT = 'ZUPGTR';
         INFOT = 1;
         zupgtr('/', 0, A, TAU, Z, 1, W, INFO );
         chkxer('ZUPGTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zupgtr('U', -1, A, TAU, Z, 1, W, INFO );
         chkxer('ZUPGTR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zupgtr('U', 2, A, TAU, Z, 1, W, INFO );
         chkxer('ZUPGTR', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // ZUPMTR

         SRNAMT = 'ZUPMTR';
         INFOT = 1;
         zupmtr('/', 'U', 'N', 0, 0, A, TAU, C, 1, W, INFO );
         chkxer('ZUPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zupmtr('L', '/', 'N', 0, 0, A, TAU, C, 1, W, INFO );
         chkxer('ZUPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zupmtr('L', 'U', '/', 0, 0, A, TAU, C, 1, W, INFO );
         chkxer('ZUPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zupmtr('L', 'U', 'N', -1, 0, A, TAU, C, 1, W, INFO );
         chkxer('ZUPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zupmtr('L', 'U', 'N', 0, -1, A, TAU, C, 1, W, INFO );
         chkxer('ZUPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zupmtr('L', 'U', 'N', 2, 0, A, TAU, C, 1, W, INFO );
         chkxer('ZUPMTR', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // ZPTEQR

         SRNAMT = 'ZPTEQR';
         INFOT = 1;
         zpteqr('/', 0, D, E, Z, 1, RW, INFO );
         chkxer('ZPTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zpteqr('N', -1, D, E, Z, 1, RW, INFO );
         chkxer('ZPTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zpteqr('V', 2, D, E, Z, 1, RW, INFO );
         chkxer('ZPTEQR', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // ZSTEIN

         SRNAMT = 'ZSTEIN';
         INFOT = 1;
         zstein(-1, D, E, 0, X, I1, I2, Z, 1, RW, IW, I3, INFO );
         chkxer('ZSTEIN', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zstein(0, D, E, -1, X, I1, I2, Z, 1, RW, IW, I3, INFO );
         chkxer('ZSTEIN', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zstein(0, D, E, 1, X, I1, I2, Z, 1, RW, IW, I3, INFO );
         chkxer('ZSTEIN', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zstein(2, D, E, 0, X, I1, I2, Z, 1, RW, IW, I3, INFO );
         chkxer('ZSTEIN', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // ZSTEQR

         SRNAMT = 'ZSTEQR';
         INFOT = 1;
         zsteqr('/', 0, D, E, Z, 1, RW, INFO );
         chkxer('ZSTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsteqr('N', -1, D, E, Z, 1, RW, INFO );
         chkxer('ZSTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zsteqr('V', 2, D, E, Z, 1, RW, INFO );
         chkxer('ZSTEQR', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // ZSTEDC

         SRNAMT = 'ZSTEDC';
         INFOT = 1;
         zstedc('/', 0, D, E, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zstedc('N', -1, D, E, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zstedc('V', 2, D, E, Z, 1, W, 4, RW, 23, IW, 28, INFO );
         chkxer('ZSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zstedc('N', 2, D, E, Z, 1, W, 0, RW, 1, IW, 1, INFO );
         chkxer('ZSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zstedc('V', 2, D, E, Z, 2, W, 0, RW, 23, IW, 28, INFO );
         chkxer('ZSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zstedc('N', 2, D, E, Z, 1, W, 1, RW, 0, IW, 1, INFO );
         chkxer('ZSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zstedc('I', 2, D, E, Z, 2, W, 1, RW, 1, IW, 12, INFO );
         chkxer('ZSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zstedc('V', 2, D, E, Z, 2, W, 4, RW, 1, IW, 28, INFO );
         chkxer('ZSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zstedc('N', 2, D, E, Z, 1, W, 1, RW, 1, IW, 0, INFO );
         chkxer('ZSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zstedc('I', 2, D, E, Z, 2, W, 1, RW, 23, IW, 0, INFO );
         chkxer('ZSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zstedc('V', 2, D, E, Z, 2, W, 4, RW, 23, IW, 0, INFO );
         chkxer('ZSTEDC', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // ZHEEVD

         SRNAMT = 'ZHEEVD';
         INFOT = 1;
         zheevd('/', 'U', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zheevd('N', '/', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zheevd('N', 'U', -1, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zheevd('N', 'U', 2, A, 1, X, W, 3, RW, 2, IW, 1, INFO );
         chkxer('ZHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zheevd('N', 'U', 1, A, 1, X, W, 0, RW, 1, IW, 1, INFO );
         chkxer('ZHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zheevd('N', 'U', 2, A, 2, X, W, 2, RW, 2, IW, 1, INFO );
         chkxer('ZHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zheevd('V', 'U', 2, A, 2, X, W, 3, RW, 25, IW, 12, INFO );
         chkxer('ZHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zheevd('N', 'U', 1, A, 1, X, W, 1, RW, 0, IW, 1, INFO );
         chkxer('ZHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zheevd('N', 'U', 2, A, 2, X, W, 3, RW, 1, IW, 1, INFO );
         chkxer('ZHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zheevd('V', 'U', 2, A, 2, X, W, 8, RW, 18, IW, 12, INFO );
         chkxer('ZHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zheevd('N', 'U', 1, A, 1, X, W, 1, RW, 1, IW, 0, INFO );
         chkxer('ZHEEVD', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zheevd('V', 'U', 2, A, 2, X, W, 8, RW, 25, IW, 11, INFO );
         chkxer('ZHEEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // ZHEEVD_2STAGE

         SRNAMT = 'ZHEEVD_2STAGE';
         INFOT = 1;
         zheevd_2stage('/', 'U', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         zheevd_2stage('V', 'U', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zheevd_2stage('N', '/', 0, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zheevd_2stage('N', 'U', -1, A, 1, X, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zheevd_2stage('N', 'U', 2, A, 1, X, W, 3, RW, 2, IW, 1, INFO );
         chkxer('ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zheevd_2stage('N', 'U', 1, A, 1, X, W, 0, RW, 1, IW, 1, INFO );
         chkxer('ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zheevd_2stage('N', 'U', 2, A, 2, X, W, 2, RW, 2, IW, 1, INFO );
         chkxer('ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 8
          // CALL ZHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 3,
      // $                            RW, 25, IW, 12, INFO )
          // CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 10;
         zheevd_2stage('N', 'U', 1, A, 1, X, W, 1, RW, 0, IW, 1, INFO );
         chkxer('ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zheevd_2stage('N', 'U', 2, A, 2, X, W, 25, RW, 1, IW, 1, INFO );
         chkxer('ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 10
          // CALL ZHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 8,
      // $                            RW, 18, IW, 12, INFO )
          // CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 12;
         zheevd_2stage('N', 'U', 1, A, 1, X, W, 1, RW, 1, IW, 0, INFO );
         chkxer('ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
          // CALL ZHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 8,
      // $                            RW, 25, IW, 11, INFO )
          // CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
         NT = NT + 10;

         // ZHEEV

         SRNAMT = 'ZHEEV ';
         INFOT = 1;
         zheev('/', 'U', 0, A, 1, X, W, 1, RW, INFO );
         chkxer('ZHEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zheev('N', '/', 0, A, 1, X, W, 1, RW, INFO );
         chkxer('ZHEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zheev('N', 'U', -1, A, 1, X, W, 1, RW, INFO );
         chkxer('ZHEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zheev('N', 'U', 2, A, 1, X, W, 3, RW, INFO );
         chkxer('ZHEEV ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zheev('N', 'U', 2, A, 2, X, W, 2, RW, INFO );
         chkxer('ZHEEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 5;

         // ZHEEV_2STAGE

         SRNAMT = 'ZHEEV_2STAGE ';
         INFOT = 1;
         zheev_2stage('/', 'U', 0, A, 1, X, W, 1, RW, INFO );
         chkxer('ZHEEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         zheev_2stage('V', 'U', 0, A, 1, X, W, 1, RW, INFO );
         chkxer('ZHEEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zheev_2stage('N', '/', 0, A, 1, X, W, 1, RW, INFO );
         chkxer('ZHEEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zheev_2stage('N', 'U', -1, A, 1, X, W, 1, RW, INFO );
         chkxer('ZHEEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zheev_2stage('N', 'U', 2, A, 1, X, W, 3, RW, INFO );
         chkxer('ZHEEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zheev_2stage('N', 'U', 2, A, 2, X, W, 2, RW, INFO );
         chkxer('ZHEEV_2STAGE ', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // ZHEEVX

         SRNAMT = 'ZHEEVX';
         INFOT = 1;
         zheevx('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('ZHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zheevx('V', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('ZHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zheevx('V', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         INFOT = 4;
         zheevx('V', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('ZHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zheevx('V', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, 3, RW, IW, I3, INFO );
         chkxer('ZHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zheevx('V', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('ZHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zheevx('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('ZHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zheevx('V', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 2, W, 3, RW, IW, I3, INFO );
         chkxer('ZHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zheevx('V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 3, RW, IW, I3, INFO );
         chkxer('ZHEEVX', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         zheevx('V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, 2, RW, IW, I1, INFO );
         chkxer('ZHEEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // ZHEEVX_2STAGE

         SRNAMT = 'ZHEEVX_2STAGE';
         INFOT = 1;
         zheevx_2stage('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         zheevx_2stage('V', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zheevx_2stage('N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zheevx_2stage('N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         INFOT = 4;
         zheevx_2stage('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zheevx_2stage('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, 3, RW, IW, I3, INFO );
         chkxer('ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zheevx_2stage('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zheevx_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, RW, IW, I3, INFO );
         chkxer('ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zheevx_2stage('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 2, W, 3, RW, IW, I3, INFO );
         chkxer('ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zheevx_2stage('N', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 0, W, 3, RW, IW, I3, INFO );
         chkxer('ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         zheevx_2stage('N', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, 0, RW, IW, I1, INFO );
         chkxer('ZHEEVX_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // ZHEEVR

         SRNAMT = 'ZHEEVR';
         N = 1;
         INFOT = 1;
         zheevr('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zheevr('V', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zheevr('V', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zheevr('V', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zheevr('V', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zheevr('V', 'V', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zheevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 10;

         zheevr('V', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zheevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 0, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         zheevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 0, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         zheevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 0, IW( 2*N-1 ), 10*N, INFO );
         chkxer('ZHEEVR', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         zheevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW, 0, INFO );
         chkxer('ZHEEVR', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // ZHEEVR_2STAGE

         SRNAMT = 'ZHEEVR_2STAGE';
         N = 1;
         INFOT = 1;
         zheevr_2stage('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         zheevr_2stage('V', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zheevr_2stage('N', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zheevr_2stage('N', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zheevr_2stage('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zheevr_2stage('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zheevr_2stage('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zheevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zheevr_2stage('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, R, Z, 1, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zheevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 0, IW, Q, 2*N, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         zheevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 0, RW, 24*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         zheevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, RW, 0, IW( 2*N-1 ), 10*N, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 22;
         zheevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, RW, 24*N, IW, 0, INFO );
         chkxer('ZHEEVR_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 13;

         // ZHPEVD

         SRNAMT = 'ZHPEVD';
         INFOT = 1;
         zhpevd('/', 'U', 0, A, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhpevd('N', '/', 0, A, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhpevd('N', 'U', -1, A, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zhpevd('V', 'U', 2, A, X, Z, 1, W, 4, RW, 25, IW, 12, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zhpevd('N', 'U', 1, A, X, Z, 1, W, 0, RW, 1, IW, 1, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zhpevd('N', 'U', 2, A, X, Z, 2, W, 1, RW, 2, IW, 1, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zhpevd('V', 'U', 2, A, X, Z, 2, W, 2, RW, 25, IW, 12, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zhpevd('N', 'U', 1, A, X, Z, 1, W, 1, RW, 0, IW, 1, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zhpevd('N', 'U', 2, A, X, Z, 2, W, 2, RW, 1, IW, 1, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zhpevd('V', 'U', 2, A, X, Z, 2, W, 4, RW, 18, IW, 12, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zhpevd('N', 'U', 1, A, X, Z, 1, W, 1, RW, 1, IW, 0, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zhpevd('N', 'U', 2, A, X, Z, 2, W, 2, RW, 2, IW, 0, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zhpevd('V', 'U', 2, A, X, Z, 2, W, 4, RW, 25, IW, 2, INFO );
         chkxer('ZHPEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 13;

         // ZHPEV

         SRNAMT = 'ZHPEV ';
         INFOT = 1;
         zhpev('/', 'U', 0, A, X, Z, 1, W, RW, INFO );
         chkxer('ZHPEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhpev('N', '/', 0, A, X, Z, 1, W, RW, INFO );
         chkxer('ZHPEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhpev('N', 'U', -1, A, X, Z, 1, W, RW, INFO );
         chkxer('ZHPEV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zhpev('V', 'U', 2, A, X, Z, 1, W, RW, INFO );
         chkxer('ZHPEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // ZHPEVX

         SRNAMT = 'ZHPEVX';
         INFOT = 1;
         zhpevx('/', 'A', 'U', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhpevx('V', '/', 'U', 0, A, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhpevx('V', 'A', '/', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zhpevx('V', 'A', 'U', -1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zhpevx('V', 'V', 'U', 1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zhpevx('V', 'I', 'U', 1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zhpevx('V', 'I', 'U', 2, A, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 2, W, RW, IW, I3, INFO );
         chkxer('ZHPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         zhpevx('V', 'A', 'U', 2, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHPEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

      // Test error exits for the HB path.

      } else if ( LSAMEN( 2, C2, 'HB' ) ) {

         // ZHBTRD

         SRNAMT = 'ZHBTRD';
         INFOT = 1;
         zhbtrd('/', 'U', 0, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('ZHBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhbtrd('N', '/', 0, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('ZHBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhbtrd('N', 'U', -1, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('ZHBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zhbtrd('N', 'U', 0, -1, A, 1, D, E, Z, 1, W, INFO );
         chkxer('ZHBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zhbtrd('N', 'U', 1, 1, A, 1, D, E, Z, 1, W, INFO );
         chkxer('ZHBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zhbtrd('V', 'U', 2, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('ZHBTRD', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // ZHETRD_HB2ST

         SRNAMT = 'ZHETRD_HB2ST';
         INFOT = 1;
         zhetrd_hb2st('/', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhetrd_hb2st('N', '/', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhetrd_hb2st('N', 'H', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhetrd_hb2st('N', 'N', '/', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zhetrd_hb2st('N', 'N', 'U', -1, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zhetrd_hb2st('N', 'N', 'U', 0, -1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zhetrd_hb2st('N', 'N', 'U', 0, 1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zhetrd_hb2st('N', 'N', 'U', 0, 0, A, 1, D, E,  C, 0, W, 1, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zhetrd_hb2st('N', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 0, INFO );
         chkxer('ZHETRD_HB2ST', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // ZHBEVD

         SRNAMT = 'ZHBEVD';
         INFOT = 1;
         zhbevd('/', 'U', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhbevd('N', '/', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhbevd('N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zhbevd('N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zhbevd('N', 'U', 2, 1, A, 1, X, Z, 1, W, 2, RW, 2, IW, 1, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zhbevd('V', 'U', 2, 1, A, 2, X, Z, 1, W, 8, RW, 25, IW, 12, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zhbevd('N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zhbevd('N', 'U', 2, 1, A, 2, X, Z, 2, W, 1, RW, 2, IW, 1, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zhbevd('V', 'U', 2, 1, A, 2, X, Z, 2, W, 2, RW, 25, IW, 12, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zhbevd('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 0, IW, 1, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zhbevd('N', 'U', 2, 1, A, 2, X, Z, 2, W, 2, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zhbevd('V', 'U', 2, 1, A, 2, X, Z, 2, W, 8, RW, 2, IW, 12, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zhbevd('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 0, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zhbevd('N', 'U', 2, 1, A, 2, X, Z, 2, W, 2, RW, 2, IW, 0, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zhbevd('V', 'U', 2, 1, A, 2, X, Z, 2, W, 8, RW, 25, IW, 2, INFO );
         chkxer('ZHBEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 15;

         // ZHBEVD_2STAGE

         SRNAMT = 'ZHBEVD_2STAGE';
         INFOT = 1;
         zhbevd_2stage('/', 'U', 0, 0, A, 1, X, Z, 1,  W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         zhbevd_2stage('V', 'U', 0, 0, A, 1, X, Z, 1,  W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhbevd_2stage('N', '/', 0, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhbevd_2stage('N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zhbevd_2stage('N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zhbevd_2stage('N', 'U', 2, 1, A, 1, X, Z, 1, W, 2, RW, 2, IW, 1, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zhbevd_2stage('N', 'U', 2, 1, A, 2, X, Z, 0, W, 8, RW, 25, IW, 12, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zhbevd_2stage('N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zhbevd_2stage('N', 'U', 2, 1, A, 2, X, Z, 2, W, 1, RW, 2, IW, 1, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 11
          // CALL ZHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
      // $                         W, 2, RW, 25, IW, 12, INFO )
          // CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 13;
         zhbevd_2stage('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 0, IW, 1, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zhbevd_2stage('N', 'U', 2, 1, A, 2, X, Z, 2, W, 25, RW, 1, IW, 1, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 13
          // CALL ZHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
      // $                          W, 25, RW, 2, IW, 12, INFO )
          // CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 15;
         zhbevd_2stage('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, RW, 1, IW, 0, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         zhbevd_2stage('N', 'U', 2, 1, A, 2, X, Z, 2, W, 25, RW, 2, IW, 0, INFO );
         chkxer('ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 15
          // CALL ZHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
      // $                          W, 25, RW, 25, IW, 2, INFO )
          // CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         NT = NT + 13;

         // ZHBEV

         SRNAMT = 'ZHBEV ';
         INFOT = 1;
         zhbev('/', 'U', 0, 0, A, 1, X, Z, 1, W, RW, INFO );
         chkxer('ZHBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhbev('N', '/', 0, 0, A, 1, X, Z, 1, W, RW, INFO );
         chkxer('ZHBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhbev('N', 'U', -1, 0, A, 1, X, Z, 1, W, RW, INFO );
         chkxer('ZHBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zhbev('N', 'U', 0, -1, A, 1, X, Z, 1, W, RW, INFO );
         chkxer('ZHBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zhbev('N', 'U', 2, 1, A, 1, X, Z, 1, W, RW, INFO );
         chkxer('ZHBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zhbev('V', 'U', 2, 0, A, 1, X, Z, 1, W, RW, INFO );
         chkxer('ZHBEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // ZHBEV_2STAGE

         SRNAMT = 'ZHBEV_2STAGE ';
         INFOT = 1;
         zhbev_2stage('/', 'U', 0, 0, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         zhbev_2stage('V', 'U', 0, 0, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhbev_2stage('N', '/', 0, 0, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhbev_2stage('N', 'U', -1, 0, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zhbev_2stage('N', 'U', 0, -1, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zhbev_2stage('N', 'U', 2, 1, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zhbev_2stage('N', 'U', 2, 0, A, 1, X, Z, 0, W, 0, RW, INFO );
         chkxer('ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zhbev_2stage('N', 'U', 2, 0, A, 1, X, Z, 1, W, 0, RW, INFO );
         chkxer('ZHBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

         // ZHBEVX

         SRNAMT = 'ZHBEVX';
         INFOT = 1;
         zhbevx('/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhbevx('V', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhbevx('V', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         INFOT = 4;
         zhbevx('V', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zhbevx('V', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zhbevx('V', 'A', 'U', 2, 1, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, RW, IW, I3, INFO );
         chkxer('ZHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zhbevx('V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, RW, IW, I3, INFO );
         chkxer('ZHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zhbevx('V', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zhbevx('V', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zhbevx('V', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         zhbevx('V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, RW, IW, I3, INFO );
         chkxer('ZHBEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // ZHBEVX_2STAGE

         SRNAMT = 'ZHBEVX_2STAGE';
         INFOT = 1;
         zhbevx_2stage('/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         INFOT = 1;
         zhbevx_2stage('V', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zhbevx_2stage('N', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zhbevx_2stage('N', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         INFOT = 4;
         zhbevx_2stage('N', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zhbevx_2stage('N', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zhbevx_2stage('N', 'A', 'U', 2, 1, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, 0, RW, IW, I3, INFO );
         chkxer('ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 9
          // CALL ZHBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 1,
      // $                       0.0, 0.0, 0, 0, 0.0,
      // $                       M, X, Z, 2, W, 0, RW, IW, I3, INFO )
          // CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 11;
         zhbevx_2stage('N', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zhbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         zhbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         zhbevx_2stage('N', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 0, W, 0, RW, IW, I3, INFO );
         chkxer('ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         zhbevx_2stage('N', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, RW, IW, I3, INFO );
         chkxer('ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 12;
      }

      // Print a summary line.

      if ( OK ) {
         WRITE( NOUT, FMT = 9999 )PATH, NT;
      } else {
         WRITE( NOUT, FMT = 9998 )PATH;
      }

 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits', ' (', I3, ' tests done)' );
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ', 'exits ***' );

      return;
      }
