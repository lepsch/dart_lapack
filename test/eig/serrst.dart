      void serrst(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

// =====================================================================

      // NMAX has to be at least 3 or LIW may be too small
      // .. Parameters ..
      int                NMAX, LIW, LW;
      const              NMAX = 3, LIW = 12*NMAX, LW = 20*NMAX ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, INFO, J, M, N, NSPLIT, NT;
      // ..
      // .. Local Arrays ..
      int                I1( NMAX ), I2( NMAX ), I3( NMAX ), IW( LIW );
      REAL               A( NMAX, NMAX ), C( NMAX, NMAX ), D( NMAX ), E( NMAX ), Q( NMAX, NMAX ), R( NMAX ), TAU( NMAX ), W( LW ), X( NMAX ), Z( NMAX, NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, SOPGTR, SOPMTR, SORGTR, SORMTR, SPTEQR, SSBEV, SSBEVD, SSBEVX, SSBTRD, SSPEV, SSPEVD, SSPEVX, SSPTRD, SSTEBZ, SSTEDC, SSTEIN, SSTEQR, SSTERF, SSTEV, SSTEVD, SSTEVR, SSTEVX, SSYEV, SSYEVD, SSYEVR, SSYEVX, SSYTRD, SSYTD2, SSYEVD_2STAGE, SSYEVR_2STAGE, SSYEVX_2STAGE, SSYEV_2STAGE, SSBEV_2STAGE, SSBEVD_2STAGE, SSBEVX_2STAGE, SSYTRD_2STAGE, SSYTRD_SY2SB, SSYTRD_SB2ST
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
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = 1. / REAL( I+J );
         } // 10
      } // 20
      for (J = 1; J <= NMAX; J++) { // 30
         D( J ) = REAL( J );
         E( J ) = 0.0;
         I1( J ) = J;
         I2( J ) = J;
         TAU( J ) = 1.;
      } // 30
      OK = true;
      NT = 0;

      // Test error exits for the ST path.

      if ( LSAMEN( 2, C2, 'ST' ) ) {

         // SSYTRD

         SRNAMT = 'SSYTRD';
         INFOT = 1;
         ssytrd('/', 0, A, 1, D, E, TAU, W, 1, INFO );
         chkxer('SSYTRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrd('U', -1, A, 1, D, E, TAU, W, 1, INFO );
         chkxer('SSYTRD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytrd('U', 2, A, 1, D, E, TAU, W, 1, INFO );
         chkxer('SSYTRD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ssytrd('U', 0, A, 1, D, E, TAU, W, 0, INFO );
         chkxer('SSYTRD', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // SSYTD2

         SRNAMT = 'SSYTD2';
         INFOT = 1;
         ssytd2('/', 0, A, 1, D, E, TAU, INFO );
         chkxer('SSYTD2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytd2('U', -1, A, 1, D, E, TAU, INFO );
         chkxer('SSYTD2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytd2('U', 2, A, 1, D, E, TAU, INFO );
         chkxer('SSYTD2', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // SSYTRD_2STAGE

         SRNAMT = 'SSYTRD_2STAGE';
         INFOT = 1;
         ssytrd_2stage('/', 'U', 0, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         ssytrd_2stage('H', 'U', 0, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrd_2stage('N', '/', 0, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssytrd_2stage('N', 'U', -1, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssytrd_2stage('N', 'U', 2, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssytrd_2stage('N', 'U', 0, A, 1, D, E, TAU,  C, 0, W, 1, INFO );
         chkxer('SSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         ssytrd_2stage('N', 'U', 0, A, 1, D, E, TAU,  C, 1, W, 0, INFO );
         chkxer('SSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // SSYTRD_SY2SB

         SRNAMT = 'SSYTRD_SY2SB';
         INFOT = 1;
         ssytrd_sy2sb('/', 0, 0, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('SSYTRD_SY2SB', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrd_sy2sb('U', -1, 0, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('SSYTRD_SY2SB', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssytrd_sy2sb('U', 0, -1, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('SSYTRD_SY2SB', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssytrd_sy2sb('U', 2, 0, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('SSYTRD_SY2SB', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ssytrd_sy2sb('U', 0, 2, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('SSYTRD_SY2SB', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssytrd_sy2sb('U', 0, 0, A, 1, C, 1, TAU, W, 0, INFO );
         chkxer('SSYTRD_SY2SB', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // SSYTRD_SB2ST

         SRNAMT = 'SSYTRD_SB2ST';
         INFOT = 1;
         ssytrd_sb2st('/', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrd_sb2st('Y', '/', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrd_sb2st('Y', 'H', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssytrd_sb2st('Y', 'N', '/', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytrd_sb2st('Y', 'N', 'U', -1, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssytrd_sb2st('Y', 'N', 'U', 0, -1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ssytrd_sb2st('Y', 'N', 'U', 0, 1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ssytrd_sb2st('Y', 'N', 'U', 0, 0, A, 1, D, E,  C, 0, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         ssytrd_sb2st('Y', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 0, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // SORGTR

         SRNAMT = 'SORGTR';
         INFOT = 1;
         sorgtr('/', 0, A, 1, TAU, W, 1, INFO );
         chkxer('SORGTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sorgtr('U', -1, A, 1, TAU, W, 1, INFO );
         chkxer('SORGTR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sorgtr('U', 2, A, 1, TAU, W, 1, INFO );
         chkxer('SORGTR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sorgtr('U', 3, A, 3, TAU, W, 1, INFO );
         chkxer('SORGTR', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // SORMTR

         SRNAMT = 'SORMTR';
         INFOT = 1;
         sormtr('/', 'U', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sormtr('L', '/', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sormtr('L', 'U', '/', 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sormtr('L', 'U', 'N', -1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sormtr('L', 'U', 'N', 0, -1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sormtr('L', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('SORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sormtr('R', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sormtr('L', 'U', 'N', 2, 0, A, 2, TAU, C, 1, W, 1, INFO );
         chkxer('SORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sormtr('L', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('SORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         sormtr('R', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('SORMTR', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // SSPTRD

         SRNAMT = 'SSPTRD';
         INFOT = 1;
         ssptrd('/', 0, A, D, E, TAU, INFO );
         chkxer('SSPTRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssptrd('U', -1, A, D, E, TAU, INFO );
         chkxer('SSPTRD', INFOT, NOUT, LERR, OK );
         NT = NT + 2;

         // SOPGTR

         SRNAMT = 'SOPGTR';
         INFOT = 1;
         sopgtr('/', 0, A, TAU, Z, 1, W, INFO );
         chkxer('SOPGTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sopgtr('U', -1, A, TAU, Z, 1, W, INFO );
         chkxer('SOPGTR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sopgtr('U', 2, A, TAU, Z, 1, W, INFO );
         chkxer('SOPGTR', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // SOPMTR

         SRNAMT = 'SOPMTR';
         INFOT = 1;
         sopmtr('/', 'U', 'N', 0, 0, A, TAU, C, 1, W, INFO );
         chkxer('SOPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sopmtr('L', '/', 'N', 0, 0, A, TAU, C, 1, W, INFO );
         chkxer('SOPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sopmtr('L', 'U', '/', 0, 0, A, TAU, C, 1, W, INFO );
         chkxer('SOPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sopmtr('L', 'U', 'N', -1, 0, A, TAU, C, 1, W, INFO );
         chkxer('SOPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sopmtr('L', 'U', 'N', 0, -1, A, TAU, C, 1, W, INFO );
         chkxer('SOPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sopmtr('L', 'U', 'N', 2, 0, A, TAU, C, 1, W, INFO );
         chkxer('SOPMTR', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // SPTEQR

         SRNAMT = 'SPTEQR';
         INFOT = 1;
         spteqr('/', 0, D, E, Z, 1, W, INFO );
         chkxer('SPTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         spteqr('N', -1, D, E, Z, 1, W, INFO );
         chkxer('SPTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         spteqr('V', 2, D, E, Z, 1, W, INFO );
         chkxer('SPTEQR', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // SSTEBZ

         SRNAMT = 'SSTEBZ';
         INFOT = 1;
         sstebz('/', 'E', 0, 0.0, 1.0, 1, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('SSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sstebz('A', '/', 0, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('SSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sstebz('A', 'E', -1, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('SSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sstebz('V', 'E', 0, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('SSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sstebz('I', 'E', 0, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('SSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sstebz('I', 'E', 1, 0.0, 0.0, 2, 1, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('SSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sstebz('I', 'E', 1, 0.0, 0.0, 1, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('SSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sstebz('I', 'E', 1, 0.0, 0.0, 1, 2, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('SSTEBZ', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

         // SSTEIN

         SRNAMT = 'SSTEIN';
         INFOT = 1;
         sstein(-1, D, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEIN', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sstein(0, D, E, -1, X, I1, I2, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEIN', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         sstein(0, D, E, 1, X, I1, I2, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEIN', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sstein(2, D, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEIN', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // SSTEQR

         SRNAMT = 'SSTEQR';
         INFOT = 1;
         ssteqr('/', 0, D, E, Z, 1, W, INFO );
         chkxer('SSTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssteqr('N', -1, D, E, Z, 1, W, INFO );
         chkxer('SSTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssteqr('V', 2, D, E, Z, 1, W, INFO );
         chkxer('SSTEQR', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // SSTERF

         SRNAMT = 'SSTERF';
         INFOT = 1;
         ssterf(-1, D, E, INFO );
         chkxer('SSTERF', INFOT, NOUT, LERR, OK );
         NT = NT + 1;

         // SSTEDC

         SRNAMT = 'SSTEDC';
         INFOT = 1;
         sstedc('/', 0, D, E, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sstedc('N', -1, D, E, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sstedc('V', 2, D, E, Z, 1, W, 23, IW, 28, INFO );
         chkxer('SSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sstedc('N', 1, D, E, Z, 1, W, 0, IW, 1, INFO );
         chkxer('SSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sstedc('I', 2, D, E, Z, 2, W, 0, IW, 12, INFO );
         chkxer('SSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sstedc('V', 2, D, E, Z, 2, W, 0, IW, 28, INFO );
         chkxer('SSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sstedc('N', 1, D, E, Z, 1, W, 1, IW, 0, INFO );
         chkxer('SSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sstedc('I', 2, D, E, Z, 2, W, 19, IW, 0, INFO );
         chkxer('SSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sstedc('V', 2, D, E, Z, 2, W, 23, IW, 0, INFO );
         chkxer('SSTEDC', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // SSTEVD

         SRNAMT = 'SSTEVD';
         INFOT = 1;
         sstevd('/', 0, D, E, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSTEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sstevd('N', -1, D, E, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSTEVD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sstevd('V', 2, D, E, Z, 1, W, 19, IW, 12, INFO );
         chkxer('SSTEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sstevd('N', 1, D, E, Z, 1, W, 0, IW, 1, INFO );
         chkxer('SSTEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sstevd('V', 2, D, E, Z, 2, W, 12, IW, 12, INFO );
         chkxer('SSTEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sstevd('N', 0, D, E, Z, 1, W, 1, IW, 0, INFO );
         chkxer('SSTEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         sstevd('V', 2, D, E, Z, 2, W, 19, IW, 11, INFO );
         chkxer('SSTEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // SSTEV

         SRNAMT = 'SSTEV ';
         INFOT = 1;
         sstev('/', 0, D, E, Z, 1, W, INFO );
         chkxer('SSTEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sstev('N', -1, D, E, Z, 1, W, INFO );
         chkxer('SSTEV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         sstev('V', 2, D, E, Z, 1, W, INFO );
         chkxer('SSTEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // SSTEVX

         SRNAMT = 'SSTEVX';
         INFOT = 1;
         sstevx('/', 'A', 0, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sstevx('N', '/', 0, D, E, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sstevx('N', 'A', -1, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sstevx('N', 'V', 1, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sstevx('N', 'I', 1, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sstevx('N', 'I', 1, D, E, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sstevx('N', 'I', 2, D, E, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sstevx('N', 'I', 1, D, E, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         sstevx('V', 'A', 2, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSTEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // SSTEVR

         N = 1;
         SRNAMT = 'SSTEVR';
         INFOT = 1;
         sstevr('/', 'A', 0, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sstevr('V', '/', 0, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sstevr('V', 'A', -1, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sstevr('V', 'V', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sstevr('V', 'I', 1, D, E, 0.0, 0.0, 0, 1, 0.0, M, W, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         N = 2;
         sstevr('V', 'I', 2, D, E, 0.0, 0.0, 2, 1, 0.0, M, W, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         N = 1;
         sstevr('V', 'I', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, W, Z, 0, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         sstevr('V', 'I', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, W, Z, 1, IW, X, 20*N-1, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 19;
         sstevr('V', 'I', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, W, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N-1, INFO );
         chkxer('SSTEVR', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // SSYEVD

         SRNAMT = 'SSYEVD';
         INFOT = 1;
         ssyevd('/', 'U', 0, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('SSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssyevd('N', '/', 0, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('SSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssyevd('N', 'U', -1, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('SSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssyevd('N', 'U', 2, A, 1, X, W, 3, IW, 1, INFO );
         chkxer('SSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssyevd('N', 'U', 1, A, 1, X, W, 0, IW, 1, INFO );
         chkxer('SSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssyevd('N', 'U', 2, A, 2, X, W, 4, IW, 1, INFO );
         chkxer('SSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssyevd('V', 'U', 2, A, 2, X, W, 20, IW, 12, INFO );
         chkxer('SSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssyevd('N', 'U', 1, A, 1, X, W, 1, IW, 0, INFO );
         chkxer('SSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssyevd('N', 'U', 2, A, 2, X, W, 5, IW, 0, INFO );
         chkxer('SSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssyevd('V', 'U', 2, A, 2, X, W, 27, IW, 11, INFO );
         chkxer('SSYEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // SSYEVD_2STAGE

         SRNAMT = 'SSYEVD_2STAGE';
         INFOT = 1;
         ssyevd_2stage('/', 'U', 0, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('SSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         ssyevd_2stage('V', 'U', 0, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('SSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssyevd_2stage('N', '/', 0, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('SSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssyevd_2stage('N', 'U', -1, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('SSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssyevd_2stage('N', 'U', 2, A, 1, X, W, 3, IW, 1, INFO );
         chkxer('SSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssyevd_2stage('N', 'U', 1, A, 1, X, W, 0, IW, 1, INFO );
         chkxer('SSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssyevd_2stage('N', 'U', 2, A, 2, X, W, 4, IW, 1, INFO );
         chkxer('SSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 8
          // CALL SSYEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 20, IW, 12, INFO )
          // CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 10;
         ssyevd_2stage('N', 'U', 1, A, 1, X, W, 1, IW, 0, INFO );
         chkxer('SSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssyevd_2stage('N', 'U', 2, A, 2, X, W, 25, IW, 0, INFO );
         chkxer('SSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 10
          // CALL SSYEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 27, IW, 11, INFO )
          // CALL CHKXER( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         NT = NT + 9;

         // SSYEVR

         SRNAMT = 'SSYEVR';
         N = 1;
         INFOT = 1;
         ssyevr('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssyevr('V', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssyevr('V', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssyevr('V', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssyevr('V', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssyevr('V', 'V', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ssyevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 10;

         ssyevr('V', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         ssyevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 0, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         ssyevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 0, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         ssyevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 0, INFO );
         chkxer('SSYEVR', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // SSYEVR_2STAGE

         SRNAMT = 'SSYEVR_2STAGE';
         N = 1;
         INFOT = 1;
         ssyevr_2stage('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         ssyevr_2stage('V', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssyevr_2stage('N', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssyevr_2stage('N', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssyevr_2stage('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssyevr_2stage('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssyevr_2stage('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ssyevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssyevr_2stage('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         ssyevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 0, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         ssyevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 0, IW( 2*N+1 ), 10*N, INFO );
         chkxer('SSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         ssyevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 0, INFO );
         chkxer('SSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // SSYEV

         SRNAMT = 'SSYEV ';
         INFOT = 1;
         ssyev('/', 'U', 0, A, 1, X, W, 1, INFO );
         chkxer('SSYEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssyev('N', '/', 0, A, 1, X, W, 1, INFO );
         chkxer('SSYEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssyev('N', 'U', -1, A, 1, X, W, 1, INFO );
         chkxer('SSYEV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssyev('N', 'U', 2, A, 1, X, W, 3, INFO );
         chkxer('SSYEV ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssyev('N', 'U', 1, A, 1, X, W, 1, INFO );
         chkxer('SSYEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 5;

         // SSYEV_2STAGE

         SRNAMT = 'SSYEV_2STAGE ';
         INFOT = 1;
         ssyev_2stage('/', 'U', 0, A, 1, X, W, 1, INFO );
         chkxer('SSYEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         ssyev_2stage('V', 'U', 0, A, 1, X, W, 1, INFO );
         chkxer('SSYEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssyev_2stage('N', '/', 0, A, 1, X, W, 1, INFO );
         chkxer('SSYEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssyev_2stage('N', 'U', -1, A, 1, X, W, 1, INFO );
         chkxer('SSYEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssyev_2stage('N', 'U', 2, A, 1, X, W, 3, INFO );
         chkxer('SSYEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssyev_2stage('N', 'U', 1, A, 1, X, W, 1, INFO );
         chkxer('SSYEV_2STAGE ', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // SSYEVX

         SRNAMT = 'SSYEVX';
         INFOT = 1;
         ssyevx('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('SSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssyevx('N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('SSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssyevx('N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         INFOT = 4;
         ssyevx('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('SSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssyevx('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 16, IW, I3, INFO );
         chkxer('SSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssyevx('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('SSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ssyevx('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('SSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ssyevx('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('SSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssyevx('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 16, IW, I3, INFO );
         chkxer('SSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssyevx('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('SSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         ssyevx('V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 16, IW, I3, INFO );
         chkxer('SSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         ssyevx('V', 'A', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSYEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // SSYEVX_2STAGE

         SRNAMT = 'SSYEVX_2STAGE';
         INFOT = 1;
         ssyevx_2stage('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('SSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         ssyevx_2stage('V', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('SSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssyevx_2stage('N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('SSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssyevx_2stage('N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         INFOT = 4;
         ssyevx_2stage('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('SSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssyevx_2stage('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 16, IW, I3, INFO );
         chkxer('SSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssyevx_2stage('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('SSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ssyevx_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('SSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ssyevx_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('SSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssyevx_2stage('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 16, IW, I3, INFO );
         chkxer('SSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssyevx_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('SSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         ssyevx_2stage('N', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 0, W, 16, IW, I3, INFO );
         chkxer('SSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         ssyevx_2stage('N', 'A', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 13;

         // SSPEVD

         SRNAMT = 'SSPEVD';
         INFOT = 1;
         sspevd('/', 'U', 0, A, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sspevd('N', '/', 0, A, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sspevd('N', 'U', -1, A, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sspevd('V', 'U', 2, A, X, Z, 1, W, 23, IW, 12, INFO );
         chkxer('SSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sspevd('N', 'U', 1, A, X, Z, 1, W, 0, IW, 1, INFO );
         chkxer('SSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sspevd('N', 'U', 2, A, X, Z, 1, W, 3, IW, 1, INFO );
         chkxer('SSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sspevd('V', 'U', 2, A, X, Z, 2, W, 16, IW, 12, INFO );
         chkxer('SSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sspevd('N', 'U', 1, A, X, Z, 1, W, 1, IW, 0, INFO );
         chkxer('SSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sspevd('N', 'U', 2, A, X, Z, 1, W, 4, IW, 0, INFO );
         chkxer('SSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         sspevd('V', 'U', 2, A, X, Z, 2, W, 23, IW, 11, INFO );
         chkxer('SSPEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // SSPEV

         SRNAMT = 'SSPEV ';
         INFOT = 1;
         sspev('/', 'U', 0, A, W, Z, 1, X, INFO );
         chkxer('SSPEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sspev('N', '/', 0, A, W, Z, 1, X, INFO );
         chkxer('SSPEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sspev('N', 'U', -1, A, W, Z, 1, X, INFO );
         chkxer('SSPEV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sspev('V', 'U', 2, A, W, Z, 1, X, INFO );
         chkxer('SSPEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // SSPEVX

         SRNAMT = 'SSPEVX';
         INFOT = 1;
         sspevx('/', 'A', 'U', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sspevx('N', '/', 'U', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         sspevx('N', 'A', '/', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         INFOT = 4;
         sspevx('N', 'A', 'U', -1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         sspevx('N', 'V', 'U', 1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sspevx('N', 'I', 'U', 1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         sspevx('N', 'I', 'U', 1, A, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sspevx('N', 'I', 'U', 2, A, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         sspevx('N', 'I', 'U', 1, A, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         sspevx('V', 'A', 'U', 2, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSPEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

      // Test error exits for the SB path.

      } else if ( LSAMEN( 2, C2, 'SB' ) ) {

         // SSBTRD

         SRNAMT = 'SSBTRD';
         INFOT = 1;
         ssbtrd('/', 'U', 0, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('SSBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssbtrd('N', '/', 0, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('SSBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssbtrd('N', 'U', -1, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('SSBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssbtrd('N', 'U', 0, -1, A, 1, D, E, Z, 1, W, INFO );
         chkxer('SSBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssbtrd('N', 'U', 1, 1, A, 1, D, E, Z, 1, W, INFO );
         chkxer('SSBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssbtrd('V', 'U', 2, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('SSBTRD', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // SSYTRD_SB2ST

         SRNAMT = 'SSYTRD_SB2ST';
         INFOT = 1;
         ssytrd_sb2st('/', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrd_sb2st('N', '/', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrd_sb2st('N', 'H', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssytrd_sb2st('N', 'N', '/', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytrd_sb2st('N', 'N', 'U', -1, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssytrd_sb2st('N', 'N', 'U', 0, -1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ssytrd_sb2st('N', 'N', 'U', 0, 1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ssytrd_sb2st('N', 'N', 'U', 0, 0, A, 1, D, E,  C, 0, W, 1, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         ssytrd_sb2st('N', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 0, INFO );
         chkxer('SSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // SSBEVD

         SRNAMT = 'SSBEVD';
         INFOT = 1;
         ssbevd('/', 'U', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssbevd('N', '/', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssbevd('N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssbevd('N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssbevd('N', 'U', 2, 1, A, 1, X, Z, 1, W, 4, IW, 1, INFO );
         chkxer('SSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ssbevd('V', 'U', 2, 1, A, 2, X, Z, 1, W, 25, IW, 12, INFO );
         chkxer('SSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ssbevd('N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, IW, 1, INFO );
         chkxer('SSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ssbevd('N', 'U', 2, 0, A, 1, X, Z, 1, W, 3, IW, 1, INFO );
         chkxer('SSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ssbevd('V', 'U', 2, 0, A, 1, X, Z, 2, W, 18, IW, 12, INFO );
         chkxer('SSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         ssbevd('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, IW, 0, INFO );
         chkxer('SSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         ssbevd('V', 'U', 2, 0, A, 1, X, Z, 2, W, 25, IW, 11, INFO );
         chkxer('SSBEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // SSBEVD_2STAGE

         SRNAMT = 'SSBEVD_2STAGE';
         INFOT = 1;
         ssbevd_2stage('/', 'U', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         ssbevd_2stage('V', 'U', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssbevd_2stage('N', '/', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssbevd_2stage('N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssbevd_2stage('N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('SSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssbevd_2stage('N', 'U', 2, 1, A, 1, X, Z, 1, W, 4, IW, 1, INFO );
         chkxer('SSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 9
          // CALL SSBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 1, W,
      // $                                      25, IW, 12, INFO )
          // CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 11;
         ssbevd_2stage('N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, IW, 1, INFO );
         chkxer('SSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ssbevd_2stage('N', 'U', 2, 0, A, 1, X, Z, 1, W, 3, IW, 1, INFO );
         chkxer('SSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 11
          // CALL SSBEVD_2STAGE( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
      // $                                      18, IW, 12, INFO )
          // CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 13;
         ssbevd_2stage('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, IW, 0, INFO );
         chkxer('SSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 13
          // CALL SSBEVD_2STAGE( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
      // $                                      25, IW, 11, INFO )
          // CALL CHKXER( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
          // NT = NT + 12
         NT = NT + 9;

         // SSBEV

         SRNAMT = 'SSBEV ';
         INFOT = 1;
         ssbev('/', 'U', 0, 0, A, 1, X, Z, 1, W, INFO );
         chkxer('SSBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssbev('N', '/', 0, 0, A, 1, X, Z, 1, W, INFO );
         chkxer('SSBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssbev('N', 'U', -1, 0, A, 1, X, Z, 1, W, INFO );
         chkxer('SSBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssbev('N', 'U', 0, -1, A, 1, X, Z, 1, W, INFO );
         chkxer('SSBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssbev('N', 'U', 2, 1, A, 1, X, Z, 1, W, INFO );
         chkxer('SSBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ssbev('V', 'U', 2, 0, A, 1, X, Z, 1, W, INFO );
         chkxer('SSBEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // SSBEV_2STAGE

         SRNAMT = 'SSBEV_2STAGE ';
         INFOT = 1;
         ssbev_2stage('/', 'U', 0, 0, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('SSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         ssbev_2stage('V', 'U', 0, 0, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('SSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssbev_2stage('N', '/', 0, 0, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('SSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssbev_2stage('N', 'U', -1, 0, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('SSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssbev_2stage('N', 'U', 0, -1, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('SSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssbev_2stage('N', 'U', 2, 1, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('SSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ssbev_2stage('N', 'U', 2, 0, A, 1, X, Z, 0, W, 0, INFO );
         chkxer('SSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ssbev_2stage('N', 'U', 0, 0, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('SSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

         // SSBEVX

         SRNAMT = 'SSBEVX';
         INFOT = 1;
         ssbevx('/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssbevx('N', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssbevx('N', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssbevx('N', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssbevx('N', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ssbevx('N', 'A', 'U', 2, 1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ssbevx('V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         ssbevx('N', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         ssbevx('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         ssbevx('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         ssbevx('N', 'I', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         ssbevx('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         ssbevx('V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('SSBEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 13;

         // SSBEVX_2STAGE

         SRNAMT = 'SSBEVX_2STAGE';
         INFOT = 1;
         ssbevx_2stage('/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         ssbevx_2stage('V', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssbevx_2stage('N', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssbevx_2stage('N', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssbevx_2stage('N', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssbevx_2stage('N', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ssbevx_2stage('N', 'A', 'U', 2, 1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 9
          // CALL SSBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0,
      // $          0.0, 0, 0, 0.0, M, X, Z, 2, W, 0, IW, I3, INFO )
          // CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 11;
         ssbevx_2stage('N', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         ssbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         ssbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         ssbevx_2stage('N', 'I', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         ssbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 18
          // CALL SSBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0,
      // $          0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO )
          // CALL CHKXER( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 20;
         ssbevx_2stage('N', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('SSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
          // NT = NT + 15
         NT = NT + 13;
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
