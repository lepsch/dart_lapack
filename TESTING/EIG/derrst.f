      SUBROUTINE DERRST( PATH, NUNIT );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

*  =====================================================================

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
      double             A( NMAX, NMAX ), C( NMAX, NMAX ), D( NMAX ), E( NMAX ), Q( NMAX, NMAX ), R( NMAX ), TAU( NMAX ), W( LW ), X( NMAX ), Z( NMAX, NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHKXER, DOPGTR, DOPMTR, DORGTR, DORMTR, DPTEQR, DSBEV, DSBEVD, DSBEVX, DSBTRD, DSPEV, DSPEVD, DSPEVX, DSPTRD, DSTEBZ, DSTEDC, DSTEIN, DSTEQR, DSTERF, DSTEV, DSTEVD, DSTEVR, DSTEVX, DSYEV, DSYEVD, DSYEVR, DSYEVX, DSYTRD, DSYTD2, DSYEVD_2STAGE, DSYEVR_2STAGE, DSYEVX_2STAGE, DSYEV_2STAGE, DSBEV_2STAGE, DSBEVD_2STAGE, DSBEVX_2STAGE, DSYTRD_2STAGE, DSYTRD_SY2SB, DSYTRD_SB2ST
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
            A( I, J ) = 1.0 / DBLE( I+J );
         } // 10
      } // 20
      for (J = 1; J <= NMAX; J++) { // 30
         D( J ) = DBLE( J );
         E( J ) = 0.0;
         I1( J ) = J;
         I2( J ) = J;
         TAU( J ) = 1.0;
      } // 30
      OK = true;
      NT = 0;

      // Test error exits for the ST path.

      if ( LSAMEN( 2, C2, 'ST' ) ) {

         // DSYTRD

         SRNAMT = 'DSYTRD';
         INFOT = 1;
         dsytrd('/', 0, A, 1, D, E, TAU, W, 1, INFO );
         chkxer('DSYTRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrd('U', -1, A, 1, D, E, TAU, W, 1, INFO );
         chkxer('DSYTRD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytrd('U', 2, A, 1, D, E, TAU, W, 1, INFO );
         chkxer('DSYTRD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsytrd('U', 0, A, 1, D, E, TAU, W, 0, INFO );
         chkxer('DSYTRD', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // DSYTD2

         SRNAMT = 'DSYTD2';
         INFOT = 1;
         dsytd2('/', 0, A, 1, D, E, TAU, INFO );
         chkxer('DSYTD2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytd2('U', -1, A, 1, D, E, TAU, INFO );
         chkxer('DSYTD2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytd2('U', 2, A, 1, D, E, TAU, INFO );
         chkxer('DSYTD2', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // DSYTRD_2STAGE

         SRNAMT = 'DSYTRD_2STAGE';
         INFOT = 1;
         dsytrd_2stage('/', 'U', 0, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         dsytrd_2stage('H', 'U', 0, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrd_2stage('N', '/', 0, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsytrd_2stage('N', 'U', -1, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsytrd_2stage('N', 'U', 2, A, 1, D, E, TAU,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsytrd_2stage('N', 'U', 0, A, 1, D, E, TAU,  C, 0, W, 1, INFO );
         chkxer('DSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dsytrd_2stage('N', 'U', 0, A, 1, D, E, TAU,  C, 1, W, 0, INFO );
         chkxer('DSYTRD_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // DSYTRD_SY2SB

         SRNAMT = 'DSYTRD_SY2SB';
         INFOT = 1;
         dsytrd_sy2sb('/', 0, 0, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('DSYTRD_SY2SB', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrd_sy2sb('U', -1, 0, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('DSYTRD_SY2SB', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsytrd_sy2sb('U', 0, -1, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('DSYTRD_SY2SB', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsytrd_sy2sb('U', 2, 0, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('DSYTRD_SY2SB', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dsytrd_sy2sb('U', 0, 2, A, 1, C, 1, TAU, W, 1, INFO );
         chkxer('DSYTRD_SY2SB', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsytrd_sy2sb('U', 0, 0, A, 1, C, 1, TAU, W, 0, INFO );
         chkxer('DSYTRD_SY2SB', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // DSYTRD_SB2ST

         SRNAMT = 'DSYTRD_SB2ST';
         INFOT = 1;
         dsytrd_sb2st('/', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrd_sb2st('Y', '/', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrd_sb2st('Y', 'H', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsytrd_sb2st('Y', 'N', '/', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytrd_sb2st('Y', 'N', 'U', -1, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsytrd_sb2st('Y', 'N', 'U', 0, -1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dsytrd_sb2st('Y', 'N', 'U', 0, 1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dsytrd_sb2st('Y', 'N', 'U', 0, 0, A, 1, D, E,  C, 0, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dsytrd_sb2st('Y', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 0, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // DORGTR

         SRNAMT = 'DORGTR';
         INFOT = 1;
         dorgtr('/', 0, A, 1, TAU, W, 1, INFO );
         chkxer('DORGTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dorgtr('U', -1, A, 1, TAU, W, 1, INFO );
         chkxer('DORGTR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dorgtr('U', 2, A, 1, TAU, W, 1, INFO );
         chkxer('DORGTR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dorgtr('U', 3, A, 3, TAU, W, 1, INFO );
         chkxer('DORGTR', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // DORMTR

         SRNAMT = 'DORMTR';
         INFOT = 1;
         dormtr('/', 'U', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dormtr('L', '/', 'N', 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dormtr('L', 'U', '/', 0, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dormtr('L', 'U', 'N', -1, 0, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dormtr('L', 'U', 'N', 0, -1, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dormtr('L', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('DORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dormtr('R', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dormtr('L', 'U', 'N', 2, 0, A, 2, TAU, C, 1, W, 1, INFO );
         chkxer('DORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dormtr('L', 'U', 'N', 0, 2, A, 1, TAU, C, 1, W, 1, INFO );
         chkxer('DORMTR', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dormtr('R', 'U', 'N', 2, 0, A, 1, TAU, C, 2, W, 1, INFO );
         chkxer('DORMTR', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // DSPTRD

         SRNAMT = 'DSPTRD';
         INFOT = 1;
         dsptrd('/', 0, A, D, E, TAU, INFO );
         chkxer('DSPTRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsptrd('U', -1, A, D, E, TAU, INFO );
         chkxer('DSPTRD', INFOT, NOUT, LERR, OK );
         NT = NT + 2;

         // DOPGTR

         SRNAMT = 'DOPGTR';
         INFOT = 1;
         dopgtr('/', 0, A, TAU, Z, 1, W, INFO );
         chkxer('DOPGTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dopgtr('U', -1, A, TAU, Z, 1, W, INFO );
         chkxer('DOPGTR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dopgtr('U', 2, A, TAU, Z, 1, W, INFO );
         chkxer('DOPGTR', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // DOPMTR

         SRNAMT = 'DOPMTR';
         INFOT = 1;
         dopmtr('/', 'U', 'N', 0, 0, A, TAU, C, 1, W, INFO );
         chkxer('DOPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dopmtr('L', '/', 'N', 0, 0, A, TAU, C, 1, W, INFO );
         chkxer('DOPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dopmtr('L', 'U', '/', 0, 0, A, TAU, C, 1, W, INFO );
         chkxer('DOPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dopmtr('L', 'U', 'N', -1, 0, A, TAU, C, 1, W, INFO );
         chkxer('DOPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dopmtr('L', 'U', 'N', 0, -1, A, TAU, C, 1, W, INFO );
         chkxer('DOPMTR', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dopmtr('L', 'U', 'N', 2, 0, A, TAU, C, 1, W, INFO );
         chkxer('DOPMTR', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // DPTEQR

         SRNAMT = 'DPTEQR';
         INFOT = 1;
         dpteqr('/', 0, D, E, Z, 1, W, INFO );
         chkxer('DPTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dpteqr('N', -1, D, E, Z, 1, W, INFO );
         chkxer('DPTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dpteqr('V', 2, D, E, Z, 1, W, INFO );
         chkxer('DPTEQR', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // DSTEBZ

         SRNAMT = 'DSTEBZ';
         INFOT = 1;
         dstebz('/', 'E', 0, 0.0, 1.0, 1, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('DSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dstebz('A', '/', 0, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('DSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dstebz('A', 'E', -1, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('DSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dstebz('V', 'E', 0, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('DSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dstebz('I', 'E', 0, 0.0, 0.0, 0, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('DSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dstebz('I', 'E', 1, 0.0, 0.0, 2, 1, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('DSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dstebz('I', 'E', 1, 0.0, 0.0, 1, 0, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('DSTEBZ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dstebz('I', 'E', 1, 0.0, 0.0, 1, 2, 0.0, D, E, M, NSPLIT, X, I1, I2, W, IW, INFO );
         chkxer('DSTEBZ', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

         // DSTEIN

         SRNAMT = 'DSTEIN';
         INFOT = 1;
         dstein(-1, D, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEIN', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dstein(0, D, E, -1, X, I1, I2, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEIN', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dstein(0, D, E, 1, X, I1, I2, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEIN', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dstein(2, D, E, 0, X, I1, I2, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEIN', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // DSTEQR

         SRNAMT = 'DSTEQR';
         INFOT = 1;
         dsteqr('/', 0, D, E, Z, 1, W, INFO );
         chkxer('DSTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsteqr('N', -1, D, E, Z, 1, W, INFO );
         chkxer('DSTEQR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsteqr('V', 2, D, E, Z, 1, W, INFO );
         chkxer('DSTEQR', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // DSTERF

         SRNAMT = 'DSTERF';
         INFOT = 1;
         dsterf(-1, D, E, INFO );
         chkxer('DSTERF', INFOT, NOUT, LERR, OK );
         NT = NT + 1;

         // DSTEDC

         SRNAMT = 'DSTEDC';
         INFOT = 1;
         dstedc('/', 0, D, E, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dstedc('N', -1, D, E, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dstedc('V', 2, D, E, Z, 1, W, 23, IW, 28, INFO );
         chkxer('DSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dstedc('N', 1, D, E, Z, 1, W, 0, IW, 1, INFO );
         chkxer('DSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dstedc('I', 2, D, E, Z, 2, W, 0, IW, 12, INFO );
         chkxer('DSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dstedc('V', 2, D, E, Z, 2, W, 0, IW, 28, INFO );
         chkxer('DSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dstedc('N', 1, D, E, Z, 1, W, 1, IW, 0, INFO );
         chkxer('DSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dstedc('I', 2, D, E, Z, 2, W, 19, IW, 0, INFO );
         chkxer('DSTEDC', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dstedc('V', 2, D, E, Z, 2, W, 23, IW, 0, INFO );
         chkxer('DSTEDC', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // DSTEVD

         SRNAMT = 'DSTEVD';
         INFOT = 1;
         dstevd('/', 0, D, E, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSTEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dstevd('N', -1, D, E, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSTEVD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dstevd('V', 2, D, E, Z, 1, W, 19, IW, 12, INFO );
         chkxer('DSTEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dstevd('N', 1, D, E, Z, 1, W, 0, IW, 1, INFO );
         chkxer('DSTEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dstevd('V', 2, D, E, Z, 2, W, 12, IW, 12, INFO );
         chkxer('DSTEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dstevd('N', 0, D, E, Z, 1, W, 1, IW, 0, INFO );
         chkxer('DSTEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dstevd('V', 2, D, E, Z, 2, W, 19, IW, 11, INFO );
         chkxer('DSTEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 7;

         // DSTEV

         SRNAMT = 'DSTEV ';
         INFOT = 1;
         dstev('/', 0, D, E, Z, 1, W, INFO );
         chkxer('DSTEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dstev('N', -1, D, E, Z, 1, W, INFO );
         chkxer('DSTEV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dstev('V', 2, D, E, Z, 1, W, INFO );
         chkxer('DSTEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 3;

         // DSTEVX

         SRNAMT = 'DSTEVX';
         INFOT = 1;
         dstevx('/', 'A', 0, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dstevx('N', '/', 0, D, E, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dstevx('N', 'A', -1, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dstevx('N', 'V', 1, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dstevx('N', 'I', 1, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dstevx('N', 'I', 1, D, E, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dstevx('N', 'I', 2, D, E, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dstevx('N', 'I', 1, D, E, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dstevx('V', 'A', 2, D, E, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSTEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // DSTEVR

         N = 1;
         SRNAMT = 'DSTEVR';
         INFOT = 1;
         dstevr('/', 'A', 0, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dstevr('V', '/', 0, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dstevr('V', 'A', -1, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dstevr('V', 'V', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dstevr('V', 'I', 1, D, E, 0.0, 0.0, 0, 1, 0.0, M, W, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         N = 2;
         dstevr('V', 'I', 2, D, E, 0.0, 0.0, 2, 1, 0.0, M, W, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         N = 1;
         dstevr('V', 'I', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, W, Z, 0, IW, X, 20*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         dstevr('V', 'I', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, W, Z, 1, IW, X, 20*N-1, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSTEVR', INFOT, NOUT, LERR, OK );
         INFOT = 19;
         dstevr('V', 'I', 1, D, E, 0.0, 0.0, 1, 1, 0.0, M, W, Z, 1, IW, X, 20*N, IW( 2*N+1 ), 10*N-1, INFO );
         chkxer('DSTEVR', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // DSYEVD

         SRNAMT = 'DSYEVD';
         INFOT = 1;
         dsyevd('/', 'U', 0, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('DSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsyevd('N', '/', 0, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('DSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsyevd('N', 'U', -1, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('DSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsyevd('N', 'U', 2, A, 1, X, W, 3, IW, 1, INFO );
         chkxer('DSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsyevd('N', 'U', 1, A, 1, X, W, 0, IW, 1, INFO );
         chkxer('DSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsyevd('N', 'U', 2, A, 2, X, W, 4, IW, 1, INFO );
         chkxer('DSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsyevd('V', 'U', 2, A, 2, X, W, 20, IW, 12, INFO );
         chkxer('DSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsyevd('N', 'U', 1, A, 1, X, W, 1, IW, 0, INFO );
         chkxer('DSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsyevd('N', 'U', 2, A, 2, X, W, 5, IW, 0, INFO );
         chkxer('DSYEVD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsyevd('V', 'U', 2, A, 2, X, W, 27, IW, 11, INFO );
         chkxer('DSYEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // DSYEVD_2STAGE

         SRNAMT = 'DSYEVD_2STAGE';
         INFOT = 1;
         dsyevd_2stage('/', 'U', 0, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('DSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         dsyevd_2stage('V', 'U', 0, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('DSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsyevd_2stage('N', '/', 0, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('DSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsyevd_2stage('N', 'U', -1, A, 1, X, W, 1, IW, 1, INFO );
         chkxer('DSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsyevd_2stage('N', 'U', 2, A, 1, X, W, 3, IW, 1, INFO );
         chkxer('DSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsyevd_2stage('N', 'U', 1, A, 1, X, W, 0, IW, 1, INFO );
         chkxer('DSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsyevd_2stage('N', 'U', 2, A, 2, X, W, 4, IW, 1, INFO );
         chkxer('DSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 8
          // CALL DSYEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 20, IW, 12, INFO )
          // CALL CHKXER( 'DSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 10;
         dsyevd_2stage('N', 'U', 1, A, 1, X, W, 1, IW, 0, INFO );
         chkxer('DSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsyevd_2stage('N', 'U', 2, A, 2, X, W, 25, IW, 0, INFO );
         chkxer('DSYEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 10
          // CALL DSYEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 27, IW, 11, INFO )
          // CALL CHKXER( 'DSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         NT = NT + 9;

         // DSYEVR

         SRNAMT = 'DSYEVR';
         N = 1;
         INFOT = 1;
         dsyevr('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsyevr('V', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsyevr('V', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsyevr('V', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsyevr('V', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsyevr('V', 'V', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsyevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 10;

         dsyevr('V', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dsyevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 0, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         dsyevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 0, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         dsyevr('V', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 0, INFO );
         chkxer('DSYEVR', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // DSYEVR_2STAGE

         SRNAMT = 'DSYEVR_2STAGE';
         N = 1;
         INFOT = 1;
         dsyevr_2stage('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         dsyevr_2stage('V', 'A', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsyevr_2stage('N', '/', 'U', 0, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsyevr_2stage('N', 'A', '/', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsyevr_2stage('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsyevr_2stage('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsyevr_2stage('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsyevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsyevr_2stage('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dsyevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 0, IW, Q, 26*N, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         dsyevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 0, IW( 2*N+1 ), 10*N, INFO );
         chkxer('DSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 20;
         dsyevr_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 1, 0.0, M, R, Z, 1, IW, Q, 26*N, IW( 2*N+1 ), 0, INFO );
         chkxer('DSYEVR_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // DSYEV

         SRNAMT = 'DSYEV ';
         INFOT = 1;
         dsyev('/', 'U', 0, A, 1, X, W, 1, INFO );
         chkxer('DSYEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsyev('N', '/', 0, A, 1, X, W, 1, INFO );
         chkxer('DSYEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsyev('N', 'U', -1, A, 1, X, W, 1, INFO );
         chkxer('DSYEV ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsyev('N', 'U', 2, A, 1, X, W, 3, INFO );
         chkxer('DSYEV ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsyev('N', 'U', 1, A, 1, X, W, 1, INFO );
         chkxer('DSYEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 5;

         // DSYEV_2STAGE

         SRNAMT = 'DSYEV_2STAGE ';
         INFOT = 1;
         dsyev_2stage('/', 'U', 0, A, 1, X, W, 1, INFO );
         chkxer('DSYEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         dsyev_2stage('V', 'U', 0, A, 1, X, W, 1, INFO );
         chkxer('DSYEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsyev_2stage('N', '/', 0, A, 1, X, W, 1, INFO );
         chkxer('DSYEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsyev_2stage('N', 'U', -1, A, 1, X, W, 1, INFO );
         chkxer('DSYEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsyev_2stage('N', 'U', 2, A, 1, X, W, 3, INFO );
         chkxer('DSYEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsyev_2stage('N', 'U', 1, A, 1, X, W, 1, INFO );
         chkxer('DSYEV_2STAGE ', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // DSYEVX

         SRNAMT = 'DSYEVX';
         INFOT = 1;
         dsyevx('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('DSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsyevx('N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('DSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsyevx('N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         INFOT = 4;
         dsyevx('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('DSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsyevx('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 16, IW, I3, INFO );
         chkxer('DSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsyevx('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('DSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsyevx('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('DSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsyevx('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('DSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsyevx('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 16, IW, I3, INFO );
         chkxer('DSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsyevx('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('DSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dsyevx('V', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 16, IW, I3, INFO );
         chkxer('DSYEVX', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         dsyevx('V', 'A', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSYEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 12;

         // DSYEVX_2STAGE

         SRNAMT = 'DSYEVX_2STAGE';
         INFOT = 1;
         dsyevx_2stage('/', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('DSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         dsyevx_2stage('V', 'A', 'U', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('DSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsyevx_2stage('N', '/', 'U', 0, A, 1, 0.0, 1.0, 1, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('DSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsyevx_2stage('N', 'A', '/', 0, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         INFOT = 4;
         dsyevx_2stage('N', 'A', 'U', -1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 1, IW, I3, INFO );
         chkxer('DSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsyevx_2stage('N', 'A', 'U', 2, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 16, IW, I3, INFO );
         chkxer('DSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsyevx_2stage('N', 'V', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('DSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsyevx_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('DSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsyevx_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('DSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsyevx_2stage('N', 'I', 'U', 2, A, 2, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 16, IW, I3, INFO );
         chkxer('DSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsyevx_2stage('N', 'I', 'U', 1, A, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, 8, IW, I3, INFO );
         chkxer('DSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 15;
         dsyevx_2stage('N', 'A', 'U', 2, A, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 0, W, 16, IW, I3, INFO );
         chkxer('DSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 17;
         dsyevx_2stage('N', 'A', 'U', 1, A, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSYEVX_2STAGE', INFOT, NOUT, LERR, OK );
         NT = NT + 13;

         // DSPEVD

         SRNAMT = 'DSPEVD';
         INFOT = 1;
         dspevd('/', 'U', 0, A, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dspevd('N', '/', 0, A, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dspevd('N', 'U', -1, A, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dspevd('V', 'U', 2, A, X, Z, 1, W, 23, IW, 12, INFO );
         chkxer('DSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dspevd('N', 'U', 1, A, X, Z, 1, W, 0, IW, 1, INFO );
         chkxer('DSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dspevd('N', 'U', 2, A, X, Z, 1, W, 3, IW, 1, INFO );
         chkxer('DSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dspevd('V', 'U', 2, A, X, Z, 2, W, 16, IW, 12, INFO );
         chkxer('DSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dspevd('N', 'U', 1, A, X, Z, 1, W, 1, IW, 0, INFO );
         chkxer('DSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dspevd('N', 'U', 2, A, X, Z, 1, W, 4, IW, 0, INFO );
         chkxer('DSPEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dspevd('V', 'U', 2, A, X, Z, 2, W, 23, IW, 11, INFO );
         chkxer('DSPEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

         // DSPEV

         SRNAMT = 'DSPEV ';
         INFOT = 1;
         dspev('/', 'U', 0, A, W, Z, 1, X, INFO );
         chkxer('DSPEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dspev('N', '/', 0, A, W, Z, 1, X, INFO );
         chkxer('DSPEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dspev('N', 'U', -1, A, W, Z, 1, X, INFO );
         chkxer('DSPEV ', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dspev('V', 'U', 2, A, W, Z, 1, X, INFO );
         chkxer('DSPEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 4;

         // DSPEVX

         SRNAMT = 'DSPEVX';
         INFOT = 1;
         dspevx('/', 'A', 'U', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dspevx('N', '/', 'U', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dspevx('N', 'A', '/', 0, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         INFOT = 4;
         dspevx('N', 'A', 'U', -1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dspevx('N', 'V', 'U', 1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dspevx('N', 'I', 'U', 1, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dspevx('N', 'I', 'U', 1, A, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dspevx('N', 'I', 'U', 2, A, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dspevx('N', 'I', 'U', 1, A, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSPEVX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dspevx('V', 'A', 'U', 2, A, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSPEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 10;

      // Test error exits for the SB path.

      } else if ( LSAMEN( 2, C2, 'SB' ) ) {

         // DSBTRD

         SRNAMT = 'DSBTRD';
         INFOT = 1;
         dsbtrd('/', 'U', 0, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('DSBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsbtrd('N', '/', 0, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('DSBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsbtrd('N', 'U', -1, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('DSBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsbtrd('N', 'U', 0, -1, A, 1, D, E, Z, 1, W, INFO );
         chkxer('DSBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsbtrd('N', 'U', 1, 1, A, 1, D, E, Z, 1, W, INFO );
         chkxer('DSBTRD', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsbtrd('V', 'U', 2, 0, A, 1, D, E, Z, 1, W, INFO );
         chkxer('DSBTRD', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // DSYTRD_SB2ST

         SRNAMT = 'DSYTRD_SB2ST';
         INFOT = 1;
         dsytrd_sb2st('/', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrd_sb2st('N', '/', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrd_sb2st('N', 'H', 'U', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsytrd_sb2st('N', 'N', '/', 0, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytrd_sb2st('N', 'N', 'U', -1, 0, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsytrd_sb2st('N', 'N', 'U', 0, -1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dsytrd_sb2st('N', 'N', 'U', 0, 1, A, 1, D, E,  C, 1, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dsytrd_sb2st('N', 'N', 'U', 0, 0, A, 1, D, E,  C, 0, W, 1, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dsytrd_sb2st('N', 'N', 'U', 0, 0, A, 1, D, E,  C, 1, W, 0, INFO );
         chkxer('DSYTRD_SB2ST', INFOT, NOUT, LERR, OK );
         NT = NT + 9;

         // DSBEVD

         SRNAMT = 'DSBEVD';
         INFOT = 1;
         dsbevd('/', 'U', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsbevd('N', '/', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsbevd('N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsbevd('N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsbevd('N', 'U', 2, 1, A, 1, X, Z, 1, W, 4, IW, 1, INFO );
         chkxer('DSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsbevd('V', 'U', 2, 1, A, 2, X, Z, 1, W, 25, IW, 12, INFO );
         chkxer('DSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dsbevd('N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, IW, 1, INFO );
         chkxer('DSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dsbevd('N', 'U', 2, 0, A, 1, X, Z, 1, W, 3, IW, 1, INFO );
         chkxer('DSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dsbevd('V', 'U', 2, 0, A, 1, X, Z, 2, W, 18, IW, 12, INFO );
         chkxer('DSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dsbevd('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, IW, 0, INFO );
         chkxer('DSBEVD', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dsbevd('V', 'U', 2, 0, A, 1, X, Z, 2, W, 25, IW, 11, INFO );
         chkxer('DSBEVD', INFOT, NOUT, LERR, OK );
         NT = NT + 11;

         // DSBEVD_2STAGE

         SRNAMT = 'DSBEVD_2STAGE';
         INFOT = 1;
         dsbevd_2stage('/', 'U', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         dsbevd_2stage('V', 'U', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsbevd_2stage('N', '/', 0, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsbevd_2stage('N', 'U', -1, 0, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsbevd_2stage('N', 'U', 0, -1, A, 1, X, Z, 1, W, 1, IW, 1, INFO );
         chkxer('DSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsbevd_2stage('N', 'U', 2, 1, A, 1, X, Z, 1, W, 4, IW, 1, INFO );
         chkxer('DSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 9
          // CALL DSBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 1, W,
      // $                                      25, IW, 12, INFO )
          // CALL CHKXER( 'DSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 11;
         dsbevd_2stage('N', 'U', 1, 0, A, 1, X, Z, 1, W, 0, IW, 1, INFO );
         chkxer('DSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dsbevd_2stage('N', 'U', 2, 0, A, 1, X, Z, 1, W, 3, IW, 1, INFO );
         chkxer('DSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 11
          // CALL DSBEVD_2STAGE( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
      // $                                      18, IW, 12, INFO )
          // CALL CHKXER( 'DSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 13;
         dsbevd_2stage('N', 'U', 1, 0, A, 1, X, Z, 1, W, 1, IW, 0, INFO );
         chkxer('DSBEVD_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 13
          // CALL DSBEVD_2STAGE( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
      // $                                      25, IW, 11, INFO )
          // CALL CHKXER( 'DSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
          // NT = NT + 12
         NT = NT + 9;

         // DSBEV

         SRNAMT = 'DSBEV ';
         INFOT = 1;
         dsbev('/', 'U', 0, 0, A, 1, X, Z, 1, W, INFO );
         chkxer('DSBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsbev('N', '/', 0, 0, A, 1, X, Z, 1, W, INFO );
         chkxer('DSBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsbev('N', 'U', -1, 0, A, 1, X, Z, 1, W, INFO );
         chkxer('DSBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsbev('N', 'U', 0, -1, A, 1, X, Z, 1, W, INFO );
         chkxer('DSBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsbev('N', 'U', 2, 1, A, 1, X, Z, 1, W, INFO );
         chkxer('DSBEV ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsbev('V', 'U', 2, 0, A, 1, X, Z, 1, W, INFO );
         chkxer('DSBEV ', INFOT, NOUT, LERR, OK );
         NT = NT + 6;

         // DSBEV_2STAGE

         SRNAMT = 'DSBEV_2STAGE ';
         INFOT = 1;
         dsbev_2stage('/', 'U', 0, 0, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('DSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         dsbev_2stage('V', 'U', 0, 0, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('DSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsbev_2stage('N', '/', 0, 0, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('DSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsbev_2stage('N', 'U', -1, 0, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('DSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsbev_2stage('N', 'U', 0, -1, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('DSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsbev_2stage('N', 'U', 2, 1, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('DSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsbev_2stage('N', 'U', 2, 0, A, 1, X, Z, 0, W, 0, INFO );
         chkxer('DSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dsbev_2stage('N', 'U', 0, 0, A, 1, X, Z, 1, W, 0, INFO );
         chkxer('DSBEV_2STAGE ', INFOT, NOUT, LERR, OK );
         NT = NT + 8;

         // DSBEVX

         SRNAMT = 'DSBEVX';
         INFOT = 1;
         dsbevx('/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsbevx('N', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsbevx('N', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsbevx('N', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsbevx('N', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dsbevx('N', 'A', 'U', 2, 1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsbevx('V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 2, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         dsbevx('N', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dsbevx('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dsbevx('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dsbevx('N', 'I', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dsbevx('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         INFOT = 18;
         dsbevx('V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, IW, I3, INFO );
         chkxer('DSBEVX', INFOT, NOUT, LERR, OK );
         NT = NT + 13;

         // DSBEVX_2STAGE

         SRNAMT = 'DSBEVX_2STAGE';
         INFOT = 1;
         dsbevx_2stage('/', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 1;
         dsbevx_2stage('V', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsbevx_2stage('N', '/', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsbevx_2stage('N', 'A', '/', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsbevx_2stage('N', 'A', 'U', -1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsbevx_2stage('N', 'A', 'U', 0, -1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dsbevx_2stage('N', 'A', 'U', 2, 1, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 9
          // CALL DSBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0,
      // $          0.0, 0, 0, 0.0, M, X, Z, 2, W, 0, IW, I3, INFO )
          // CALL CHKXER( 'DSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 11;
         dsbevx_2stage('N', 'V', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dsbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dsbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dsbevx_2stage('N', 'I', 'U', 2, 0, A, 1, Q, 1, 0.0, 0.0, 2, 1, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 13;
         dsbevx_2stage('N', 'I', 'U', 1, 0, A, 1, Q, 1, 0.0, 0.0, 1, 2, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
          // INFOT = 18
          // CALL DSBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0,
      // $          0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO )
          // CALL CHKXER( 'DSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFOT = 20;
         dsbevx_2stage('N', 'A', 'U', 0, 0, A, 1, Q, 1, 0.0, 0.0, 0, 0, 0.0, M, X, Z, 1, W, 0, IW, I3, INFO );
         chkxer('DSBEVX_2STAGE', INFOT, NOUT, LERR, OK );
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

      // End of DERRST

      }
