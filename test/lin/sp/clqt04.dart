      void clqt04(M,N,NB,RESULT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     M, N, NB;
      // .. Return values ..
      double RESULT(6);

// =====================================================================

      // ..
      // .. Local allocatable arrays
      Complex, ALLOCATABLE :: AF(:,:), Q(:,:), L(:,:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:);
      double, ALLOCATABLE :: RWORK(:);

      // .. Parameters ..
      double       ZERO;
      Complex    ONE, CZERO;
      const    ZERO = 0.0;
      const    ONE = (1.0,0.0), CZERO=(0.0,0.0) ;
      int     INFO, J, K, LL, LWORK, LDT;
      double    ANORM, EPS, RESID, CNORM, DNORM;
      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      //- REAL     SLAMCH;
      //- REAL     CLANGE, CLANSY;
      //- bool     lsame;
      // EXTERNAL SLAMCH, CLANGE, CLANSY, lsame
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Data statements ..
      const ISEED = [ 1988, 1989, 1990, 1991 ];

      EPS = SLAMCH( 'Epsilon' );
      K = min(M,N);
      LL = max(M,N);
      LWORK = max(2,LL)*max(2,LL)*NB;

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(N,N), L(LL,N), RWORK(LL), WORK(LWORK), T(NB,N), C(M,N), CF(M,N), D(N,M), DF(N,M) );

      // Put random numbers into A and copy to AF

      LDT=NB;
      for (J = 1; J <= N; J++) {
         clarnv(2, ISEED, M, A( 1, J ) );
      }
      clacpy('Full', M, N, A, M, AF, M );

      // Factor the matrix A in the array AF.

      cgelqt(M, N, NB, AF, M, T, LDT, WORK, INFO );

      // Generate the n-by-n matrix Q

      claset('Full', N, N, CZERO, ONE, Q, N );
      cgemlqt('R', 'N', N, N, K, NB, AF, M, T, LDT, Q, N, WORK, INFO );

      // Copy L

      claset('Full', LL, N, CZERO, CZERO, L, LL );
      clacpy('Lower', M, N, AF, M, L, LL );

      // Compute |L - A*Q'| / |A| and store in RESULT(1)

      cgemm('N', 'C', M, N, N, -ONE, A, M, Q, N, ONE, L, LL );
      ANORM = CLANGE( '1', M, N, A, M, RWORK );
      RESID = CLANGE( '1', M, N, L, LL, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = RESID / (EPS*max(1,M)*ANORM);
      } else {
         RESULT[1] = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      claset('Full', N, N, CZERO, ONE, L, LL );
      cherk('U', 'C', N, N, REAL(-ONE), Q, N, double(ONE), L, LL);
      RESID = CLANSY( '1', 'Upper', N, L, LL, RWORK );
      RESULT[2] = RESID / (EPS*max(1,N));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= M; J++) {
         clarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = CLANGE( '1', N, M, D, N, RWORK);
      clacpy('Full', N, M, D, N, DF, N );

      // Apply Q to C as Q*C

      cgemlqt('L', 'N', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

      // Compute |Q*D - Q*D| / |D|

      cgemm('N', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N );
      RESID = CLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT[3] = RESID / (EPS*max(1,M)*DNORM);
      } else {
         RESULT[3] = ZERO;
      }

      // Copy D into DF again

      clacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as QT*D

      cgemlqt('L', 'C', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

      // Compute |QT*D - QT*D| / |D|

      cgemm('C', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N );
      RESID = CLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT[4] = RESID / (EPS*max(1,M)*DNORM);
      } else {
         RESULT[4] = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= N; J++) {
         clarnv(2, ISEED, M, C( 1, J ) );
      }
      CNORM = CLANGE( '1', M, N, C, M, RWORK);
      clacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as C*Q

      cgemlqt('R', 'N', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

      // Compute |C*Q - C*Q| / |C|

      cgemm('N', 'N', M, N, N, -ONE, C, M, Q, N, ONE, CF, M );
      RESID = CLANGE( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[5] = RESID / (EPS*max(1,M)*DNORM);
      } else {
         RESULT[5] = ZERO;
      }

      // Copy C into CF again

      clacpy('Full', M, N, C, M, CF, M );

      // Apply Q to D as D*QT

      cgemlqt('R', 'C', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

      // Compute |C*QT - C*QT| / |C|

      cgemm('N', 'C', M, N, N, -ONE, C, M, Q, N, ONE, CF, M );
      RESID = CLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[6] = RESID / (EPS*max(1,M)*DNORM);
      } else {
         RESULT[6] = ZERO;
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, L, RWORK, WORK, T, C, D, CF, DF);

      }
