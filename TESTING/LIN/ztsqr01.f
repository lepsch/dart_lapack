      void ztsqr01(TSSW, M, N, MB, NB, RESULT) {
      // IMPLICIT NONE

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String            TSSW;
      int               M, N, MB, NB;
      // .. Return values ..
      double            RESULT(6);

// =====================================================================

      // ..
      // .. Local allocatable arrays
      COMPLEX*16, ALLOCATABLE :: AF(:,:), Q(:,:), R(:,:), WORK( : ), T(:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:), LQ(:,:);
      double          , ALLOCATABLE :: RWORK(:);

      // .. Parameters ..
      double           ZERO;
      COMPLEX*16 ONE, CZERO;
      const    ZERO = 0.0, ONE = (1.0,0.0), CZERO=(0.0,0.0) ;
      // ..
      // .. Local Scalars ..
      bool    TESTZEROS, TS;
      int     INFO, J, K, L, LWORK, TSIZE, MNB;
      double             ANORM, EPS, RESID, CNORM, DNORM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      COMPLEX*16         TQUERY( 5 ), WORKQUERY( 1 );
      // ..
      // .. External Functions ..
      double           DLAMCH, ZLANGE, ZLANSY;
      bool     LSAME;
      int     ILAENV;
      // EXTERNAL DLAMCH, ZLANGE, ZLANSY, LSAME, ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // .. Scalars in Common ..
      String             srnamt;
      // ..
      // .. Common blocks ..
      // COMMON / srnamc / srnamt
      // ..
      // .. Data statements ..
      const ISEED = [ 1988, 1989, 1990, 1991 ];

      // TEST TALL SKINNY OR SHORT WIDE

      TS = LSAME(TSSW, 'TS');

      // TEST MATRICES WITH HALF OF MATRIX BEING ZEROS

      TESTZEROS = false;

      EPS = DLAMCH( 'Epsilon' );
      K = min(M,N);
      L = max(M,N,1);
      MNB = max( MB, NB);
      LWORK = max(3,L)*MNB;

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(L,L), R(M,L), RWORK(L), C(M,N), CF(M,N), D(N,M), DF(N,M), LQ(L,N) );

      // Put random numbers into A and copy to AF

      for (J = 1; J <= N; J++) {
         zlarnv(2, ISEED, M, A( 1, J ) );
      }
      if (TESTZEROS) {
         if (M >= 4) {
            for (J = 1; J <= N; J++) {
               zlarnv(2, ISEED, M/2, A( M/4, J ) );
            }
         }
      }
      zlacpy('Full', M, N, A, M, AF, M );

      if (TS) {

      // Factor the matrix A in the array AF.

      zgeqr(M, N, AF, M, TQUERY, -1, WORKQUERY, -1, INFO );
      TSIZE = INT( TQUERY( 1 ) );
      LWORK = INT( WORKQUERY( 1 ) );
      zgemqr('L', 'N', M, M, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      zgemqr('L', 'N', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      zgemqr('L', 'C', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      zgemqr('R', 'N', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      zgemqr('R', 'C', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      ALLOCATE ( T( TSIZE ) );
      ALLOCATE ( WORK( LWORK ) );
      srnamt = 'ZGEQR';
      zgeqr(M, N, AF, M, T, TSIZE, WORK, LWORK, INFO );

      // Generate the m-by-m matrix Q

      zlaset('Full', M, M, CZERO, ONE, Q, M );
      srnamt = 'ZGEMQR';
      zgemqr('L', 'N', M, M, K, AF, M, T, TSIZE, Q, M, WORK, LWORK, INFO );

      // Copy R

      zlaset('Full', M, N, CZERO, CZERO, R, M );
      zlacpy('Upper', M, N, AF, M, R, M );

      // Compute |R - Q'*A| / |A| and store in RESULT(1)

      zgemm('C', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M );
      ANORM = ZLANGE( '1', M, N, A, M, RWORK );
      RESID = ZLANGE( '1', M, N, R, M, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = RESID / (EPS*max(1,M)*ANORM);
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      zlaset('Full', M, M, CZERO, ONE, R, M );
      zherk('U', 'C', M, M, DREAL(-ONE), Q, M, DREAL(ONE), R, M );
      RESID = ZLANSY( '1', 'Upper', M, R, M, RWORK );
      RESULT( 2 ) = RESID / (EPS*max(1,M));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= N; J++) {
         zlarnv(2, ISEED, M, C( 1, J ) );
      }
      CNORM = ZLANGE( '1', M, N, C, M, RWORK);
      zlacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as Q*C

      srnamt = 'ZGEMQR';
      zgemqr('L', 'N', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

      // Compute |Q*C - Q*C| / |C|

      zgemm('N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M );
      RESID = ZLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 3 ) = RESID / (EPS*max(1,M)*CNORM);
      } else {
         RESULT( 3 ) = ZERO;
      }

      // Copy C into CF again

      zlacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as QT*C

      srnamt = 'ZGEMQR';
      zgemqr('L', 'C', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

      // Compute |QT*C - QT*C| / |C|

      zgemm('C', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M );
      RESID = ZLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 4 ) = RESID / (EPS*max(1,M)*CNORM);
      } else {
         RESULT( 4 ) = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= M; J++) {
         zlarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = ZLANGE( '1', N, M, D, N, RWORK);
      zlacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*Q

      srnamt = 'ZGEMQR';
      zgemqr('R', 'N', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

      // Compute |D*Q - D*Q| / |D|

      zgemm('N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N );
      RESID = ZLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT( 5 ) = RESID / (EPS*max(1,M)*DNORM);
      } else {
         RESULT( 5 ) = ZERO;
      }

      // Copy D into DF again

      zlacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*QT

      zgemqr('R', 'C', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

      // Compute |D*QT - D*QT| / |D|

      zgemm('N', 'C', N, M, M, -ONE, D, N, Q, M, ONE, DF, N );
      RESID = ZLANGE( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 6 ) = RESID / (EPS*max(1,M)*DNORM);
      } else {
         RESULT( 6 ) = ZERO;
      }

      // Short and wide

      } else {
      zgelq(M, N, AF, M, TQUERY, -1, WORKQUERY, -1, INFO );
      TSIZE = INT( TQUERY( 1 ) );
      LWORK = INT( WORKQUERY( 1 ) );
      zgemlq('R', 'N', N, N, K, AF, M, TQUERY, TSIZE, Q, N, WORKQUERY, -1, INFO );
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      zgemlq('L', 'N', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      zgemlq('L', 'C', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      zgemlq('R', 'N', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      zgemlq('R', 'C', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      ALLOCATE ( T( TSIZE ) );
      ALLOCATE ( WORK( LWORK ) );
      srnamt = 'ZGELQ';
      zgelq(M, N, AF, M, T, TSIZE, WORK, LWORK, INFO );


      // Generate the n-by-n matrix Q

      zlaset('Full', N, N, CZERO, ONE, Q, N );
      srnamt = 'ZGEMLQ';
      zgemlq('R', 'N', N, N, K, AF, M, T, TSIZE, Q, N, WORK, LWORK, INFO );

      // Copy R

      zlaset('Full', M, N, CZERO, CZERO, LQ, L );
      zlacpy('Lower', M, N, AF, M, LQ, L );

      // Compute |L - A*Q'| / |A| and store in RESULT(1)

      zgemm('N', 'C', M, N, N, -ONE, A, M, Q, N, ONE, LQ, L );
      ANORM = ZLANGE( '1', M, N, A, M, RWORK );
      RESID = ZLANGE( '1', M, N, LQ, L, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = RESID / (EPS*max(1,N)*ANORM);
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      zlaset('Full', N, N, CZERO, ONE, LQ, L );
      zherk('U', 'C', N, N, DREAL(-ONE), Q, N, DREAL(ONE), LQ, L);
      RESID = ZLANSY( '1', 'Upper', N, LQ, L, RWORK );
      RESULT( 2 ) = RESID / (EPS*max(1,N));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= M; J++) {
         zlarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = ZLANGE( '1', N, M, D, N, RWORK);
      zlacpy('Full', N, M, D, N, DF, N );

      // Apply Q to C as Q*C

      zgemlq('L', 'N', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

      // Compute |Q*D - Q*D| / |D|

      zgemm('N', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N );
      RESID = ZLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT( 3 ) = RESID / (EPS*max(1,N)*DNORM);
      } else {
         RESULT( 3 ) = ZERO;
      }

      // Copy D into DF again

      zlacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as QT*D

      zgemlq('L', 'C', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

      // Compute |QT*D - QT*D| / |D|

      zgemm('C', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N );
      RESID = ZLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT( 4 ) = RESID / (EPS*max(1,N)*DNORM);
      } else {
         RESULT( 4 ) = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= N; J++) {
         zlarnv(2, ISEED, M, C( 1, J ) );
      }
      CNORM = ZLANGE( '1', M, N, C, M, RWORK);
      zlacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as C*Q

      zgemlq('R', 'N', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

      // Compute |C*Q - C*Q| / |C|

      zgemm('N', 'N', M, N, N, -ONE, C, M, Q, N, ONE, CF, M );
      RESID = ZLANGE( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 5 ) = RESID / (EPS*max(1,N)*CNORM);
      } else {
         RESULT( 5 ) = ZERO;
      }

      // Copy C into CF again

      zlacpy('Full', M, N, C, M, CF, M );

      // Apply Q to D as D*QT

      zgemlq('R', 'C', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

      // Compute |C*QT - C*QT| / |C|

      zgemm('N', 'C', M, N, N, -ONE, C, M, Q, N, ONE, CF, M );
      RESID = ZLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 6 ) = RESID / (EPS*max(1,N)*CNORM);
      } else {
         RESULT( 6 ) = ZERO;
      }

      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF);

      return;
      }
