      void dtsqr01(final int TSSW, final int M, final int N, final int MB, final int NB, final int RESULT,) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String            TSSW;
      int               M, N, MB, NB;
      // .. Return values ..
      double            RESULT(6);

// =====================================================================

      // ..
      // .. Local allocatable arrays
      double          , ALLOCATABLE :: AF(:,:), Q(:,:), R(:,:), RWORK(:), WORK( : ), T(:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:), LQ(:,:);

      // .. Parameters ..
      double           ONE, ZERO;
      const    ZERO = 0.0, ONE = 1.0 ;
      bool    TESTZEROS, TS;
      int     INFO, J, K, L, LWORK, TSIZE, MNB;
      double             ANORM, EPS, RESID, CNORM, DNORM;
      int                ISEED( 4 );
      double             TQUERY( 5 ), WORKQUERY( 1 );
      // ..
      // .. External Functions ..
      //- double           DLAMCH, DLANGE, DLANSY;
      //- bool     lsame;
      //- int     ILAENV;
      // EXTERNAL DLAMCH, DLANGE, DLANSY, lsame, ILAENV
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

      TS = lsame(TSSW, 'TS');

      // TEST MATRICES WITH HALF OF MATRIX BEING ZEROS

      TESTZEROS = false;

      EPS = dlamch( 'Epsilon' );
      K = min(M,N);
      L = max(M,N,1);
      MNB = max( MB, NB);
      LWORK = max(3,L)*MNB;

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(L,L), R(M,L), RWORK(L), C(M,N), CF(M,N), D(N,M), DF(N,M), LQ(L,N) );

      // Put random numbers into A and copy to AF

      for (J = 1; J <= N; J++) {
         dlarnv(2, ISEED, M, A( 1, J ) );
      }
      if (TESTZEROS) {
         if (M >= 4) {
            for (J = 1; J <= N; J++) {
               dlarnv(2, ISEED, M/2, A( M/4, J ) );
            }
         }
      }
      dlacpy('Full', M, N, A, M, AF, M );

      if (TS) {

      // Factor the matrix A in the array AF.

      dgeqr(M, N, AF, M, TQUERY, -1, WORKQUERY, -1, INFO );
      TSIZE = INT( TQUERY( 1 ) );
      LWORK = INT( WORKQUERY( 1 ) );
      dgemqr('L', 'N', M, M, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      dgemqr('L', 'N', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      dgemqr('L', 'T', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      dgemqr('R', 'N', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      dgemqr('R', 'T', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      ALLOCATE ( T( TSIZE ) );
      ALLOCATE ( WORK( LWORK ) );
      srnamt = 'DGEQR';
      dgeqr(M, N, AF, M, T, TSIZE, WORK, LWORK, INFO );

      // Generate the m-by-m matrix Q

      dlaset('Full', M, M, ZERO, ONE, Q, M );
      srnamt = 'DGEMQR';
      dgemqr('L', 'N', M, M, K, AF, M, T, TSIZE, Q, M, WORK, LWORK, INFO );

      // Copy R

      dlaset('Full', M, N, ZERO, ZERO, R, M );
      dlacpy('Upper', M, N, AF, M, R, M );

      // Compute |R - Q'*A| / |A| and store in RESULT(1)

      dgemm('T', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M );
      ANORM = dlange( '1', M, N, A, M, RWORK );
      RESID = dlange( '1', M, N, R, M, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = RESID / (EPS*max(1,M)*ANORM);
      } else {
         RESULT[1] = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      dlaset('Full', M, M, ZERO, ONE, R, M );
      dsyrk('U', 'C', M, M, -ONE, Q, M, ONE, R, M );
      RESID = dlansy( '1', 'Upper', M, R, M, RWORK );
      RESULT[2] = RESID / (EPS*max(1,M));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= N; J++) {
         dlarnv(2, ISEED, M, C( 1, J ) );
      }
      CNORM = dlange( '1', M, N, C, M, RWORK);
      dlacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as Q*C

      srnamt = 'DGEMQR';
      dgemqr('L', 'N', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

      // Compute |Q*C - Q*C| / |C|

      dgemm('N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M );
      RESID = dlange( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[3] = RESID / (EPS*max(1,M)*CNORM);
      } else {
         RESULT[3] = ZERO;
      }

      // Copy C into CF again

      dlacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as QT*C

      srnamt = 'DGEMQR';
      dgemqr('L', 'T', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

      // Compute |QT*C - QT*C| / |C|

      dgemm('T', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M );
      RESID = dlange( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[4] = RESID / (EPS*max(1,M)*CNORM);
      } else {
         RESULT[4] = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= M; J++) {
         dlarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = dlange( '1', N, M, D, N, RWORK);
      dlacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*Q

      srnamt = 'DGEMQR';
      dgemqr('R', 'N', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

      // Compute |D*Q - D*Q| / |D|

      dgemm('N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N );
      RESID = dlange( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT[5] = RESID / (EPS*max(1,M)*DNORM);
      } else {
         RESULT[5] = ZERO;
      }

      // Copy D into DF again

      dlacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*QT

      dgemqr('R', 'T', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

      // Compute |D*QT - D*QT| / |D|

      dgemm('N', 'T', N, M, M, -ONE, D, N, Q, M, ONE, DF, N );
      RESID = dlange( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[6] = RESID / (EPS*max(1,M)*DNORM);
      } else {
         RESULT[6] = ZERO;
      }

      // Short and wide

      } else {
      dgelq(M, N, AF, M, TQUERY, -1, WORKQUERY, -1, INFO );
      TSIZE = INT( TQUERY( 1 ) );
      LWORK = INT( WORKQUERY( 1 ) );
      dgemlq('R', 'N', N, N, K, AF, M, TQUERY, TSIZE, Q, N, WORKQUERY, -1, INFO );
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      dgemlq('L', 'N', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      dgemlq('L', 'T', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      dgemlq('R', 'N', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      dgemlq('R', 'T', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );
      ALLOCATE ( T( TSIZE ) );
      ALLOCATE ( WORK( LWORK ) );
      srnamt = 'DGELQ';
      dgelq(M, N, AF, M, T, TSIZE, WORK, LWORK, INFO );


      // Generate the n-by-n matrix Q

      dlaset('Full', N, N, ZERO, ONE, Q, N );
      srnamt = 'DGEMLQ';
      dgemlq('R', 'N', N, N, K, AF, M, T, TSIZE, Q, N, WORK, LWORK, INFO );

      // Copy R

      dlaset('Full', M, N, ZERO, ZERO, LQ, L );
      dlacpy('Lower', M, N, AF, M, LQ, L );

      // Compute |L - A*Q'| / |A| and store in RESULT(1)

      dgemm('N', 'T', M, N, N, -ONE, A, M, Q, N, ONE, LQ, L );
      ANORM = dlange( '1', M, N, A, M, RWORK );
      RESID = dlange( '1', M, N, LQ, L, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = RESID / (EPS*max(1,N)*ANORM);
      } else {
         RESULT[1] = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      dlaset('Full', N, N, ZERO, ONE, LQ, L );
      dsyrk('U', 'C', N, N, -ONE, Q, N, ONE, LQ, L );
      RESID = dlansy( '1', 'Upper', N, LQ, L, RWORK );
      RESULT[2] = RESID / (EPS*max(1,N));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= M; J++) {
         dlarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = dlange( '1', N, M, D, N, RWORK);
      dlacpy('Full', N, M, D, N, DF, N );

      // Apply Q to C as Q*C

      dgemlq('L', 'N', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

      // Compute |Q*D - Q*D| / |D|

      dgemm('N', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N );
      RESID = dlange( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT[3] = RESID / (EPS*max(1,N)*DNORM);
      } else {
         RESULT[3] = ZERO;
      }

      // Copy D into DF again

      dlacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as QT*D

      dgemlq('L', 'T', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

      // Compute |QT*D - QT*D| / |D|

      dgemm('T', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N );
      RESID = dlange( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT[4] = RESID / (EPS*max(1,N)*DNORM);
      } else {
         RESULT[4] = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= N; J++) {
         dlarnv(2, ISEED, M, C( 1, J ) );
      }
      CNORM = dlange( '1', M, N, C, M, RWORK);
      dlacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as C*Q

      dgemlq('R', 'N', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

      // Compute |C*Q - C*Q| / |C|

      dgemm('N', 'N', M, N, N, -ONE, C, M, Q, N, ONE, CF, M );
      RESID = dlange( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[5] = RESID / (EPS*max(1,N)*CNORM);
      } else {
         RESULT[5] = ZERO;
      }

      // Copy C into CF again

      dlacpy('Full', M, N, C, M, CF, M );

      // Apply Q to D as D*QT

      dgemlq('R', 'T', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

      // Compute |C*QT - C*QT| / |C|

      dgemm('N', 'T', M, N, N, -ONE, C, M, Q, N, ONE, CF, M );
      RESID = dlange( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[6] = RESID / (EPS*max(1,N)*CNORM);
      } else {
         RESULT[6] = ZERO;
      }

      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF);

      }
