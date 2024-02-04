      void dlqt04(M,N,NB,RESULT) {
      // IMPLICIT NONE

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     M, N, NB, LDT;
      // .. Return values ..
      double           RESULT(6);

// =====================================================================

      // ..
      // .. Local allocatable arrays
      double          , ALLOCATABLE :: AF(:,:), Q(:,:), L(:,:), RWORK(:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:);

      // .. Parameters ..
      double           ONE, ZERO;
      const    ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int     INFO, J, K, LL, LWORK;
      double             ANORM, EPS, RESID, CNORM, DNORM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      //- double           DLAMCH, DLANGE, DLANSY;
      //- bool     LSAME;
      // EXTERNAL DLAMCH, DLANGE, DLANSY, LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Data statements ..
      const ISEED = [ 1988, 1989, 1990, 1991 ];

      EPS = DLAMCH( 'Epsilon' );
      K = min(M,N);
      LL = max(M,N);
      LWORK = max(2,LL)*max(2,LL)*NB;

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(N,N), L(LL,N), RWORK(LL), WORK(LWORK), T(NB,N), C(M,N), CF(M,N), D(N,M), DF(N,M) );

      // Put random numbers into A and copy to AF

      LDT=NB;
      for (J = 1; J <= N; J++) {
         dlarnv(2, ISEED, M, A( 1, J ) );
      }
      dlacpy('Full', M, N, A, M, AF, M );

      // Factor the matrix A in the array AF.

      dgelqt(M, N, NB, AF, M, T, LDT, WORK, INFO );

      // Generate the n-by-n matrix Q

      dlaset('Full', N, N, ZERO, ONE, Q, N );
      dgemlqt('R', 'N', N, N, K, NB, AF, M, T, LDT, Q, N, WORK, INFO );

      // Copy R

      dlaset('Full', M, N, ZERO, ZERO, L, LL );
      dlacpy('Lower', M, N, AF, M, L, LL );

      // Compute |L - A*Q'| / |A| and store in RESULT(1)

      dgemm('N', 'T', M, N, N, -ONE, A, M, Q, N, ONE, L, LL );
      ANORM = DLANGE( '1', M, N, A, M, RWORK );
      RESID = DLANGE( '1', M, N, L, LL, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = RESID / (EPS*max(1,M)*ANORM);
      } else {
         RESULT[1] = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      dlaset('Full', N, N, ZERO, ONE, L, LL );
      dsyrk('U', 'C', N, N, -ONE, Q, N, ONE, L, LL );
      RESID = DLANSY( '1', 'Upper', N, L, LL, RWORK );
      RESULT[2] = RESID / (EPS*max(1,N));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= M; J++) {
         dlarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = DLANGE( '1', N, M, D, N, RWORK);
      dlacpy('Full', N, M, D, N, DF, N );

      // Apply Q to C as Q*C

      dgemlqt('L', 'N', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

      // Compute |Q*D - Q*D| / |D|

      dgemm('N', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N );
      RESID = DLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT[3] = RESID / (EPS*max(1,M)*DNORM);
      } else {
         RESULT[3] = ZERO;
      }

      // Copy D into DF again

      dlacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as QT*D

      dgemlqt('L', 'T', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

      // Compute |QT*D - QT*D| / |D|

      dgemm('T', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N );
      RESID = DLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT[4] = RESID / (EPS*max(1,M)*DNORM);
      } else {
         RESULT[4] = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= N; J++) {
         dlarnv(2, ISEED, M, C( 1, J ) );
      }
      CNORM = DLANGE( '1', M, N, C, M, RWORK);
      dlacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as C*Q

      dgemlqt('R', 'N', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

      // Compute |C*Q - C*Q| / |C|

      dgemm('N', 'N', M, N, N, -ONE, C, M, Q, N, ONE, CF, M );
      RESID = DLANGE( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[5] = RESID / (EPS*max(1,M)*DNORM);
      } else {
         RESULT[5] = ZERO;
      }

      // Copy C into CF again

      dlacpy('Full', M, N, C, M, CF, M );

      // Apply Q to D as D*QT

      dgemlqt('R', 'T', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

      // Compute |C*QT - C*QT| / |C|

      dgemm('N', 'T', M, N, N, -ONE, C, M, Q, N, ONE, CF, M );
      RESID = DLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[6] = RESID / (EPS*max(1,M)*DNORM);
      } else {
         RESULT[6] = ZERO;
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, L, RWORK, WORK, T, C, D, CF, DF);

      return;
      }
