      SUBROUTINE SLQT04(M,N,NB,RESULT);
      // IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     M, N, NB, LDT;
      // .. Return values ..
      REAL RESULT(6);

*  =====================================================================

      // ..
      // .. Local allocatable arrays
      REAL, ALLOCATABLE :: AF(:,:), Q(:,:), L(:,:), RWORK(:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:);

      // .. Parameters ..
      REAL ONE, ZERO;
      const    ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int     INFO, J, K, LL, LWORK;
      REAL   ANORM, EPS, RESID, CNORM, DNORM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      REAL SLAMCH, SLANGE, SLANSY;
      bool     LSAME;
      // EXTERNAL SLAMCH, SLANGE, SLANSY, LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Data statements ..
      DATA ISEED / 1988, 1989, 1990, 1991 /;

      EPS = SLAMCH( 'Epsilon' );
      K = MIN(M,N);
      LL = MAX(M,N);
      LWORK = MAX(2,LL)*MAX(2,LL)*NB;

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(N,N), L(LL,N), RWORK(LL), WORK(LWORK), T(NB,N), C(M,N), CF(M,N), D(N,M), DF(N,M) );

      // Put random numbers into A and copy to AF

      LDT=NB;
      for (J = 1; J <= N; J++) {
         slarnv(2, ISEED, M, A( 1, J ) );
      }
      slacpy('Full', M, N, A, M, AF, M );

      // Factor the matrix A in the array AF.

      sgelqt(M, N, NB, AF, M, T, LDT, WORK, INFO );

      // Generate the n-by-n matrix Q

      slaset('Full', N, N, ZERO, ONE, Q, N );
      sgemlqt('R', 'N', N, N, K, NB, AF, M, T, LDT, Q, N, WORK, INFO );

      // Copy R

      slaset('Full', M, N, ZERO, ZERO, L, LL );
      slacpy('Lower', M, N, AF, M, L, LL );

      // Compute |L - A*Q'| / |A| and store in RESULT(1)

      sgemm('N', 'T', M, N, N, -ONE, A, M, Q, N, ONE, L, LL );
      ANORM = SLANGE( '1', M, N, A, M, RWORK );
      RESID = SLANGE( '1', M, N, L, LL, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = RESID / (EPS*MAX(1,M)*ANORM);
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      slaset('Full', N, N, ZERO, ONE, L, LL );
      ssyrk('U', 'C', N, N, -ONE, Q, N, ONE, L, LL );
      RESID = SLANSY( '1', 'Upper', N, L, LL, RWORK );
      RESULT( 2 ) = RESID / (EPS*MAX(1,N));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= M; J++) {
         slarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = SLANGE( '1', N, M, D, N, RWORK);
      slacpy('Full', N, M, D, N, DF, N );

      // Apply Q to C as Q*C

      sgemlqt('L', 'N', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

      // Compute |Q*D - Q*D| / |D|

      sgemm('N', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N );
      RESID = SLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT( 3 ) = RESID / (EPS*MAX(1,M)*DNORM);
      } else {
         RESULT( 3 ) = ZERO;
      }

      // Copy D into DF again

      slacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as QT*D

      sgemlqt('L', 'T', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

      // Compute |QT*D - QT*D| / |D|

      sgemm('T', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N );
      RESID = SLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT( 4 ) = RESID / (EPS*MAX(1,M)*DNORM);
      } else {
         RESULT( 4 ) = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= N; J++) {
         slarnv(2, ISEED, M, C( 1, J ) );
      }
      CNORM = SLANGE( '1', M, N, C, M, RWORK);
      slacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as C*Q

      sgemlqt('R', 'N', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

      // Compute |C*Q - C*Q| / |C|

      sgemm('N', 'N', M, N, N, -ONE, C, M, Q, N, ONE, CF, M );
      RESID = SLANGE( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 5 ) = RESID / (EPS*MAX(1,M)*DNORM);
      } else {
         RESULT( 5 ) = ZERO;
      }

      // Copy C into CF again

      slacpy('Full', M, N, C, M, CF, M );

      // Apply Q to D as D*QT

      sgemlqt('R', 'T', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

      // Compute |C*QT - C*QT| / |C|

      sgemm('N', 'T', M, N, N, -ONE, C, M, Q, N, ONE, CF, M );
      RESID = SLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 6 ) = RESID / (EPS*MAX(1,M)*DNORM);
      } else {
         RESULT( 6 ) = ZERO;
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, L, RWORK, WORK, T, C, D, CF, DF);

      RETURN;
      }
