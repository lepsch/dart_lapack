      SUBROUTINE DQRT04(M,N,NB,RESULT);
      // IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     M, N, NB, LDT;
      // .. Return values ..
      double           RESULT(6);

*  =====================================================================

      // ..
      // .. Local allocatable arrays
      double          , ALLOCATABLE :: AF(:,:), Q(:,:), R(:,:), RWORK(:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:);

      // .. Parameters ..
      double           ONE, ZERO;
      const    ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int     INFO, J, K, L, LWORK;
      double             ANORM, EPS, RESID, CNORM, DNORM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      double           DLAMCH, DLANGE, DLANSY;
      bool     LSAME;
      // EXTERNAL DLAMCH, DLANGE, DLANSY, LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Data statements ..
      DATA ISEED / 1988, 1989, 1990, 1991 /;

      EPS = DLAMCH( 'Epsilon' );
      K = MIN(M,N);
      L = MAX(M,N);
      LWORK = MAX(2,L)*MAX(2,L)*NB;

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(M,M), R(M,L), RWORK(L), WORK(LWORK), T(NB,N), C(M,N), CF(M,N), D(N,M), DF(N,M) );

      // Put random numbers into A and copy to AF

      LDT=NB;
      for (J = 1; J <= N; J++) {
         dlarnv(2, ISEED, M, A( 1, J ) );
      }
      dlacpy('Full', M, N, A, M, AF, M );

      // Factor the matrix A in the array AF.

      dgeqrt(M, N, NB, AF, M, T, LDT, WORK, INFO );

      // Generate the m-by-m matrix Q

      dlaset('Full', M, M, ZERO, ONE, Q, M );
      dgemqrt('R', 'N', M, M, K, NB, AF, M, T, LDT, Q, M, WORK, INFO );

      // Copy R

      dlaset('Full', M, N, ZERO, ZERO, R, M );
      dlacpy('Upper', M, N, AF, M, R, M );

      // Compute |R - Q'*A| / |A| and store in RESULT(1)

      dgemm('T', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M );
      ANORM = DLANGE( '1', M, N, A, M, RWORK );
      RESID = DLANGE( '1', M, N, R, M, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = RESID / (EPS*MAX(1,M)*ANORM);
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      dlaset('Full', M, M, ZERO, ONE, R, M );
      dsyrk('U', 'C', M, M, -ONE, Q, M, ONE, R, M );
      RESID = DLANSY( '1', 'Upper', M, R, M, RWORK );
      RESULT( 2 ) = RESID / (EPS*MAX(1,M));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= N; J++) {
         dlarnv(2, ISEED, M, C( 1, J ) );
      }
      CNORM = DLANGE( '1', M, N, C, M, RWORK);
      dlacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as Q*C

      dgemqrt('L', 'N', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

      // Compute |Q*C - Q*C| / |C|

      dgemm('N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M );
      RESID = DLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 3 ) = RESID / (EPS*MAX(1,M)*CNORM);
      } else {
         RESULT( 3 ) = ZERO;
      }

      // Copy C into CF again

      dlacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as QT*C

      dgemqrt('L', 'T', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

      // Compute |QT*C - QT*C| / |C|

      dgemm('T', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M );
      RESID = DLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 4 ) = RESID / (EPS*MAX(1,M)*CNORM);
      } else {
         RESULT( 4 ) = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= M; J++) {
         dlarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = DLANGE( '1', N, M, D, N, RWORK);
      dlacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*Q

      dgemqrt('R', 'N', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

      // Compute |D*Q - D*Q| / |D|

      dgemm('N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N );
      RESID = DLANGE( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 5 ) = RESID / (EPS*MAX(1,M)*DNORM);
      } else {
         RESULT( 5 ) = ZERO;
      }

      // Copy D into DF again

      dlacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*QT

      dgemqrt('R', 'T', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

      // Compute |D*QT - D*QT| / |D|

      dgemm('N', 'T', N, M, M, -ONE, D, N, Q, M, ONE, DF, N );
      RESID = DLANGE( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 6 ) = RESID / (EPS*MAX(1,M)*DNORM);
      } else {
         RESULT( 6 ) = ZERO;
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF);

      RETURN;
      }
