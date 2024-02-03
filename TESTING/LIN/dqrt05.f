      SUBROUTINE DQRT05(M,N,L,NB,RESULT);
      // IMPLICIT NONE

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     LWORK, M, N, L, NB, LDT;
      // .. Return values ..
      double           RESULT(6);

// =====================================================================

      // ..
      // .. Local allocatable arrays
      double          , ALLOCATABLE :: AF(:,:), Q(:,:), R(:,:), RWORK(:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:);

      // .. Parameters ..
      double           ONE, ZERO;
      const    ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int     INFO, J, K, M2, NP1;
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
      // .. Data statements ..
      const ISEED = [ 1988, 1989, 1990, 1991 ];

      EPS = DLAMCH( 'Epsilon' );
      K = N;
      M2 = M+N;
      if ( M > 0 ) {
         NP1 = N+1;
      } else {
         NP1 = 1;
      }
      LWORK = M2*M2*NB;

      // Dynamically allocate all arrays

      ALLOCATE(A(M2,N),AF(M2,N),Q(M2,M2),R(M2,M2),RWORK(M2), WORK(LWORK),T(NB,N),C(M2,N),CF(M2,N), D(N,M2),DF(N,M2) );

      // Put random stuff into A

      LDT=NB;
      dlaset('Full', M2, N, ZERO, ZERO, A, M2 );
      dlaset('Full', NB, N, ZERO, ZERO, T, NB );
      for (J = 1; J <= N; J++) {
         dlarnv(2, ISEED, J, A( 1, J ) );
      }
      if ( M > 0 ) {
         for (J = 1; J <= N; J++) {
            dlarnv(2, ISEED, M-L, A( min(N+M,N+1), J ) );
         }
      }
      if ( L > 0 ) {
         for (J = 1; J <= N; J++) {
            dlarnv(2, ISEED, min(J,L), A( min(N+M,N+M-L+1), J ) );
         }
      }

      // Copy the matrix A to the array AF.

      dlacpy('Full', M2, N, A, M2, AF, M2 );

      // Factor the matrix A in the array AF.

      dtpqrt(M,N,L,NB,AF,M2,AF(NP1,1),M2,T,LDT,WORK,INFO);

      // Generate the (M+N)-by-(M+N) matrix Q by applying H to I

      dlaset('Full', M2, M2, ZERO, ONE, Q, M2 );
      dgemqrt('R', 'N', M2, M2, K, NB, AF, M2, T, LDT, Q, M2, WORK, INFO );

      // Copy R

      dlaset('Full', M2, N, ZERO, ZERO, R, M2 );
      dlacpy('Upper', M2, N, AF, M2, R, M2 );

      // Compute |R - Q'*A| / |A| and store in RESULT(1)

      dgemm('T', 'N', M2, N, M2, -ONE, Q, M2, A, M2, ONE, R, M2 );
      ANORM = DLANGE( '1', M2, N, A, M2, RWORK );
      RESID = DLANGE( '1', M2, N, R, M2, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = RESID / (EPS*ANORM*max(1,M2));
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      dlaset('Full', M2, M2, ZERO, ONE, R, M2 );
      dsyrk('U', 'C', M2, M2, -ONE, Q, M2, ONE, R, M2 );
      RESID = DLANSY( '1', 'Upper', M2, R, M2, RWORK );
      RESULT( 2 ) = RESID / (EPS*max(1,M2));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= N; J++) {
         dlarnv(2, ISEED, M2, C( 1, J ) );
      }
      CNORM = DLANGE( '1', M2, N, C, M2, RWORK);
      dlacpy('Full', M2, N, C, M2, CF, M2 );

      // Apply Q to C as Q*C

      dtpmqrt('L','N', M,N,K,L,NB,AF(NP1,1),M2,T,LDT,CF,M2, CF(NP1,1),M2,WORK,INFO);

      // Compute |Q*C - Q*C| / |C|

      dgemm('N', 'N', M2, N, M2, -ONE, Q, M2, C, M2, ONE, CF, M2 );
      RESID = DLANGE( '1', M2, N, CF, M2, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 3 ) = RESID / (EPS*max(1,M2)*CNORM);
      } else {
         RESULT( 3 ) = ZERO;
      }

      // Copy C into CF again

      dlacpy('Full', M2, N, C, M2, CF, M2 );

      // Apply Q to C as QT*C

      dtpmqrt('L','T',M,N,K,L,NB,AF(NP1,1),M2,T,LDT,CF,M2, CF(NP1,1),M2,WORK,INFO);

      // Compute |QT*C - QT*C| / |C|

      dgemm('T','N',M2,N,M2,-ONE,Q,M2,C,M2,ONE,CF,M2);
      RESID = DLANGE( '1', M2, N, CF, M2, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 4 ) = RESID / (EPS*max(1,M2)*CNORM);
      } else {
         RESULT( 4 ) = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= M2; J++) {
         dlarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = DLANGE( '1', N, M2, D, N, RWORK);
      dlacpy('Full', N, M2, D, N, DF, N );

      // Apply Q to D as D*Q

      dtpmqrt('R','N',N,M,N,L,NB,AF(NP1,1),M2,T,LDT,DF,N, DF(1,NP1),N,WORK,INFO);

      // Compute |D*Q - D*Q| / |D|

      dgemm('N','N',N,M2,M2,-ONE,D,N,Q,M2,ONE,DF,N);
      RESID = DLANGE('1',N, M2,DF,N,RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 5 ) = RESID / (EPS*max(1,M2)*DNORM);
      } else {
         RESULT( 5 ) = ZERO;
      }

      // Copy D into DF again

      dlacpy('Full',N,M2,D,N,DF,N );

      // Apply Q to D as D*QT

      dtpmqrt('R','T',N,M,N,L,NB,AF(NP1,1),M2,T,LDT,DF,N, DF(1,NP1),N,WORK,INFO);


      // Compute |D*QT - D*QT| / |D|

      dgemm('N', 'T', N, M2, M2, -ONE, D, N, Q, M2, ONE, DF, N );
      RESID = DLANGE( '1', N, M2, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 6 ) = RESID / (EPS*max(1,M2)*DNORM);
      } else {
         RESULT( 6 ) = ZERO;
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF);
      return;
      }
