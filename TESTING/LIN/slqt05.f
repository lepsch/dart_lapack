      SUBROUTINE SLQT05(M,N,L,NB,RESULT);
      // IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int      LWORK, M, N, L, NB, LDT;
      // .. Return values ..
      REAL     RESULT(6);

*  =====================================================================

      // ..
      // .. Local allocatable arrays
      REAL, ALLOCATABLE :: AF(:,:), Q(:,:), R(:,:), RWORK(:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:);

      // .. Parameters ..
      REAL ONE, ZERO;
      const    ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int         INFO, J, K, N2, NP1,i;
      REAL        ANORM, EPS, RESID, CNORM, DNORM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      REAL        SLAMCH, SLANGE, SLANSY;
      bool        LSAME;
      // EXTERNAL SLAMCH, SLANGE, SLANSY, LSAME
      // ..
      // .. Data statements ..
      DATA ISEED / 1988, 1989, 1990, 1991 /;

      EPS = SLAMCH( 'Epsilon' );
      K = M;
      N2 = M+N;
      if ( N > 0 ) {
         NP1 = M+1;
      } else {
         NP1 = 1;
      }
      LWORK = N2*N2*NB;

      // Dynamically allocate all arrays

      ALLOCATE(A(M,N2),AF(M,N2),Q(N2,N2),R(N2,N2),RWORK(N2), WORK(LWORK),T(NB,M),C(N2,M),CF(N2,M), D(M,N2),DF(M,N2) );

      // Put random stuff into A

      LDT=NB;
      slaset('Full', M, N2, ZERO, ZERO, A, M );
      slaset('Full', NB, M, ZERO, ZERO, T, NB );
      for (J = 1; J <= M; J++) {
         slarnv(2, ISEED, M-J+1, A( J, J ) );
      }
      if ( N > 0 ) {
         for (J = 1; J <= N-L; J++) {
            slarnv(2, ISEED, M, A( 1, MIN(N+M,M+1) + J - 1 ) );
         }
      }
      if ( L > 0 ) {
         for (J = 1; J <= L; J++) {
            slarnv(2, ISEED, M-J+1, A( J, MIN(N+M,N+M-L+1) + J - 1 ) );
         }
      }

      // Copy the matrix A to the array AF.

      slacpy('Full', M, N2, A, M, AF, M );

      // Factor the matrix A in the array AF.

      stplqt(M,N,L,NB,AF,M,AF(1,NP1),M,T,LDT,WORK,INFO);

      // Generate the (M+N)-by-(M+N) matrix Q by applying H to I

      slaset('Full', N2, N2, ZERO, ONE, Q, N2 );
      sgemlqt('L', 'N', N2, N2, K, NB, AF, M, T, LDT, Q, N2, WORK, INFO );

      // Copy L

      slaset('Full', N2, N2, ZERO, ZERO, R, N2 );
      slacpy('Lower', M, N2, AF, M, R, N2 );

      // Compute |L - A*Q*T| / |A| and store in RESULT(1)

      sgemm('N', 'T', M, N2, N2, -ONE,  A, M, Q, N2, ONE, R, N2);
      ANORM = SLANGE( '1', M, N2, A, M, RWORK );
      RESID = SLANGE( '1', M, N2, R, N2, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = RESID / (EPS*ANORM*MAX(1,N2));
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute |I - Q*Q'| and store in RESULT(2)

      slaset('Full', N2, N2, ZERO, ONE, R, N2 );
      ssyrk('U', 'N', N2, N2, -ONE, Q, N2, ONE, R, N2 );
      RESID = SLANSY( '1', 'Upper', N2, R, N2, RWORK );
      RESULT( 2 ) = RESID / (EPS*MAX(1,N2));

      // Generate random m-by-n matrix C and a copy CF

      slaset('Full', N2, M, ZERO, ONE, C, N2 );
      for (J = 1; J <= M; J++) {
         slarnv(2, ISEED, N2, C( 1, J ) );
      }
      CNORM = SLANGE( '1', N2, M, C, N2, RWORK);
      slacpy('Full', N2, M, C, N2, CF, N2 );

      // Apply Q to C as Q*C

      stpmlqt('L','N', N,M,K,L,NB,AF(1, NP1),M,T,LDT,CF,N2, CF(NP1,1),N2,WORK,INFO);

      // Compute |Q*C - Q*C| / |C|

      sgemm('N', 'N', N2, M, N2, -ONE, Q, N2, C, N2, ONE, CF, N2 );
      RESID = SLANGE( '1', N2, M, CF, N2, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 3 ) = RESID / (EPS*MAX(1,N2)*CNORM);
      } else {
         RESULT( 3 ) = ZERO;
      }


      // Copy C into CF again

      slacpy('Full', N2, M, C, N2, CF, N2 );

      // Apply Q to C as QT*C

      stpmlqt('L','T',N,M,K,L,NB,AF(1,NP1),M,T,LDT,CF,N2, CF(NP1,1),N2,WORK,INFO);

      // Compute |QT*C - QT*C| / |C|

      sgemm('T','N',N2,M,N2,-ONE,Q,N2,C,N2,ONE,CF,N2);
      RESID = SLANGE( '1', N2, M, CF, N2, RWORK );

      if ( CNORM > ZERO ) {
         RESULT( 4 ) = RESID / (EPS*MAX(1,N2)*CNORM);
      } else {
         RESULT( 4 ) = ZERO;
      }

      // Generate random m-by-n matrix D and a copy DF

      for (J = 1; J <= N2; J++) {
         slarnv(2, ISEED, M, D( 1, J ) );
      }
      DNORM = SLANGE( '1', M, N2, D, M, RWORK);
      slacpy('Full', M, N2, D, M, DF, M );

      // Apply Q to D as D*Q

      stpmlqt('R','N',M,N,K,L,NB,AF(1,NP1),M,T,LDT,DF,M, DF(1,NP1),M,WORK,INFO);

      // Compute |D*Q - D*Q| / |D|

      sgemm('N','N',M,N2,N2,-ONE,D,M,Q,N2,ONE,DF,M);
      RESID = SLANGE('1',M, N2,DF,M,RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 5 ) = RESID / (EPS*MAX(1,N2)*DNORM);
      } else {
         RESULT( 5 ) = ZERO;
      }

      // Copy D into DF again

      slacpy('Full',M,N2,D,M,DF,M );

      // Apply Q to D as D*QT

      stpmlqt('R','T',M,N,K,L,NB,AF(1,NP1),M,T,LDT,DF,M, DF(1,NP1),M,WORK,INFO);


      // Compute |D*QT - D*QT| / |D|

      sgemm('N', 'T', M, N2, N2, -ONE, D, M, Q, N2, ONE, DF, M );
      RESID = SLANGE( '1', M, N2, DF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 6 ) = RESID / (EPS*MAX(1,N2)*DNORM);
      } else {
         RESULT( 6 ) = ZERO;
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF);
      return;
      }
