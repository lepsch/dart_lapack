      SUBROUTINE ZLQT05(M,N,L,NB,RESULT);
      // IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     LWORK, M, N, L, NB, LDT;
      // .. Return values ..
      double           RESULT(6);

*  =====================================================================

      // ..
      // .. Local allocatable arrays
      COMPLEX*16, ALLOCATABLE :: AF(:,:), Q(:,:), R(:,:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:);
      double          , ALLOCATABLE :: RWORK(:);

      // .. Parameters ..
      double           ZERO;
      COMPLEX*16       ONE, CZERO;
      const    ZERO = 0.0, ONE = (1.0,0.0), CZERO=(0.0,0.0) ;
      // ..
      // .. Local Scalars ..
      int     INFO, J, K, N2, NP1,i;
      double             ANORM, EPS, RESID, CNORM, DNORM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      double           DLAMCH;
      double           ZLANGE, ZLANSY;
      bool     LSAME;
      // EXTERNAL DLAMCH, ZLANGE, ZLANSY, LSAME
      // ..
      // .. Data statements ..
      DATA ISEED / 1988, 1989, 1990, 1991 /;

      EPS = DLAMCH( 'Epsilon' );
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
      zlaset('Full', M, N2, CZERO, CZERO, A, M );
      zlaset('Full', NB, M, CZERO, CZERO, T, NB );
      for (J = 1; J <= M; J++) {
         zlarnv(2, ISEED, M-J+1, A( J, J ) );
      }
      if ( N > 0 ) {
         for (J = 1; J <= N-L; J++) {
            zlarnv(2, ISEED, M, A( 1, MIN(N+M,M+1) + J - 1 ) );
         }
      }
      if ( L > 0 ) {
         for (J = 1; J <= L; J++) {
            zlarnv(2, ISEED, M-J+1, A( J, MIN(N+M,N+M-L+1) + J - 1 ) );
         }
      }

      // Copy the matrix A to the array AF.

      zlacpy('Full', M, N2, A, M, AF, M );

      // Factor the matrix A in the array AF.

      ztplqt(M,N,L,NB,AF,M,AF(1,NP1),M,T,LDT,WORK,INFO);

      // Generate the (M+N)-by-(M+N) matrix Q by applying H to I

      zlaset('Full', N2, N2, CZERO, ONE, Q, N2 );
      zgemlqt('L', 'N', N2, N2, K, NB, AF, M, T, LDT, Q, N2, WORK, INFO );

      // Copy L

      zlaset('Full', N2, N2, CZERO, CZERO, R, N2 );
      zlacpy('Lower', M, N2, AF, M, R, N2 );

      // Compute |L - A*Q*C| / |A| and store in RESULT(1)

      zgemm('N', 'C', M, N2, N2, -ONE,  A, M, Q, N2, ONE, R, N2);
      ANORM = ZLANGE( '1', M, N2, A, M, RWORK );
      RESID = ZLANGE( '1', M, N2, R, N2, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = RESID / (EPS*ANORM*MAX(1,N2));
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute |I - Q*Q'| and store in RESULT(2)

      zlaset('Full', N2, N2, CZERO, ONE, R, N2 );
      zherk('U', 'N', N2, N2, DREAL(-ONE), Q, N2, DREAL(ONE), R, N2 );
      RESID = ZLANSY( '1', 'Upper', N2, R, N2, RWORK );
      RESULT( 2 ) = RESID / (EPS*MAX(1,N2));

      // Generate random m-by-n matrix C and a copy CF

      zlaset('Full', N2, M, CZERO, ONE, C, N2 );
      for (J = 1; J <= M; J++) {
         zlarnv(2, ISEED, N2, C( 1, J ) );
      }
      CNORM = ZLANGE( '1', N2, M, C, N2, RWORK);
      zlacpy('Full', N2, M, C, N2, CF, N2 );

      // Apply Q to C as Q*C

      ztpmlqt('L','N', N,M,K,L,NB,AF(1, NP1),M,T,LDT,CF,N2, CF(NP1,1),N2,WORK,INFO);

      // Compute |Q*C - Q*C| / |C|

      zgemm('N', 'N', N2, M, N2, -ONE, Q, N2, C, N2, ONE, CF, N2 );
      RESID = ZLANGE( '1', N2, M, CF, N2, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 3 ) = RESID / (EPS*MAX(1,N2)*CNORM);
      } else {
         RESULT( 3 ) = ZERO;
      }


      // Copy C into CF again

      zlacpy('Full', N2, M, C, N2, CF, N2 );

      // Apply Q to C as QT*C

      ztpmlqt('L','C',N,M,K,L,NB,AF(1,NP1),M,T,LDT,CF,N2, CF(NP1,1),N2,WORK,INFO);

      // Compute |QT*C - QT*C| / |C|

      zgemm('C','N',N2,M,N2,-ONE,Q,N2,C,N2,ONE,CF,N2);
      RESID = ZLANGE( '1', N2, M, CF, N2, RWORK );

      if ( CNORM > ZERO ) {
         RESULT( 4 ) = RESID / (EPS*MAX(1,N2)*CNORM);
      } else {
         RESULT( 4 ) = ZERO;
      }

      // Generate random m-by-n matrix D and a copy DF

      for (J = 1; J <= N2; J++) {
         zlarnv(2, ISEED, M, D( 1, J ) );
      }
      DNORM = ZLANGE( '1', M, N2, D, M, RWORK);
      zlacpy('Full', M, N2, D, M, DF, M );

      // Apply Q to D as D*Q

      ztpmlqt('R','N',M,N,K,L,NB,AF(1,NP1),M,T,LDT,DF,M, DF(1,NP1),M,WORK,INFO);

      // Compute |D*Q - D*Q| / |D|

      zgemm('N','N',M,N2,N2,-ONE,D,M,Q,N2,ONE,DF,M);
      RESID = ZLANGE('1',M, N2,DF,M,RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 5 ) = RESID / (EPS*MAX(1,N2)*DNORM);
      } else {
         RESULT( 5 ) = ZERO;
      }

      // Copy D into DF again

      zlacpy('Full',M,N2,D,M,DF,M );

      // Apply Q to D as D*QT

      ztpmlqt('R','C',M,N,K,L,NB,AF(1,NP1),M,T,LDT,DF,M, DF(1,NP1),M,WORK,INFO);


      // Compute |D*QT - D*QT| / |D|

      zgemm('N', 'C', M, N2, N2, -ONE, D, M, Q, N2, ONE, DF, M );
      RESID = ZLANGE( '1', M, N2, DF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 6 ) = RESID / (EPS*MAX(1,N2)*DNORM);
      } else {
         RESULT( 6 ) = ZERO;
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF);
      return;
      }
