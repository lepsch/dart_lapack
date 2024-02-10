      void clqt05(M,N,L,NB,RESULT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     LWORK, M, N, L, NB, LDT;
      // .. Return values ..
      double RESULT(6);

// =====================================================================

      // ..
      // .. Local allocatable arrays
      Complex, ALLOCATABLE :: AF(:,:), Q(:,:), R(:,:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:);
      double, ALLOCATABLE :: RWORK(:);

      // .. Parameters ..
      double ZERO;
      Complex       ONE, CZERO;
      const    ZERO = 0.0, ONE = (1.0,0.0), CZERO=(0.0,0.0) ;
      int     INFO, J, K, N2, NP1,i;
      double   ANORM, EPS, RESID, CNORM, DNORM;
      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      //- REAL SLAMCH;
      //- REAL CLANGE, CLANSY;
      //- bool     lsame;
      // EXTERNAL SLAMCH, CLANGE, CLANSY, lsame
      // ..
      // .. Data statements ..
      const ISEED = [ 1988, 1989, 1990, 1991 ];

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
      claset('Full', M, N2, CZERO, CZERO, A, M );
      claset('Full', NB, M, CZERO, CZERO, T, NB );
      for (J = 1; J <= M; J++) {
         clarnv(2, ISEED, M-J+1, A( J, J ) );
      }
      if ( N > 0 ) {
         for (J = 1; J <= N-L; J++) {
            clarnv(2, ISEED, M, A( 1, min(N+M,M+1) + J - 1 ) );
         }
      }
      if ( L > 0 ) {
         for (J = 1; J <= L; J++) {
            clarnv(2, ISEED, M-J+1, A( J, min(N+M,N+M-L+1) + J - 1 ) );
         }
      }

      // Copy the matrix A to the array AF.

      clacpy('Full', M, N2, A, M, AF, M );

      // Factor the matrix A in the array AF.

      ctplqt(M,N,L,NB,AF,M,AF(1,NP1),M,T,LDT,WORK,INFO);

      // Generate the (M+N)-by-(M+N) matrix Q by applying H to I

      claset('Full', N2, N2, CZERO, ONE, Q, N2 );
      cgemlqt('L', 'N', N2, N2, K, NB, AF, M, T, LDT, Q, N2, WORK, INFO );

      // Copy L

      claset('Full', N2, N2, CZERO, CZERO, R, N2 );
      clacpy('Lower', M, N2, AF, M, R, N2 );

      // Compute |L - A*Q*C| / |A| and store in RESULT(1)

      cgemm('N', 'C', M, N2, N2, -ONE,  A, M, Q, N2, ONE, R, N2);
      ANORM = CLANGE( '1', M, N2, A, M, RWORK );
      RESID = CLANGE( '1', M, N2, R, N2, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = RESID / (EPS*ANORM*max(1,N2));
      } else {
         RESULT[1] = ZERO;
      }

      // Compute |I - Q*Q'| and store in RESULT(2)

      claset('Full', N2, N2, CZERO, ONE, R, N2 );
      cherk('U', 'N', N2, N2, REAL(-ONE), Q, N2, double(ONE), R, N2 );
      RESID = CLANSY( '1', 'Upper', N2, R, N2, RWORK );
      RESULT[2] = RESID / (EPS*max(1,N2));

      // Generate random m-by-n matrix C and a copy CF

      claset('Full', N2, M, CZERO, ONE, C, N2 );
      for (J = 1; J <= M; J++) {
         clarnv(2, ISEED, N2, C( 1, J ) );
      }
      CNORM = CLANGE( '1', N2, M, C, N2, RWORK);
      clacpy('Full', N2, M, C, N2, CF, N2 );

      // Apply Q to C as Q*C

      ctpmlqt('L','N', N,M,K,L,NB,AF(1, NP1),M,T,LDT,CF,N2, CF(NP1,1),N2,WORK,INFO);

      // Compute |Q*C - Q*C| / |C|

      cgemm('N', 'N', N2, M, N2, -ONE, Q, N2, C, N2, ONE, CF, N2 );
      RESID = CLANGE( '1', N2, M, CF, N2, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[3] = RESID / (EPS*max(1,N2)*CNORM);
      } else {
         RESULT[3] = ZERO;
      }


      // Copy C into CF again

      clacpy('Full', N2, M, C, N2, CF, N2 );

      // Apply Q to C as QT*C

      ctpmlqt('L','C',N,M,K,L,NB,AF(1,NP1),M,T,LDT,CF,N2, CF(NP1,1),N2,WORK,INFO);

      // Compute |QT*C - QT*C| / |C|

      cgemm('C','N',N2,M,N2,-ONE,Q,N2,C,N2,ONE,CF,N2);
      RESID = CLANGE( '1', N2, M, CF, N2, RWORK );

      if ( CNORM > ZERO ) {
         RESULT[4] = RESID / (EPS*max(1,N2)*CNORM);
      } else {
         RESULT[4] = ZERO;
      }

      // Generate random m-by-n matrix D and a copy DF

      for (J = 1; J <= N2; J++) {
         clarnv(2, ISEED, M, D( 1, J ) );
      }
      DNORM = CLANGE( '1', M, N2, D, M, RWORK);
      clacpy('Full', M, N2, D, M, DF, M );

      // Apply Q to D as D*Q

      ctpmlqt('R','N',M,N,K,L,NB,AF(1,NP1),M,T,LDT,DF,M, DF(1,NP1),M,WORK,INFO);

      // Compute |D*Q - D*Q| / |D|

      cgemm('N','N',M,N2,N2,-ONE,D,M,Q,N2,ONE,DF,M);
      RESID = CLANGE('1',M, N2,DF,M,RWORK );
      if ( CNORM > ZERO ) {
         RESULT[5] = RESID / (EPS*max(1,N2)*DNORM);
      } else {
         RESULT[5] = ZERO;
      }

      // Copy D into DF again

      clacpy('Full',M,N2,D,M,DF,M );

      // Apply Q to D as D*QT

      ctpmlqt('R','C',M,N,K,L,NB,AF(1,NP1),M,T,LDT,DF,M, DF(1,NP1),M,WORK,INFO);


      // Compute |D*QT - D*QT| / |D|

      cgemm('N', 'C', M, N2, N2, -ONE, D, M, Q, N2, ONE, DF, M );
      RESID = CLANGE( '1', M, N2, DF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[6] = RESID / (EPS*max(1,N2)*DNORM);
      } else {
         RESULT[6] = ZERO;
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF);
      }
