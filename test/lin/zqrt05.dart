      void zqrt05(final int M, final int N, final int L, final int NB, final int RESULT,) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     LWORK, M, N, L, NB, LDT;
      // .. Return values ..
      double           RESULT(6);

// =====================================================================

      // ..
      // .. Local allocatable arrays
      Complex, ALLOCATABLE :: AF(:,:), Q(:,:), R(:,:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:);
      double          , ALLOCATABLE :: RWORK(:);

      // .. Parameters ..
      double           ZERO;
      Complex ONE, CZERO;
      const    ZERO = 0.0, ONE = (1.0,0.0), CZERO=(0.0,0.0) ;
      int     INFO, J, K, M2, NP1;
      double           ANORM, EPS, RESID, CNORM, DNORM;
      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      //- double           DLAMCH;
      //- double           ZLANGE, ZLANSY;
      //- bool     lsame;
      // EXTERNAL DLAMCH, ZLANGE, ZLANSY, lsame
      // ..
      // .. Data statements ..
      const ISEED = [ 1988, 1989, 1990, 1991 ];

      EPS = dlamch( 'Epsilon' );
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
      zlaset('Full', M2, N, CZERO, CZERO, A, M2 );
      zlaset('Full', NB, N, CZERO, CZERO, T, NB );
      for (J = 1; J <= N; J++) {
         zlarnv(2, ISEED, J, A( 1, J ) );
      }
      if ( M > 0 ) {
         for (J = 1; J <= N; J++) {
            zlarnv(2, ISEED, M-L, A( min(N+M,N+1), J ) );
         }
      }
      if ( L > 0 ) {
         for (J = 1; J <= N; J++) {
            zlarnv(2, ISEED, min(J,L), A( min(N+M,N+M-L+1), J ) );
         }
      }

      // Copy the matrix A to the array AF.

      zlacpy('Full', M2, N, A, M2, AF, M2 );

      // Factor the matrix A in the array AF.

      ztpqrt(M,N,L,NB,AF,M2,AF(NP1,1),M2,T,LDT,WORK,INFO);

      // Generate the (M+N)-by-(M+N) matrix Q by applying H to I

      zlaset('Full', M2, M2, CZERO, ONE, Q, M2 );
      zgemqrt('R', 'N', M2, M2, K, NB, AF, M2, T, LDT, Q, M2, WORK, INFO );

      // Copy R

      zlaset('Full', M2, N, CZERO, CZERO, R, M2 );
      zlacpy('Upper', M2, N, AF, M2, R, M2 );

      // Compute |R - Q'*A| / |A| and store in RESULT(1)

      zgemm('C', 'N', M2, N, M2, -ONE, Q, M2, A, M2, ONE, R, M2 );
      ANORM = ZLANGE( '1', M2, N, A, M2, RWORK );
      RESID = ZLANGE( '1', M2, N, R, M2, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = RESID / (EPS*ANORM*max(1,M2));
      } else {
         RESULT[1] = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      zlaset('Full', M2, M2, CZERO, ONE, R, M2 );
      zherk('U', 'C', M2, M2, DREAL(-ONE), Q, M2, DREAL(ONE), R, M2 );
      RESID = ZLANSY( '1', 'Upper', M2, R, M2, RWORK );
      RESULT[2] = RESID / (EPS*max(1,M2));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= N; J++) {
         zlarnv(2, ISEED, M2, C( 1, J ) );
      }
      CNORM = ZLANGE( '1', M2, N, C, M2, RWORK);
      zlacpy('Full', M2, N, C, M2, CF, M2 );

      // Apply Q to C as Q*C

      ztpmqrt('L','N', M,N,K,L,NB,AF(NP1,1),M2,T,LDT,CF,M2, CF(NP1,1),M2,WORK,INFO);

      // Compute |Q*C - Q*C| / |C|

      zgemm('N', 'N', M2, N, M2, -ONE, Q, M2, C, M2, ONE, CF, M2 );
      RESID = ZLANGE( '1', M2, N, CF, M2, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[3] = RESID / (EPS*max(1,M2)*CNORM);
      } else {
         RESULT[3] = ZERO;
      }

      // Copy C into CF again

      zlacpy('Full', M2, N, C, M2, CF, M2 );

      // Apply Q to C as QT*C

      ztpmqrt('L','C',M,N,K,L,NB,AF(NP1,1),M2,T,LDT,CF,M2, CF(NP1,1),M2,WORK,INFO);

      // Compute |QT*C - QT*C| / |C|

      zgemm('C','N',M2,N,M2,-ONE,Q,M2,C,M2,ONE,CF,M2);
      RESID = ZLANGE( '1', M2, N, CF, M2, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[4] = RESID / (EPS*max(1,M2)*CNORM);
      } else {
         RESULT[4] = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= M2; J++) {
         zlarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = ZLANGE( '1', N, M2, D, N, RWORK);
      zlacpy('Full', N, M2, D, N, DF, N );

      // Apply Q to D as D*Q

      ztpmqrt('R','N',N,M,N,L,NB,AF(NP1,1),M2,T,LDT,DF,N, DF(1,NP1),N,WORK,INFO);

      // Compute |D*Q - D*Q| / |D|

      zgemm('N','N',N,M2,M2,-ONE,D,N,Q,M2,ONE,DF,N);
      RESID = ZLANGE('1',N, M2,DF,N,RWORK );
      if ( CNORM > ZERO ) {
         RESULT[5] = RESID / (EPS*max(1,M2)*DNORM);
      } else {
         RESULT[5] = ZERO;
      }

      // Copy D into DF again

      zlacpy('Full',N,M2,D,N,DF,N );

      // Apply Q to D as D*QT

      ztpmqrt('R','C',N,M,N,L,NB,AF(NP1,1),M2,T,LDT,DF,N, DF(1,NP1),N,WORK,INFO);


      // Compute |D*QT - D*QT| / |D|

      zgemm('N', 'C', N, M2, M2, -ONE, D, N, Q, M2, ONE, DF, N );
      RESID = ZLANGE( '1', N, M2, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[6] = RESID / (EPS*max(1,M2)*DNORM);
      } else {
         RESULT[6] = ZERO;
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF);
      }
