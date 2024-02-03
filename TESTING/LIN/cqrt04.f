      SUBROUTINE CQRT04(M,N,NB,RESULT);
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
      COMPLEX, ALLOCATABLE :: AF(:,:), Q(:,:), R(:,:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:);
      REAL, ALLOCATABLE :: RWORK(:);

      // .. Parameters ..
      REAL ZERO;
      COMPLEX ONE, CZERO;
      const    ZERO = 0.0, ONE = (1.0,0.0), CZERO=(0.0,0.0) ;
      // ..
      // .. Local Scalars ..
      int     INFO, J, K, L, LWORK;
      REAL   ANORM, EPS, RESID, CNORM, DNORM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      REAL SLAMCH;
      REAL CLANGE, CLANSY;
      bool     LSAME;
      // EXTERNAL SLAMCH, CLANGE, CLANSY, LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Data statements ..
      DATA ISEED / 1988, 1989, 1990, 1991 /;

      EPS = SLAMCH( 'Epsilon' );
      K = MIN(M,N);
      L = MAX(M,N);
      LWORK = MAX(2,L)*MAX(2,L)*NB;

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(M,M), R(M,L), RWORK(L), WORK(LWORK), T(NB,N), C(M,N), CF(M,N), D(N,M), DF(N,M) );

      // Put random numbers into A and copy to AF

      LDT=NB;
      for (J = 1; J <= N; J++) {
         clarnv(2, ISEED, M, A( 1, J ) );
      }
      clacpy('Full', M, N, A, M, AF, M );

      // Factor the matrix A in the array AF.

      cgeqrt(M, N, NB, AF, M, T, LDT, WORK, INFO );

      // Generate the m-by-m matrix Q

      claset('Full', M, M, CZERO, ONE, Q, M );
      cgemqrt('R', 'N', M, M, K, NB, AF, M, T, LDT, Q, M, WORK, INFO );

      // Copy R

      claset('Full', M, N, CZERO, CZERO, R, M );
      clacpy('Upper', M, N, AF, M, R, M );

      // Compute |R - Q'*A| / |A| and store in RESULT(1)

      cgemm('C', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M );
      ANORM = CLANGE( '1', M, N, A, M, RWORK );
      RESID = CLANGE( '1', M, N, R, M, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = RESID / (EPS*MAX(1,M)*ANORM);
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      claset('Full', M, M, CZERO, ONE, R, M );
      cherk('U', 'C', M, M, REAL(-ONE), Q, M, REAL(ONE), R, M );
      RESID = CLANSY( '1', 'Upper', M, R, M, RWORK );
      RESULT( 2 ) = RESID / (EPS*MAX(1,M));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= N; J++) {
         clarnv(2, ISEED, M, C( 1, J ) );
      }
      CNORM = CLANGE( '1', M, N, C, M, RWORK);
      clacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as Q*C

      cgemqrt('L', 'N', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

      // Compute |Q*C - Q*C| / |C|

      cgemm('N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M );
      RESID = CLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 3 ) = RESID / (EPS*MAX(1,M)*CNORM);
      } else {
         RESULT( 3 ) = ZERO;
      }

      // Copy C into CF again

      clacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as QT*C

      cgemqrt('L', 'C', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO);

      // Compute |QT*C - QT*C| / |C|

      cgemm('C', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M );
      RESID = CLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 4 ) = RESID / (EPS*MAX(1,M)*CNORM);
      } else {
         RESULT( 4 ) = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= M; J++) {
         clarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = CLANGE( '1', N, M, D, N, RWORK);
      clacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*Q

      cgemqrt('R', 'N', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

      // Compute |D*Q - D*Q| / |D|

      cgemm('N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N );
      RESID = CLANGE( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 5 ) = RESID / (EPS*MAX(1,M)*DNORM);
      } else {
         RESULT( 5 ) = ZERO;
      }

      // Copy D into DF again

      clacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*QT

      cgemqrt('R', 'C', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO);

      // Compute |D*QT - D*QT| / |D|

      cgemm('N', 'C', N, M, M, -ONE, D, N, Q, M, ONE, DF, N );
      RESID = CLANGE( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 6 ) = RESID / (EPS*MAX(1,M)*DNORM);
      } else {
         RESULT( 6 ) = ZERO;
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF);

      RETURN;
      }
