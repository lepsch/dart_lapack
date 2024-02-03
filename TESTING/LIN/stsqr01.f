      SUBROUTINE STSQR01(TSSW, M, N, MB, NB, RESULT);
      // IMPLICIT NONE

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String            TSSW;
      int               M, N, MB, NB;
      // .. Return values ..
      REAL              RESULT(6);

// =====================================================================

      // ..
      // .. Local allocatable arrays
      REAL, ALLOCATABLE :: AF(:,:), Q(:,:), R(:,:), RWORK(:), WORK( : ), T(:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:), LQ(:,:);

      // .. Parameters ..
      REAL     ONE, ZERO;
      const    ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool    TESTZEROS, TS;
      int     INFO, J, K, L, LWORK, TSIZE, MNB;
      REAL   ANORM, EPS, RESID, CNORM, DNORM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      REAL               TQUERY( 5 ), WORKQUERY( 1 );
      // ..
      // .. External Functions ..
      REAL     SLAMCH, SLANGE, SLANSY;
      bool     LSAME;
      int      ILAENV;
      // EXTERNAL SLAMCH, SLARNV, SLANGE, SLANSY, LSAME, ILAENV
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
      DATA ISEED / 1988, 1989, 1990, 1991 /;

      // TEST TALL SKINNY OR SHORT WIDE

      TS = LSAME(TSSW, 'TS');

      // TEST MATRICES WITH HALF OF MATRIX BEING ZEROS

      TESTZEROS = false;

      EPS = SLAMCH( 'Epsilon' );
      K = MIN(M,N);
      L = MAX(M,N,1);
      MNB = MAX ( MB, NB);
      LWORK = MAX(3,L)*MNB;

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(L,L), R(M,L), RWORK(L), C(M,N), CF(M,N), D(N,M), DF(N,M), LQ(L,N) );

      // Put random numbers into A and copy to AF

      for (J = 1; J <= N; J++) {
         slarnv(2, ISEED, M, A( 1, J ) );
      }
      if (TESTZEROS) {
         if (M >= 4) {
            for (J = 1; J <= N; J++) {
               slarnv(2, ISEED, M/2, A( M/4, J ) );
            }
         }
      }
      slacpy('Full', M, N, A, M, AF, M );

      if (TS) {

      // Factor the matrix A in the array AF.

      sgeqr(M, N, AF, M, TQUERY, -1, WORKQUERY, -1, INFO );
      TSIZE = INT( TQUERY( 1 ) );
      LWORK = INT( WORKQUERY( 1 ) );
      sgemqr('L', 'N', M, M, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) );
      sgemqr('L', 'N', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) );
      sgemqr('L', 'T', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) );
      sgemqr('R', 'N', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) );
      sgemqr('R', 'T', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) );
      ALLOCATE ( T( TSIZE ) );
      ALLOCATE ( WORK( LWORK ) );
      srnamt = 'SGEQR';
      sgeqr(M, N, AF, M, T, TSIZE, WORK, LWORK, INFO );

      // Generate the m-by-m matrix Q

      slaset('Full', M, M, ZERO, ONE, Q, M );
      srnamt = 'SGEMQR';
      sgemqr('L', 'N', M, M, K, AF, M, T, TSIZE, Q, M, WORK, LWORK, INFO );

      // Copy R

      slaset('Full', M, N, ZERO, ZERO, R, M );
      slacpy('Upper', M, N, AF, M, R, M );

      // Compute |R - Q'*A| / |A| and store in RESULT(1)

      sgemm('T', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M );
      ANORM = SLANGE( '1', M, N, A, M, RWORK );
      RESID = SLANGE( '1', M, N, R, M, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = RESID / (EPS*MAX(1,M)*ANORM);
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      slaset('Full', M, M, ZERO, ONE, R, M );
      ssyrk('U', 'C', M, M, -ONE, Q, M, ONE, R, M );
      RESID = SLANSY( '1', 'Upper', M, R, M, RWORK );
      RESULT( 2 ) = RESID / (EPS*MAX(1,M));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= N; J++) {
         slarnv(2, ISEED, M, C( 1, J ) );
      }
      CNORM = SLANGE( '1', M, N, C, M, RWORK);
      slacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as Q*C

      srnamt = 'DGEQR';
      sgemqr('L', 'N', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

      // Compute |Q*C - Q*C| / |C|

      sgemm('N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M );
      RESID = SLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 3 ) = RESID / (EPS*MAX(1,M)*CNORM);
      } else {
         RESULT( 3 ) = ZERO;
      }

      // Copy C into CF again

      slacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as QT*C

      srnamt = 'DGEQR';
      sgemqr('L', 'T', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

      // Compute |QT*C - QT*C| / |C|

      sgemm('T', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M );
      RESID = SLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 4 ) = RESID / (EPS*MAX(1,M)*CNORM);
      } else {
         RESULT( 4 ) = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= M; J++) {
         slarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = SLANGE( '1', N, M, D, N, RWORK);
      slacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*Q

      srnamt = 'DGEQR';
      sgemqr('R', 'N', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

      // Compute |D*Q - D*Q| / |D|

      sgemm('N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N );
      RESID = SLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT( 5 ) = RESID / (EPS*MAX(1,M)*DNORM);
      } else {
         RESULT( 5 ) = ZERO;
      }

      // Copy D into DF again

      slacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*QT

      sgemqr('R', 'T', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

      // Compute |D*QT - D*QT| / |D|

      sgemm('N', 'T', N, M, M, -ONE, D, N, Q, M, ONE, DF, N );
      RESID = SLANGE( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 6 ) = RESID / (EPS*MAX(1,M)*DNORM);
      } else {
         RESULT( 6 ) = ZERO;
      }

      // Short and wide

      } else {
      sgelq(M, N, AF, M, TQUERY, -1, WORKQUERY, -1, INFO );
      TSIZE = INT( TQUERY( 1 ) );
      LWORK = INT( WORKQUERY( 1 ));
      sgemlq('R', 'N', N, N, K, AF, M, TQUERY, TSIZE, Q, N, WORKQUERY, -1, INFO );
      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) );
      sgemlq('L', 'N', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) );
      sgemlq('L', 'T', N, M, K, AF, M, TQUERY, TSIZE, DF, N, WORKQUERY, -1, INFO);
      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) );
      sgemlq('R', 'N', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) );
      sgemlq('R', 'T', M, N, K, AF, M, TQUERY, TSIZE, CF, M, WORKQUERY, -1, INFO);
      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) );
      ALLOCATE ( T( TSIZE ) );
      ALLOCATE ( WORK( LWORK ) );
      srnamt = 'SGELQ';
      sgelq(M, N, AF, M, T, TSIZE, WORK, LWORK, INFO );


      // Generate the n-by-n matrix Q

      slaset('Full', N, N, ZERO, ONE, Q, N );
      srnamt = 'SGEMLQ';
      sgemlq('R', 'N', N, N, K, AF, M, T, TSIZE, Q, N, WORK, LWORK, INFO );

      // Copy R

      slaset('Full', M, N, ZERO, ZERO, LQ, L );
      slacpy('Lower', M, N, AF, M, LQ, L );

      // Compute |L - A*Q'| / |A| and store in RESULT(1)

      sgemm('N', 'T', M, N, N, -ONE, A, M, Q, N, ONE, LQ, L );
      ANORM = SLANGE( '1', M, N, A, M, RWORK );
      RESID = SLANGE( '1', M, N, LQ, L, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = RESID / (EPS*MAX(1,N)*ANORM);
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      slaset('Full', N, N, ZERO, ONE, LQ, L );
      ssyrk('U', 'C', N, N, -ONE, Q, N, ONE, LQ, L );
      RESID = SLANSY( '1', 'Upper', N, LQ, L, RWORK );
      RESULT( 2 ) = RESID / (EPS*MAX(1,N));

      // Generate random m-by-n matrix C and a copy CF

      for (J = 1; J <= M; J++) {
         slarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = SLANGE( '1', N, M, D, N, RWORK);
      slacpy('Full', N, M, D, N, DF, N );

      // Apply Q to C as Q*C

      sgemlq('L', 'N', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

      // Compute |Q*D - Q*D| / |D|

      sgemm('N', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N );
      RESID = SLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT( 3 ) = RESID / (EPS*MAX(1,N)*DNORM);
      } else {
         RESULT( 3 ) = ZERO;
      }

      // Copy D into DF again

      slacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as QT*D

      sgemlq('L', 'T', N, M, K, AF, M, T, TSIZE, DF, N, WORK, LWORK, INFO);

      // Compute |QT*D - QT*D| / |D|

      sgemm('T', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N );
      RESID = SLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT( 4 ) = RESID / (EPS*MAX(1,N)*DNORM);
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

      sgemlq('R', 'N', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

      // Compute |C*Q - C*Q| / |C|

      sgemm('N', 'N', M, N, N, -ONE, C, M, Q, N, ONE, CF, M );
      RESID = SLANGE( '1', N, M, DF, N, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 5 ) = RESID / (EPS*MAX(1,N)*CNORM);
      } else {
         RESULT( 5 ) = ZERO;
      }

      // Copy C into CF again

      slacpy('Full', M, N, C, M, CF, M );

      // Apply Q to D as D*QT

      sgemlq('R', 'T', M, N, K, AF, M, T, TSIZE, CF, M, WORK, LWORK, INFO);

      // Compute |C*QT - C*QT| / |C|

      sgemm('N', 'T', M, N, N, -ONE, C, M, Q, N, ONE, CF, M );
      RESID = SLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT( 6 ) = RESID / (EPS*MAX(1,N)*CNORM);
      } else {
         RESULT( 6 ) = ZERO;
      }

      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF);

      return;
      }
