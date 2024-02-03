      SUBROUTINE CGGSVD( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, RWORK, IWORK, INFO );

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               ALPHA( * ), BETA( * ), RWORK( * );
      COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               WANTQ, WANTU, WANTV;
      int                I, IBND, ISUB, J, NCYCLE;
      REAL               ANORM, BNORM, SMAX, TEMP, TOLA, TOLB, ULP, UNFL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, SLAMCH;
      // EXTERNAL LSAME, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGGSVP, CTGSJA, SCOPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      WANTU = LSAME( JOBU, 'U' );
      WANTV = LSAME( JOBV, 'V' );
      WANTQ = LSAME( JOBQ, 'Q' );

      INFO = 0;
      if ( !( WANTU || LSAME( JOBU, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( WANTV || LSAME( JOBV, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( WANTQ || LSAME( JOBQ, 'N' ) ) ) {
         INFO = -3;
      } else if ( M < 0 ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( P < 0 ) {
         INFO = -6;
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -10;
      } else if ( LDB < MAX( 1, P ) ) {
         INFO = -12;
      } else if ( LDU < 1 || ( WANTU && LDU < M ) ) {
         INFO = -16;
      } else if ( LDV < 1 || ( WANTV && LDV < P ) ) {
         INFO = -18;
      } else if ( LDQ < 1 || ( WANTQ && LDQ < N ) ) {
         INFO = -20;
      }
      if ( INFO != 0 ) {
         xerbla('CGGSVD', -INFO );
         RETURN;
      }

      // Compute the Frobenius norm of matrices A and B

      ANORM = CLANGE( '1', M, N, A, LDA, RWORK );
      BNORM = CLANGE( '1', P, N, B, LDB, RWORK );

      // Get machine precision and set up threshold for determining
      // the effective numerical rank of the matrices A and B.

      ULP = SLAMCH( 'Precision' );
      UNFL = SLAMCH( 'Safe Minimum' );
      TOLA = MAX( M, N )*MAX( ANORM, UNFL )*ULP;
      TOLB = MAX( P, N )*MAX( BNORM, UNFL )*ULP;

      cggsvp(JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, RWORK, WORK, WORK( N+1 ), INFO );

      // Compute the GSVD of two upper "triangular" matrices

      ctgsja(JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, NCYCLE, INFO );

      // Sort the singular values and store the pivot indices in IWORK
      // Copy ALPHA to RWORK, then sort ALPHA in RWORK

      scopy(N, ALPHA, 1, RWORK, 1 );
      IBND = MIN( L, M-K );
      for (I = 1; I <= IBND; I++) { // 20

         // Scan for largest ALPHA(K+I)

         ISUB = I;
         SMAX = RWORK( K+I );
         for (J = I + 1; J <= IBND; J++) { // 10
            TEMP = RWORK( K+J );
            if ( TEMP > SMAX ) {
               ISUB = J;
               SMAX = TEMP;
            }
         } // 10
         if ( ISUB != I ) {
            RWORK( K+ISUB ) = RWORK( K+I );
            RWORK( K+I ) = SMAX;
            IWORK( K+I ) = K + ISUB;
         } else {
            IWORK( K+I ) = K + I;
         }
      } // 20

      RETURN;

      // End of CGGSVD

      }
