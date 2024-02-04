      void ctgsja(JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, NCYCLE, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, NCYCLE, P;
      REAL               TOLA, TOLB;
      // ..
      // .. Array Arguments ..
      REAL               ALPHA( * ), BETA( * );
      COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                MAXIT;
      const              MAXIT = 40 ;
      REAL               ZERO, ONE, HUGENUM;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..

      bool               INITQ, INITU, INITV, UPPER, WANTQ, WANTU, WANTV;
      int                I, J, KCYCLE;
      REAL               A1, A3, B1, B3, CSQ, CSU, CSV, ERROR, GAMMA, RWK, SSMIN;
      COMPLEX            A2, B2, SNQ, SNU, SNV;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLAGS2, CLAPLL, CLASET, CROT, CSSCAL, SLARTG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, MIN, REAL, HUGE
      const              HUGENUM = HUGE(ZERO) ;
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      INITU = LSAME( JOBU, 'I' );
      WANTU = INITU || LSAME( JOBU, 'U' );

      INITV = LSAME( JOBV, 'I' );
      WANTV = INITV || LSAME( JOBV, 'V' );

      INITQ = LSAME( JOBQ, 'I' );
      WANTQ = INITQ || LSAME( JOBQ, 'Q' );

      INFO = 0;
      if ( !( INITU || WANTU || LSAME( JOBU, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( INITV || WANTV || LSAME( JOBV, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( INITQ || WANTQ || LSAME( JOBQ, 'N' ) ) ) {
         INFO = -3;
      } else if ( M < 0 ) {
         INFO = -4;
      } else if ( P < 0 ) {
         INFO = -5;
      } else if ( N < 0 ) {
         INFO = -6;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -10;
      } else if ( LDB < max( 1, P ) ) {
         INFO = -12;
      } else if ( LDU < 1 || ( WANTU && LDU < M ) ) {
         INFO = -18;
      } else if ( LDV < 1 || ( WANTV && LDV < P ) ) {
         INFO = -20;
      } else if ( LDQ < 1 || ( WANTQ && LDQ < N ) ) {
         INFO = -22;
      }
      if ( INFO != 0 ) {
         xerbla('CTGSJA', -INFO );
         return;
      }

      // Initialize U, V and Q, if necessary

      if (INITU) claset( 'Full', M, M, CZERO, CONE, U, LDU );
      if( INITV ) claset( 'Full', P, P, CZERO, CONE, V, LDV );
      IF( INITQ ) claset( 'Full', N, N, CZERO, CONE, Q, LDQ );

      // Loop until convergence

      UPPER = false;
      for (KCYCLE = 1; KCYCLE <= MAXIT; KCYCLE++) { // 40

         UPPER = !UPPER;

         for (I = 1; I <= L - 1; I++) { // 20
            for (J = I + 1; J <= L; J++) { // 10

               A1 = ZERO;
               A2 = CZERO;
               A3 = ZERO;
               if (K+I <= M) A1 = REAL( A( K+I, N-L+I ) );
               IF( K+J <= M ) A3 = REAL( A( K+J, N-L+J ) );

               B1 = REAL( B( I, N-L+I ) );
               B3 = REAL( B( J, N-L+J ) );

               if ( UPPER ) {
                  if (K+I <= M) A2 = A( K+I, N-L+J );
                  B2 = B( I, N-L+J );
               } else {
                  if (K+J <= M) A2 = A( K+J, N-L+I );
                  B2 = B( J, N-L+I );
               }

               clags2(UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ );

               // Update (K+I)-th and (K+J)-th rows of matrix A: U**H *A

               if (K+J <= M) crot( L, A( K+J, N-L+1 ), LDA, A( K+I, N-L+1 ), LDA, CSU, CONJG( SNU ) );

               // Update I-th and J-th rows of matrix B: V**H *B

               crot(L, B( J, N-L+1 ), LDB, B( I, N-L+1 ), LDB, CSV, CONJG( SNV ) );

               // Update (N-L+I)-th and (N-L+J)-th columns of matrices
               // A and B: A*Q and B*Q

               crot(min( K+L, M ), A( 1, N-L+J ), 1, A( 1, N-L+I ), 1, CSQ, SNQ );

               crot(L, B( 1, N-L+J ), 1, B( 1, N-L+I ), 1, CSQ, SNQ );

               if ( UPPER ) {
                  if (K+I <= M) A( K+I, N-L+J ) = CZERO;
                  B[I, N-L+J] = CZERO;
               } else {
                  if (K+J <= M) A( K+J, N-L+I ) = CZERO;
                  B[J, N-L+I] = CZERO;
               }

               // Ensure that the diagonal elements of A and B are real.

               if (K+I <= M) A( K+I, N-L+I ) = REAL( A( K+I, N-L+I ) );
               IF( K+J <= M ) A( K+J, N-L+J ) = REAL( A( K+J, N-L+J ) );
               B[I, N-L+I] = REAL( B( I, N-L+I ) );
               B[J, N-L+J] = REAL( B( J, N-L+J ) );

               // Update unitary matrices U, V, Q, if desired.

               if (WANTU && K+J <= M) crot( M, U( 1, K+J ), 1, U( 1, K+I ), 1, CSU, SNU );

               if (WANTV) crot( P, V( 1, J ), 1, V( 1, I ), 1, CSV, SNV );

               if (WANTQ) crot( N, Q( 1, N-L+J ), 1, Q( 1, N-L+I ), 1, CSQ, SNQ );

            } // 10
         } // 20

         if ( !UPPER ) {

            // The matrices A13 and B13 were lower triangular at the start
            // of the cycle, and are now upper triangular.

            // Convergence test: test the parallelism of the corresponding
            // rows of A and B.

            ERROR = ZERO;
            for (I = 1; I <= min( L, M-K ); I++) { // 30
               ccopy(L-I+1, A( K+I, N-L+I ), LDA, WORK, 1 );
               ccopy(L-I+1, B( I, N-L+I ), LDB, WORK( L+1 ), 1 );
               clapll(L-I+1, WORK, 1, WORK( L+1 ), 1, SSMIN );
               ERROR = max( ERROR, SSMIN );
            } // 30

            if( ( ERROR ).abs() <= min( TOLA, TOLB ) ) GO TO 50;
         }

         // End of cycle loop

      } // 40

      // The algorithm has not converged after MAXIT cycles.

      INFO = 1;
      GO TO 100;

      } // 50

      // If ERROR <= min(TOLA,TOLB), then the algorithm has converged.
      // Compute the generalized singular value pairs (ALPHA, BETA), and
      // set the triangular matrix R to array A.

      for (I = 1; I <= K; I++) { // 60
         ALPHA[I] = ONE;
         BETA[I] = ZERO;
      } // 60

      for (I = 1; I <= min( L, M-K ); I++) { // 70

         A1 = REAL( A( K+I, N-L+I ) );
         B1 = REAL( B( I, N-L+I ) );
         GAMMA = B1 / A1;

         if ( (GAMMA <= HUGENUM) && (GAMMA >= -HUGENUM) ) {

            if ( GAMMA < ZERO ) {
               csscal(L-I+1, -ONE, B( I, N-L+I ), LDB );
               if (WANTV) csscal( P, -ONE, V( 1, I ), 1 );
            }

            slartg(( GAMMA ).abs(), ONE, BETA( K+I ), ALPHA( K+I ), RWK );

            if ( ALPHA( K+I ) >= BETA( K+I ) ) {
               csscal(L-I+1, ONE / ALPHA( K+I ), A( K+I, N-L+I ), LDA );
            } else {
               csscal(L-I+1, ONE / BETA( K+I ), B( I, N-L+I ), LDB );
               ccopy(L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ), LDA );
            }

         } else {
            ALPHA[K+I] = ZERO;
            BETA[K+I] = ONE;
            ccopy(L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ), LDA );
         }
      } // 70

      // Post-assignment

      for (I = M + 1; I <= K + L; I++) { // 80
         ALPHA[I] = ZERO;
         BETA[I] = ONE;
      } // 80

      if ( K+L < N ) {
         for (I = K + L + 1; I <= N; I++) { // 90
            ALPHA[I] = ZERO;
            BETA[I] = ZERO;
         } // 90
      }

      } // 100
      NCYCLE = KCYCLE;

      return;
      }
