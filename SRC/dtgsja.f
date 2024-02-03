      SUBROUTINE DTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, NCYCLE, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, NCYCLE, P;
      double             TOLA, TOLB;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), Q( LDQ, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                MAXIT;
      const              MAXIT = 40 ;
      double             ZERO, ONE, HUGENUM;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..

      bool               INITQ, INITU, INITV, UPPER, WANTQ, WANTU, WANTV;
      int                I, J, KCYCLE;
      double             A1, A2, A3, B1, B2, B3, CSQ, CSU, CSV, ERROR, GAMMA, RWK, SNQ, SNU, SNV, SSMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLAGS2, DLAPLL, DLARTG, DLASET, DROT, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, HUGE
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
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -10;
      } else if ( LDB < MAX( 1, P ) ) {
         INFO = -12;
      } else if ( LDU < 1 || ( WANTU && LDU < M ) ) {
         INFO = -18;
      } else if ( LDV < 1 || ( WANTV && LDV < P ) ) {
         INFO = -20;
      } else if ( LDQ < 1 || ( WANTQ && LDQ < N ) ) {
         INFO = -22;
      }
      if ( INFO != 0 ) {
         xerbla('DTGSJA', -INFO );
         RETURN;
      }

      // Initialize U, V and Q, if necessary

      if (INITU) CALL DLASET( 'Full', M, M, ZERO, ONE, U, LDU )       IF( INITV ) CALL DLASET( 'Full', P, P, ZERO, ONE, V, LDV )       IF( INITQ ) CALL DLASET( 'Full', N, N, ZERO, ONE, Q, LDQ );

      // Loop until convergence

      UPPER = false;
      for (KCYCLE = 1; KCYCLE <= MAXIT; KCYCLE++) { // 40

         UPPER = !UPPER;

         for (I = 1; I <= L - 1; I++) { // 20
            for (J = I + 1; J <= L; J++) { // 10

               A1 = ZERO;
               A2 = ZERO;
               A3 = ZERO;
               if (K+I <= M) A1 = A( K+I, N-L+I )                IF( K+J <= M ) A3 = A( K+J, N-L+J );

               B1 = B( I, N-L+I );
               B3 = B( J, N-L+J );

               if ( UPPER ) {
                  if (K+I <= M) A2 = A( K+I, N-L+J );
                  B2 = B( I, N-L+J );
               } else {
                  if (K+J <= M) A2 = A( K+J, N-L+I );
                  B2 = B( J, N-L+I );
               }

               dlags2(UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ );

               // Update (K+I)-th and (K+J)-th rows of matrix A: U**T *A

               if (K+J <= M) CALL DROT( L, A( K+J, N-L+1 ), LDA, A( K+I, N-L+1 ), LDA, CSU, SNU );

               // Update I-th and J-th rows of matrix B: V**T *B

               drot(L, B( J, N-L+1 ), LDB, B( I, N-L+1 ), LDB, CSV, SNV );

               // Update (N-L+I)-th and (N-L+J)-th columns of matrices
               // A and B: A*Q and B*Q

               drot(MIN( K+L, M ), A( 1, N-L+J ), 1, A( 1, N-L+I ), 1, CSQ, SNQ );

               drot(L, B( 1, N-L+J ), 1, B( 1, N-L+I ), 1, CSQ, SNQ );

               if ( UPPER ) {
                  if (K+I <= M) A( K+I, N-L+J ) = ZERO;
                  B( I, N-L+J ) = ZERO;
               } else {
                  if (K+J <= M) A( K+J, N-L+I ) = ZERO;
                  B( J, N-L+I ) = ZERO;
               }

               // Update orthogonal matrices U, V, Q, if desired.

               if (WANTU && K+J <= M) CALL DROT( M, U( 1, K+J ), 1, U( 1, K+I ), 1, CSU, SNU );

               if (WANTV) CALL DROT( P, V( 1, J ), 1, V( 1, I ), 1, CSV, SNV );

               if (WANTQ) CALL DROT( N, Q( 1, N-L+J ), 1, Q( 1, N-L+I ), 1, CSQ, SNQ );

            } // 10
         } // 20

         if ( !UPPER ) {

            // The matrices A13 and B13 were lower triangular at the start
            // of the cycle, and are now upper triangular.

            // Convergence test: test the parallelism of the corresponding
            // rows of A and B.

            ERROR = ZERO;
            DO 30 I = 1, MIN( L, M-K );
               dcopy(L-I+1, A( K+I, N-L+I ), LDA, WORK, 1 );
               dcopy(L-I+1, B( I, N-L+I ), LDB, WORK( L+1 ), 1 );
               dlapll(L-I+1, WORK, 1, WORK( L+1 ), 1, SSMIN );
               ERROR = MAX( ERROR, SSMIN );
            } // 30

            IF( ABS( ERROR ) <= MIN( TOLA, TOLB ) ) GO TO 50;
         }

         // End of cycle loop

      } // 40

      // The algorithm has not converged after MAXIT cycles.

      INFO = 1;
      GO TO 100;

      } // 50

      // If ERROR <= MIN(TOLA,TOLB), then the algorithm has converged.
      // Compute the generalized singular value pairs (ALPHA, BETA), and
      // set the triangular matrix R to array A.

      for (I = 1; I <= K; I++) { // 60
         ALPHA( I ) = ONE;
         BETA( I ) = ZERO;
      } // 60

      DO 70 I = 1, MIN( L, M-K );

         A1 = A( K+I, N-L+I );
         B1 = B( I, N-L+I );
         GAMMA = B1 / A1;

         if ( (GAMMA <= HUGENUM) && (GAMMA >= -HUGENUM) ) {

            // change sign if necessary

            if ( GAMMA < ZERO ) {
               dscal(L-I+1, -ONE, B( I, N-L+I ), LDB );
               if (WANTV) CALL DSCAL( P, -ONE, V( 1, I ), 1 );
            }

            dlartg(ABS( GAMMA ), ONE, BETA( K+I ), ALPHA( K+I ), RWK );

            if ( ALPHA( K+I ) >= BETA( K+I ) ) {
               dscal(L-I+1, ONE / ALPHA( K+I ), A( K+I, N-L+I ), LDA );
            } else {
               dscal(L-I+1, ONE / BETA( K+I ), B( I, N-L+I ), LDB );
               dcopy(L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ), LDA );
            }

         } else {

            ALPHA( K+I ) = ZERO;
            BETA( K+I ) = ONE;
            dcopy(L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ), LDA );

         }

      } // 70

      // Post-assignment

      for (I = M + 1; I <= K + L; I++) { // 80
         ALPHA( I ) = ZERO;
         BETA( I ) = ONE;
      } // 80

      if ( K+L < N ) {
         for (I = K + L + 1; I <= N; I++) { // 90
            ALPHA( I ) = ZERO;
            BETA( I ) = ZERO;
         } // 90
      }

      } // 100
      NCYCLE = KCYCLE;
      RETURN;

      // End of DTGSJA

      }
