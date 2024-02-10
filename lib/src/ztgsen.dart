      void ztgsen(IJOB, WANTQ, WANTZ, SELECT, N, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, ALPHA, BETA, final Matrix<double> Q, final int LDQ, final Matrix<double> Z, final int LDZ, M, PL, PR, DIF, final Array<double> WORK, final int LWORK, final Array<int> IWORK, final int LIWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               WANTQ, WANTZ;
      int                IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK, M, N;
      double             PL, PR;
      bool               SELECT( * );
      int                IWORK( * );
      double             DIF( * );
      Complex         A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * );
      // ..

      int                IDIFJB;
      const              IDIFJB = 3 ;
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               LQUERY, SWAP, WANTD, WANTD1, WANTD2, WANTP;
      int                I, IERR, IJB, K, KASE, KS, LIWMIN, LWMIN, MN2, N1, N2;
      double             DSCALE, DSUM, RDSCAL, SAFMIN;
      Complex         TEMP1, TEMP2;
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACN2, ZLACPY, ZLASSQ, ZSCAL, ZTGEXC, ZTGSYL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCMPLX, DCONJG, MAX, SQRT
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH

      // Decode and test the input parameters

      INFO = 0;
      LQUERY = ( LWORK == -1 || LIWORK == -1 );

      if ( IJOB < 0 || IJOB > 5 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( LDQ < 1 || ( WANTQ && LDQ < N ) ) {
         INFO = -13;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -15;
      }

      if ( INFO != 0 ) {
         xerbla('ZTGSEN', -INFO );
         return;
      }

      IERR = 0;

      WANTP = IJOB == 1 || IJOB >= 4;
      WANTD1 = IJOB == 2 || IJOB == 4;
      WANTD2 = IJOB == 3 || IJOB == 5;
      WANTD = WANTD1 || WANTD2;

      // Set M to the dimension of the specified pair of deflating
      // subspaces.

      M = 0;
      if ( !LQUERY || IJOB != 0 ) {
      for (K = 1; K <= N; K++) { // 10
         ALPHA[K] = A( K, K );
         BETA[K] = B( K, K );
         if ( K < N ) {
            if( SELECT( K ) ) M = M + 1;
         } else {
            if( SELECT( N ) ) M = M + 1;
         }
      } // 10
      }

      if ( IJOB == 1 || IJOB == 2 || IJOB == 4 ) {
         LWMIN = max( 1, 2*M*( N-M ) );
         LIWMIN = max( 1, N+2 );
      } else if ( IJOB == 3 || IJOB == 5 ) {
         LWMIN = max( 1, 4*M*( N-M ) );
         LIWMIN = max( 1, 2*M*( N-M ), N+2 );
      } else {
         LWMIN = 1;
         LIWMIN = 1;
      }

      WORK[1] = LWMIN;
      IWORK[1] = LIWMIN;

      if ( LWORK < LWMIN && !LQUERY ) {
         INFO = -21;
      } else if ( LIWORK < LIWMIN && !LQUERY ) {
         INFO = -23;
      }

      if ( INFO != 0 ) {
         xerbla('ZTGSEN', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible.

      if ( M == N || M == 0 ) {
         if ( WANTP ) {
            PL = ONE;
            PR = ONE;
         }
         if ( WANTD ) {
            DSCALE = ZERO;
            DSUM = ONE;
            for (I = 1; I <= N; I++) { // 20
               zlassq(N, A( 1, I ), 1, DSCALE, DSUM );
               zlassq(N, B( 1, I ), 1, DSCALE, DSUM );
            } // 20
            DIF[1] = DSCALE*sqrt( DSUM );
            DIF[2] = DIF( 1 );
         }
         GO TO 70;
      }

      // Get machine constant

      SAFMIN = dlamch( 'S' );

      // Collect the selected blocks at the top-left corner of (A, B).

      KS = 0;
      for (K = 1; K <= N; K++) { // 30
         SWAP = SELECT( K );
         if ( SWAP ) {
            KS = KS + 1;

            // Swap the K-th block to position KS. Compute unitary Q
            // and Z that will swap adjacent diagonal blocks in (A, B).

            if (K != KS) ztgexc( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, K, KS, IERR );

            if ( IERR > 0 ) {

               // Swap is rejected: exit.

               INFO = 1;
               if ( WANTP ) {
                  PL = ZERO;
                  PR = ZERO;
               }
               if ( WANTD ) {
                  DIF[1] = ZERO;
                  DIF[2] = ZERO;
               }
               GO TO 70;
            }
         }
      } // 30
      if ( WANTP ) {

         // Solve generalized Sylvester equation for R and L:
         //            A11 * R - L * A22 = A12
         //            B11 * R - L * B22 = B12

         N1 = M;
         N2 = N - M;
         I = N1 + 1;
         zlacpy('Full', N1, N2, A( 1, I ), LDA, WORK, N1 );
         zlacpy('Full', N1, N2, B( 1, I ), LDB, WORK( N1*N2+1 ), N1 );
         IJB = 0;
         ztgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );

         // Estimate the reciprocal of norms of "projections" onto
         // left and right eigenspaces

         RDSCAL = ZERO;
         DSUM = ONE;
         zlassq(N1*N2, WORK, 1, RDSCAL, DSUM );
         PL = RDSCAL*sqrt( DSUM );
         if ( PL == ZERO ) {
            PL = ONE;
         } else {
            PL = DSCALE / ( sqrt( DSCALE*DSCALE / PL+PL )*sqrt( PL ) );
         }
         RDSCAL = ZERO;
         DSUM = ONE;
         zlassq(N1*N2, WORK( N1*N2+1 ), 1, RDSCAL, DSUM );
         PR = RDSCAL*sqrt( DSUM );
         if ( PR == ZERO ) {
            PR = ONE;
         } else {
            PR = DSCALE / ( sqrt( DSCALE*DSCALE / PR+PR )*sqrt( PR ) );
         }
      }
      if ( WANTD ) {

         // Compute estimates Difu and Difl.

         if ( WANTD1 ) {
            N1 = M;
            N2 = N - M;
            I = N1 + 1;
            IJB = IDIFJB;

            // Frobenius norm-based Difu estimate.

            ztgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );

            // Frobenius norm-based Difl estimate.

            ztgsyl('N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
         } else {

            // Compute 1-norm-based estimates of Difu and Difl using
            // reversed communication with ZLACN2. In each step a
            // generalized Sylvester equation or a transposed variant
            // is solved.

            KASE = 0;
            N1 = M;
            N2 = N - M;
            I = N1 + 1;
            IJB = 0;
            MN2 = 2*N1*N2;

            // 1-norm-based estimate of Difu.

            } // 40
            zlacn2(MN2, WORK( MN2+1 ), WORK, DIF( 1 ), KASE, ISAVE );
            if ( KASE != 0 ) {
               if ( KASE == 1 ) {

                  // Solve generalized Sylvester equation

                  ztgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               } else {

                  // Solve the transposed variant.

                  ztgsyl('C', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               }
               GO TO 40;
            }
            DIF[1] = DSCALE / DIF( 1 );

            // 1-norm-based estimate of Difl.

            } // 50
            zlacn2(MN2, WORK( MN2+1 ), WORK, DIF( 2 ), KASE, ISAVE );
            if ( KASE != 0 ) {
               if ( KASE == 1 ) {

                  // Solve generalized Sylvester equation

                  ztgsyl('N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               } else {

                  // Solve the transposed variant.

                  ztgsyl('C', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               }
               GO TO 50;
            }
            DIF[2] = DSCALE / DIF( 2 );
         }
      }

      // If B(K,K) is complex, make it real and positive (normalization
      // of the generalized Schur form) and Store the generalized
      // eigenvalues of reordered pair (A, B)

      for (K = 1; K <= N; K++) { // 60
         DSCALE = ( B( K, K ) ).abs();
         if ( DSCALE > SAFMIN ) {
            TEMP1 = DCONJG( B( K, K ) / DSCALE );
            TEMP2 = B( K, K ) / DSCALE;
            B[K][K] = DSCALE;
            zscal(N-K, TEMP1, B( K, K+1 ), LDB );
            zscal(N-K+1, TEMP1, A( K, K ), LDA );
            if (WANTQ) zscal( N, TEMP2, Q( 1, K ), 1 );
         } else {
            B[K][K] = DCMPLX( ZERO, ZERO );
         }

         ALPHA[K] = A( K, K );
         BETA[K] = B( K, K );

      } // 60

      } // 70

      WORK[1] = LWMIN;
      IWORK[1] = LIWMIN;

      }
