      void dtgsen(IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF, WORK, LWORK, IWORK, LIWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               WANTQ, WANTZ;
      int                IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK, M, N;
      double             PL, PR;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      int                IWORK( * );
      double             A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), DIF( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                IDIFJB;
      const              IDIFJB = 3 ;
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, PAIR, SWAP, WANTD, WANTD1, WANTD2, WANTP;
      int                I, IERR, IJB, K, KASE, KK, KS, LIWMIN, LWMIN, MN2, N1, N2;
      double             DSCALE, DSUM, EPS, RDSCAL, SMLNUM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DLACPY, DLAG2, DLASSQ, DTGEXC, DTGSYL, XERBLA
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN, SQRT
      // ..
      // .. Executable Statements ..

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
         INFO = -14;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -16;
      }

      if ( INFO != 0 ) {
         xerbla('DTGSEN', -INFO );
         return;
      }

      // Get machine constants

      EPS = DLAMCH( 'P' );
      SMLNUM = DLAMCH( 'S' ) / EPS;
      IERR = 0;

      WANTP = IJOB == 1 || IJOB >= 4;
      WANTD1 = IJOB == 2 || IJOB == 4;
      WANTD2 = IJOB == 3 || IJOB == 5;
      WANTD = WANTD1 || WANTD2;

      // Set M to the dimension of the specified pair of deflating
      // subspaces.

      M = 0;
      PAIR = false;
      if ( !LQUERY || IJOB != 0 ) {
      for (K = 1; K <= N; K++) { // 10
         if ( PAIR ) {
            PAIR = false;
         } else {
            if ( K < N ) {
               if ( A( K+1, K ) == ZERO ) {
                  if( SELECT( K ) ) M = M + 1;
               } else {
                  PAIR = true;
                  if( SELECT( K ) || SELECT( K+1 ) ) M = M + 2;
               }
            } else {
               if( SELECT( N ) ) M = M + 1;
            }
         }
      } // 10
      }

      if ( IJOB == 1 || IJOB == 2 || IJOB == 4 ) {
         LWMIN = max( 1, 4*N+16, 2*M*( N-M ) );
         LIWMIN = max( 1, N+6 );
      } else if ( IJOB == 3 || IJOB == 5 ) {
         LWMIN = max( 1, 4*N+16, 4*M*( N-M ) );
         LIWMIN = max( 1, 2*M*( N-M ), N+6 );
      } else {
         LWMIN = max( 1, 4*N+16 );
         LIWMIN = 1;
      }

      WORK( 1 ) = LWMIN;
      IWORK( 1 ) = LIWMIN;

      if ( LWORK < LWMIN && !LQUERY ) {
         INFO = -22;
      } else if ( LIWORK < LIWMIN && !LQUERY ) {
         INFO = -24;
      }

      if ( INFO != 0 ) {
         xerbla('DTGSEN', -INFO );
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
               dlassq(N, A( 1, I ), 1, DSCALE, DSUM );
               dlassq(N, B( 1, I ), 1, DSCALE, DSUM );
            } // 20
            DIF( 1 ) = DSCALE*sqrt( DSUM );
            DIF( 2 ) = DIF( 1 );
         }
         GO TO 60;
      }

      // Collect the selected blocks at the top-left corner of (A, B).

      KS = 0;
      PAIR = false;
      for (K = 1; K <= N; K++) { // 30
         if ( PAIR ) {
            PAIR = false;
         } else {

            SWAP = SELECT( K );
            if ( K < N ) {
               if ( A( K+1, K ) != ZERO ) {
                  PAIR = true;
                  SWAP = SWAP || SELECT( K+1 );
               }
            }

            if ( SWAP ) {
               KS = KS + 1;

               // Swap the K-th block to position KS.
               // Perform the reordering of diagonal blocks in (A, B)
               // by orthogonal transformation matrices and update
               // Q and Z accordingly (if requested):

               KK = K;
               if (K != KS) CALL DTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, KK, KS, WORK, LWORK, IERR );

               if ( IERR > 0 ) {

                  // Swap is rejected: exit.

                  INFO = 1;
                  if ( WANTP ) {
                     PL = ZERO;
                     PR = ZERO;
                  }
                  if ( WANTD ) {
                     DIF( 1 ) = ZERO;
                     DIF( 2 ) = ZERO;
                  }
                  GO TO 60;
               }

               if (PAIR) KS = KS + 1;
            }
         }
      } // 30
      if ( WANTP ) {

         // Solve generalized Sylvester equation for R and L
         // and compute PL and PR.

         N1 = M;
         N2 = N - M;
         I = N1 + 1;
         IJB = 0;
         dlacpy('Full', N1, N2, A( 1, I ), LDA, WORK, N1 );
         dlacpy('Full', N1, N2, B( 1, I ), LDB, WORK( N1*N2+1 ), N1 );
         dtgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );

         // Estimate the reciprocal of norms of "projections" onto left
         // and right eigenspaces.

         RDSCAL = ZERO;
         DSUM = ONE;
         dlassq(N1*N2, WORK, 1, RDSCAL, DSUM );
         PL = RDSCAL*sqrt( DSUM );
         if ( PL == ZERO ) {
            PL = ONE;
         } else {
            PL = DSCALE / ( sqrt( DSCALE*DSCALE / PL+PL )*sqrt( PL ) );
         }
         RDSCAL = ZERO;
         DSUM = ONE;
         dlassq(N1*N2, WORK( N1*N2+1 ), 1, RDSCAL, DSUM );
         PR = RDSCAL*sqrt( DSUM );
         if ( PR == ZERO ) {
            PR = ONE;
         } else {
            PR = DSCALE / ( sqrt( DSCALE*DSCALE / PR+PR )*sqrt( PR ) );
         }
      }

      if ( WANTD ) {

         // Compute estimates of Difu and Difl.

         if ( WANTD1 ) {
            N1 = M;
            N2 = N - M;
            I = N1 + 1;
            IJB = IDIFJB;

            // Frobenius norm-based Difu-estimate.

            dtgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR );

            // Frobenius norm-based Difl-estimate.

            dtgsyl('N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR );
         } else {


            // Compute 1-norm-based estimates of Difu and Difl using
            // reversed communication with DLACN2. In each step a
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
            dlacn2(MN2, WORK( MN2+1 ), WORK, IWORK, DIF( 1 ), KASE, ISAVE );
            if ( KASE != 0 ) {
               if ( KASE == 1 ) {

                  // Solve generalized Sylvester equation.

                  dtgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               } else {

                  // Solve the transposed variant.

                  dtgsyl('T', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               }
               GO TO 40;
            }
            DIF( 1 ) = DSCALE / DIF( 1 );

            // 1-norm-based estimate of Difl.

            } // 50
            dlacn2(MN2, WORK( MN2+1 ), WORK, IWORK, DIF( 2 ), KASE, ISAVE );
            if ( KASE != 0 ) {
               if ( KASE == 1 ) {

                  // Solve generalized Sylvester equation.

                  dtgsyl('N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               } else {

                  // Solve the transposed variant.

                  dtgsyl('T', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               }
               GO TO 50;
            }
            DIF( 2 ) = DSCALE / DIF( 2 );

         }
      }

      } // 60

      // Compute generalized eigenvalues of reordered pair (A, B) and
      // normalize the generalized Schur form.

      PAIR = false;
      for (K = 1; K <= N; K++) { // 80
         if ( PAIR ) {
            PAIR = false;
         } else {

            if ( K < N ) {
               if ( A( K+1, K ) != ZERO ) {
                  PAIR = true;
               }
            }

            if ( PAIR ) {

              // Compute the eigenvalue(s) at position K.

               WORK( 1 ) = A( K, K );
               WORK( 2 ) = A( K+1, K );
               WORK( 3 ) = A( K, K+1 );
               WORK( 4 ) = A( K+1, K+1 );
               WORK( 5 ) = B( K, K );
               WORK( 6 ) = B( K+1, K );
               WORK( 7 ) = B( K, K+1 );
               WORK( 8 ) = B( K+1, K+1 );
               dlag2(WORK, 2, WORK( 5 ), 2, SMLNUM*EPS, BETA( K ), BETA( K+1 ), ALPHAR( K ), ALPHAR( K+1 ), ALPHAI( K ) );
               ALPHAI( K+1 ) = -ALPHAI( K );

            } else {

               if ( SIGN( ONE, B( K, K ) ) < ZERO ) {

                  // If B(K,K) is negative, make it positive

                  for (I = 1; I <= N; I++) { // 70
                     A( K, I ) = -A( K, I );
                     B( K, I ) = -B( K, I );
                     if (WANTQ) Q( I, K ) = -Q( I, K );
                  } // 70
               }

               ALPHAR( K ) = A( K, K );
               ALPHAI( K ) = ZERO;
               BETA( K ) = B( K, K );

            }
         }
      } // 80

      WORK( 1 ) = LWMIN;
      IWORK( 1 ) = LIWMIN;

      return;

      // End of DTGSEN

      }
