      SUBROUTINE STGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               WANTQ, WANTZ;
      int                IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK, M, N;
      REAL               PL, PR
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      int                IWORK( * );
      REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), DIF( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                IDIFJB;
      const              IDIFJB = 3 ;
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, PAIR, SWAP, WANTD, WANTD1, WANTD2, WANTP;
      int                I, IERR, IJB, K, KASE, KK, KS, LIWMIN, LWMIN, MN2, N1, N2;
      REAL               DSCALE, DSUM, EPS, RDSCAL, SMLNUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SLACPY, SLAG2, SLASSQ, STGEXC, STGSYL, XERBLA
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SROUNDUP_LWORK
      // EXTERNAL SLAMCH, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      INFO = 0
      LQUERY = ( LWORK == -1 .OR. LIWORK == -1 )

      if ( IJOB.LT.0 .OR. IJOB.GT.5 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -9
      } else if ( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) {
         INFO = -14
      } else if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
         INFO = -16
      }

      if ( INFO.NE.0 ) {
         xerbla('STGSEN', -INFO );
         RETURN
      }

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      IERR = 0

      WANTP = IJOB == 1 .OR. IJOB.GE.4
      WANTD1 = IJOB == 2 .OR. IJOB == 4
      WANTD2 = IJOB == 3 .OR. IJOB == 5
      WANTD = WANTD1 .OR. WANTD2

      // Set M to the dimension of the specified pair of deflating
      // subspaces.

      M = 0
      PAIR = false;
      if ( .NOT.LQUERY .OR. IJOB.NE.0 ) {
      for (K = 1; K <= N; K++) { // 10
         if ( PAIR ) {
            PAIR = false;
         } else {
            if ( K.LT.N ) {
               if ( A( K+1, K ) == ZERO ) {
                  IF( SELECT( K ) ) M = M + 1
               } else {
                  PAIR = true;
                  IF( SELECT( K ) .OR. SELECT( K+1 ) ) M = M + 2
               }
            } else {
               IF( SELECT( N ) ) M = M + 1
            }
         }
      } // 10
      }

      if ( IJOB == 1 .OR. IJOB == 2 .OR. IJOB == 4 ) {
         LWMIN = MAX( 1, 4*N+16, 2*M*(N-M) )
         LIWMIN = MAX( 1, N+6 )
      } else if ( IJOB == 3 .OR. IJOB == 5 ) {
         LWMIN = MAX( 1, 4*N+16, 4*M*(N-M) )
         LIWMIN = MAX( 1, 2*M*(N-M), N+6 )
      } else {
         LWMIN = MAX( 1, 4*N+16 )
         LIWMIN = 1
      }

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN

      if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
         INFO = -22
      } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
         INFO = -24
      }

      if ( INFO.NE.0 ) {
         xerbla('STGSEN', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible.

      if ( M == N .OR. M == 0 ) {
         if ( WANTP ) {
            PL = ONE
            PR = ONE
         }
         if ( WANTD ) {
            DSCALE = ZERO
            DSUM = ONE
            for (I = 1; I <= N; I++) { // 20
               slassq(N, A( 1, I ), 1, DSCALE, DSUM );
               slassq(N, B( 1, I ), 1, DSCALE, DSUM );
            } // 20
            DIF( 1 ) = DSCALE*SQRT( DSUM )
            DIF( 2 ) = DIF( 1 )
         }
         GO TO 60
      }

      // Collect the selected blocks at the top-left corner of (A, B).

      KS = 0
      PAIR = false;
      for (K = 1; K <= N; K++) { // 30
         if ( PAIR ) {
            PAIR = false;
         } else {

            SWAP = SELECT( K )
            if ( K.LT.N ) {
               if ( A( K+1, K ).NE.ZERO ) {
                  PAIR = true;
                  SWAP = SWAP .OR. SELECT( K+1 )
               }
            }

            if ( SWAP ) {
               KS = KS + 1

               // Swap the K-th block to position KS.
               // Perform the reordering of diagonal blocks in (A, B)
               // by orthogonal transformation matrices and update
               // Q and Z accordingly (if requested):

               KK = K
               if (K.NE.KS) CALL STGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, KK, KS, WORK, LWORK, IERR );

               if ( IERR.GT.0 ) {

                  // Swap is rejected: exit.

                  INFO = 1
                  if ( WANTP ) {
                     PL = ZERO
                     PR = ZERO
                  }
                  if ( WANTD ) {
                     DIF( 1 ) = ZERO
                     DIF( 2 ) = ZERO
                  }
                  GO TO 60
               }

               if (PAIR) KS = KS + 1;
            }
         }
      } // 30
      if ( WANTP ) {

         // Solve generalized Sylvester equation for R and L
         // and compute PL and PR.

         N1 = M
         N2 = N - M
         I = N1 + 1
         IJB = 0
         slacpy('Full', N1, N2, A( 1, I ), LDA, WORK, N1 );
         slacpy('Full', N1, N2, B( 1, I ), LDB, WORK( N1*N2+1 ), N1 );
         stgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );

         // Estimate the reciprocal of norms of "projections" onto left
         // and right eigenspaces.

         RDSCAL = ZERO
         DSUM = ONE
         slassq(N1*N2, WORK, 1, RDSCAL, DSUM );
         PL = RDSCAL*SQRT( DSUM )
         if ( PL == ZERO ) {
            PL = ONE
         } else {
            PL = DSCALE / ( SQRT( DSCALE*DSCALE / PL+PL )*SQRT( PL ) )
         }
         RDSCAL = ZERO
         DSUM = ONE
         slassq(N1*N2, WORK( N1*N2+1 ), 1, RDSCAL, DSUM );
         PR = RDSCAL*SQRT( DSUM )
         if ( PR == ZERO ) {
            PR = ONE
         } else {
            PR = DSCALE / ( SQRT( DSCALE*DSCALE / PR+PR )*SQRT( PR ) )
         }
      }

      if ( WANTD ) {

         // Compute estimates of Difu and Difl.

         if ( WANTD1 ) {
            N1 = M
            N2 = N - M
            I = N1 + 1
            IJB = IDIFJB

            // Frobenius norm-based Difu-estimate.

            stgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR );

            // Frobenius norm-based Difl-estimate.

            stgsyl('N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR );
         } else {


            // Compute 1-norm-based estimates of Difu and Difl using
            // reversed communication with SLACN2. In each step a
            // generalized Sylvester equation or a transposed variant
            // is solved.

            KASE = 0
            N1 = M
            N2 = N - M
            I = N1 + 1
            IJB = 0
            MN2 = 2*N1*N2

            // 1-norm-based estimate of Difu.

            } // 40
            slacn2(MN2, WORK( MN2+1 ), WORK, IWORK, DIF( 1 ), KASE, ISAVE );
            if ( KASE.NE.0 ) {
               if ( KASE == 1 ) {

                  // Solve generalized Sylvester equation.

                  stgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               } else {

                  // Solve the transposed variant.

                  stgsyl('T', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               }
               GO TO 40
            }
            DIF( 1 ) = DSCALE / DIF( 1 )

            // 1-norm-based estimate of Difl.

            } // 50
            slacn2(MN2, WORK( MN2+1 ), WORK, IWORK, DIF( 2 ), KASE, ISAVE );
            if ( KASE.NE.0 ) {
               if ( KASE == 1 ) {

                  // Solve generalized Sylvester equation.

                  stgsyl('N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               } else {

                  // Solve the transposed variant.

                  stgsyl('T', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               }
               GO TO 50
            }
            DIF( 2 ) = DSCALE / DIF( 2 )

         }
      }

      } // 60

      // Compute generalized eigenvalues of reordered pair (A, B) and
      // normalize the generalized Schur form.

      PAIR = false;
      for (K = 1; K <= N; K++) { // 70
         if ( PAIR ) {
            PAIR = false;
         } else {

            if ( K.LT.N ) {
               if ( A( K+1, K ).NE.ZERO ) {
                  PAIR = true;
               }
            }

            if ( PAIR ) {

              // Compute the eigenvalue(s) at position K.

               WORK( 1 ) = A( K, K )
               WORK( 2 ) = A( K+1, K )
               WORK( 3 ) = A( K, K+1 )
               WORK( 4 ) = A( K+1, K+1 )
               WORK( 5 ) = B( K, K )
               WORK( 6 ) = B( K+1, K )
               WORK( 7 ) = B( K, K+1 )
               WORK( 8 ) = B( K+1, K+1 )
               slag2(WORK, 2, WORK( 5 ), 2, SMLNUM*EPS, BETA( K ), BETA( K+1 ), ALPHAR( K ), ALPHAR( K+1 ), ALPHAI( K ) );
               ALPHAI( K+1 ) = -ALPHAI( K )

            } else {

               if ( SIGN( ONE, B( K, K ) ).LT.ZERO ) {

                  // If B(K,K) is negative, make it positive

                  for (I = 1; I <= N; I++) { // 80
                     A( K, I ) = -A( K, I )
                     B( K, I ) = -B( K, I )
                     if (WANTQ) Q( I, K ) = -Q( I, K );
                  } // 80
               }

               ALPHAR( K ) = A( K, K )
               ALPHAI( K ) = ZERO
               BETA( K ) = B( K, K )

            }
         }
      } // 70

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of STGSEN

      }
