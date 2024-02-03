      SUBROUTINE CTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB, ALPHA, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF, WORK, LWORK, IWORK, LIWORK, INFO )

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
      REAL               DIF( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                IDIFJB;
      const              IDIFJB = 3 ;
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, SWAP, WANTD, WANTD1, WANTD2, WANTP;
      int                I, IERR, IJB, K, KASE, KS, LIWMIN, LWMIN, MN2, N1, N2;
      REAL               DSCALE, DSUM, RDSCAL, SAFMIN
      COMPLEX            TEMP1, TEMP2
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      REAL               SROUNDUP_LWORK
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      REAL               SLAMCH
      // EXTERNAL CLACN2, CLACPY, CLASSQ, CSCAL, CTGEXC, CTGSYL, SLAMCH, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, CONJG, MAX, SQRT
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
      } else if ( LDQ.LT.1 .OR. ( WANTQ && LDQ.LT.N ) ) {
         INFO = -13
      } else if ( LDZ.LT.1 .OR. ( WANTZ && LDZ.LT.N ) ) {
         INFO = -15
      }

      if ( INFO != 0 ) {
         xerbla('CTGSEN', -INFO );
         RETURN
      }

      IERR = 0

      WANTP = IJOB == 1 .OR. IJOB.GE.4
      WANTD1 = IJOB == 2 .OR. IJOB == 4
      WANTD2 = IJOB == 3 .OR. IJOB == 5
      WANTD = WANTD1 .OR. WANTD2

      // Set M to the dimension of the specified pair of deflating
      // subspaces.

      M = 0
      if ( .NOT.LQUERY .OR. IJOB != 0 ) {
      for (K = 1; K <= N; K++) { // 10
         ALPHA( K ) = A( K, K )
         BETA( K ) = B( K, K )
         if ( K.LT.N ) {
            IF( SELECT( K ) ) M = M + 1
         } else {
            IF( SELECT( N ) ) M = M + 1
         }
      } // 10
      }

      if ( IJOB == 1 .OR. IJOB == 2 .OR. IJOB == 4 ) {
         LWMIN = MAX( 1, 2*M*(N-M) )
         LIWMIN = MAX( 1, N+2 )
      } else if ( IJOB == 3 .OR. IJOB == 5 ) {
         LWMIN = MAX( 1, 4*M*(N-M) )
         LIWMIN = MAX( 1, 2*M*(N-M), N+2 )
      } else {
         LWMIN = 1
         LIWMIN = 1
      }

      WORK( 1 ) =  SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN

      if ( LWORK.LT.LWMIN && .NOT.LQUERY ) {
         INFO = -21
      } else if ( LIWORK.LT.LIWMIN && .NOT.LQUERY ) {
         INFO = -23
      }

      if ( INFO != 0 ) {
         xerbla('CTGSEN', -INFO );
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
               classq(N, A( 1, I ), 1, DSCALE, DSUM );
               classq(N, B( 1, I ), 1, DSCALE, DSUM );
            } // 20
            DIF( 1 ) = DSCALE*SQRT( DSUM )
            DIF( 2 ) = DIF( 1 )
         }
         GO TO 70
      }

      // Get machine constant

      SAFMIN = SLAMCH( 'S' )

      // Collect the selected blocks at the top-left corner of (A, B).

      KS = 0
      for (K = 1; K <= N; K++) { // 30
         SWAP = SELECT( K )
         if ( SWAP ) {
            KS = KS + 1

            // Swap the K-th block to position KS. Compute unitary Q
            // and Z that will swap adjacent diagonal blocks in (A, B).

            if (K != KS) CALL CTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, K, KS, IERR );

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
               GO TO 70
            }
         }
      } // 30
      if ( WANTP ) {

         // Solve generalized Sylvester equation for R and L:
                    // A11 * R - L * A22 = A12
                    // B11 * R - L * B22 = B12

         N1 = M
         N2 = N - M
         I = N1 + 1
         clacpy('Full', N1, N2, A( 1, I ), LDA, WORK, N1 );
         clacpy('Full', N1, N2, B( 1, I ), LDB, WORK( N1*N2+1 ), N1 );
         IJB = 0
         ctgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );

         // Estimate the reciprocal of norms of "projections" onto
         // left and right eigenspaces

         RDSCAL = ZERO
         DSUM = ONE
         classq(N1*N2, WORK, 1, RDSCAL, DSUM );
         PL = RDSCAL*SQRT( DSUM )
         if ( PL == ZERO ) {
            PL = ONE
         } else {
            PL = DSCALE / ( SQRT( DSCALE*DSCALE / PL+PL )*SQRT( PL ) )
         }
         RDSCAL = ZERO
         DSUM = ONE
         classq(N1*N2, WORK( N1*N2+1 ), 1, RDSCAL, DSUM );
         PR = RDSCAL*SQRT( DSUM )
         if ( PR == ZERO ) {
            PR = ONE
         } else {
            PR = DSCALE / ( SQRT( DSCALE*DSCALE / PR+PR )*SQRT( PR ) )
         }
      }
      if ( WANTD ) {

         // Compute estimates Difu and Difl.

         if ( WANTD1 ) {
            N1 = M
            N2 = N - M
            I = N1 + 1
            IJB = IDIFJB

            // Frobenius norm-based Difu estimate.

            ctgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );

            // Frobenius norm-based Difl estimate.

            ctgsyl('N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
         } else {

            // Compute 1-norm-based estimates of Difu and Difl using
            // reversed communication with CLACN2. In each step a
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
            clacn2(MN2, WORK( MN2+1 ), WORK, DIF( 1 ), KASE, ISAVE );
            if ( KASE != 0 ) {
               if ( KASE == 1 ) {

                  // Solve generalized Sylvester equation

                  ctgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               } else {

                  // Solve the transposed variant.

                  ctgsyl('C', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               }
               GO TO 40
            }
            DIF( 1 ) = DSCALE / DIF( 1 )

            // 1-norm-based estimate of Difl.

            } // 50
            clacn2(MN2, WORK( MN2+1 ), WORK, DIF( 2 ), KASE, ISAVE );
            if ( KASE != 0 ) {
               if ( KASE == 1 ) {

                  // Solve generalized Sylvester equation

                  ctgsyl('N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               } else {

                  // Solve the transposed variant.

                  ctgsyl('C', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               }
               GO TO 50
            }
            DIF( 2 ) = DSCALE / DIF( 2 )
         }
      }

      // If B(K,K) is complex, make it real and positive (normalization
      // of the generalized Schur form) and Store the generalized
      // eigenvalues of reordered pair (A, B)

      for (K = 1; K <= N; K++) { // 60
         DSCALE = ABS( B( K, K ) )
         if ( DSCALE.GT.SAFMIN ) {
            TEMP1 = CONJG( B( K, K ) / DSCALE )
            TEMP2 = B( K, K ) / DSCALE
            B( K, K ) = DSCALE
            cscal(N-K, TEMP1, B( K, K+1 ), LDB );
            cscal(N-K+1, TEMP1, A( K, K ), LDA );
            if (WANTQ) CALL CSCAL( N, TEMP2, Q( 1, K ), 1 );
         } else {
            B( K, K ) = CMPLX( ZERO, ZERO )
         }

         ALPHA( K ) = A( K, K )
         BETA( K ) = B( K, K )

      } // 60

      } // 70

      WORK( 1 ) =  SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of CTGSEN

      }
