      SUBROUTINE ZTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB, ALPHA, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               WANTQ, WANTZ;
      int                IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK, M, N;
      double             PL, PR;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      int                IWORK( * );
      double             DIF( * );
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                IDIFJB;
      const              IDIFJB = 3 ;
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, SWAP, WANTD, WANTD1, WANTD2, WANTP;
      int                I, IERR, IJB, K, KASE, KS, LIWMIN, LWMIN, MN2, N1, N2;
      double             DSCALE, DSUM, RDSCAL, SAFMIN;
      COMPLEX*16         TEMP1, TEMP2
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACN2, ZLACPY, ZLASSQ, ZSCAL, ZTGEXC, ZTGSYL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCMPLX, DCONJG, MAX, SQRT
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

      if ( IJOB.LT.0 .OR. IJOB.GT.5 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -9
      } else if ( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) {
         INFO = -13
      } else if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
         INFO = -15
      }

      if ( INFO.NE.0 ) {
         xerbla('ZTGSEN', -INFO );
         RETURN
      }

      IERR = 0

      WANTP = IJOB.EQ.1 .OR. IJOB.GE.4
      WANTD1 = IJOB.EQ.2 .OR. IJOB.EQ.4
      WANTD2 = IJOB.EQ.3 .OR. IJOB.EQ.5
      WANTD = WANTD1 .OR. WANTD2

      // Set M to the dimension of the specified pair of deflating
      // subspaces.

      M = 0
      if ( .NOT.LQUERY .OR. IJOB.NE.0 ) {
      for (K = 1; K <= N; K++) { // 10
         ALPHA( K ) = A( K, K )
         BETA( K ) = B( K, K )
         if ( K.LT.N ) {
            IF( SELECT( K ) ) M = M + 1
         } else {
            IF( SELECT( N ) ) M = M + 1
         }
   10 CONTINUE
      }

      if ( IJOB.EQ.1 .OR. IJOB.EQ.2 .OR. IJOB.EQ.4 ) {
         LWMIN = MAX( 1, 2*M*( N-M ) )
         LIWMIN = MAX( 1, N+2 )
      } else if ( IJOB.EQ.3 .OR. IJOB.EQ.5 ) {
         LWMIN = MAX( 1, 4*M*( N-M ) )
         LIWMIN = MAX( 1, 2*M*( N-M ), N+2 )
      } else {
         LWMIN = 1
         LIWMIN = 1
      }

      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN

      if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
         INFO = -21
      } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
         INFO = -23
      }

      if ( INFO.NE.0 ) {
         xerbla('ZTGSEN', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible.

      if ( M.EQ.N .OR. M.EQ.0 ) {
         if ( WANTP ) {
            PL = ONE
            PR = ONE
         }
         if ( WANTD ) {
            DSCALE = ZERO
            DSUM = ONE
            for (I = 1; I <= N; I++) { // 20
               zlassq(N, A( 1, I ), 1, DSCALE, DSUM );
               zlassq(N, B( 1, I ), 1, DSCALE, DSUM );
   20       CONTINUE
            DIF( 1 ) = DSCALE*SQRT( DSUM )
            DIF( 2 ) = DIF( 1 )
         }
         GO TO 70
      }

      // Get machine constant

      SAFMIN = DLAMCH( 'S' )

      // Collect the selected blocks at the top-left corner of (A, B).

      KS = 0
      for (K = 1; K <= N; K++) { // 30
         SWAP = SELECT( K )
         if ( SWAP ) {
            KS = KS + 1

            // Swap the K-th block to position KS. Compute unitary Q
            // and Z that will swap adjacent diagonal blocks in (A, B).

            IF( K.NE.KS ) CALL ZTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, K, KS, IERR )

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
   30 CONTINUE
      if ( WANTP ) {

         // Solve generalized Sylvester equation for R and L:
                    // A11 * R - L * A22 = A12
                    // B11 * R - L * B22 = B12

         N1 = M
         N2 = N - M
         I = N1 + 1
         zlacpy('Full', N1, N2, A( 1, I ), LDA, WORK, N1 );
         zlacpy('Full', N1, N2, B( 1, I ), LDB, WORK( N1*N2+1 ), N1 );
         IJB = 0
         ztgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );

         // Estimate the reciprocal of norms of "projections" onto
         // left and right eigenspaces

         RDSCAL = ZERO
         DSUM = ONE
         zlassq(N1*N2, WORK, 1, RDSCAL, DSUM );
         PL = RDSCAL*SQRT( DSUM )
         if ( PL.EQ.ZERO ) {
            PL = ONE
         } else {
            PL = DSCALE / ( SQRT( DSCALE*DSCALE / PL+PL )*SQRT( PL ) )
         }
         RDSCAL = ZERO
         DSUM = ONE
         zlassq(N1*N2, WORK( N1*N2+1 ), 1, RDSCAL, DSUM );
         PR = RDSCAL*SQRT( DSUM )
         if ( PR.EQ.ZERO ) {
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

            ztgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );

            // Frobenius norm-based Difl estimate.

            ztgsyl('N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
         } else {

            // Compute 1-norm-based estimates of Difu and Difl using
            // reversed communication with ZLACN2. In each step a
            // generalized Sylvester equation or a transposed variant
            // is solved.

            KASE = 0
            N1 = M
            N2 = N - M
            I = N1 + 1
            IJB = 0
            MN2 = 2*N1*N2

            // 1-norm-based estimate of Difu.

   40       CONTINUE
            zlacn2(MN2, WORK( MN2+1 ), WORK, DIF( 1 ), KASE, ISAVE );
            if ( KASE.NE.0 ) {
               if ( KASE.EQ.1 ) {

                  // Solve generalized Sylvester equation

                  ztgsyl('N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               } else {

                  // Solve the transposed variant.

                  ztgsyl('C', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               }
               GO TO 40
            }
            DIF( 1 ) = DSCALE / DIF( 1 )

            // 1-norm-based estimate of Difl.

   50       CONTINUE
            zlacn2(MN2, WORK( MN2+1 ), WORK, DIF( 2 ), KASE, ISAVE );
            if ( KASE.NE.0 ) {
               if ( KASE.EQ.1 ) {

                  // Solve generalized Sylvester equation

                  ztgsyl('N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
               } else {

                  // Solve the transposed variant.

                  ztgsyl('C', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR );
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
            TEMP1 = DCONJG( B( K, K ) / DSCALE )
            TEMP2 = B( K, K ) / DSCALE
            B( K, K ) = DSCALE
            zscal(N-K, TEMP1, B( K, K+1 ), LDB );
            zscal(N-K+1, TEMP1, A( K, K ), LDA );
            IF( WANTQ ) CALL ZSCAL( N, TEMP2, Q( 1, K ), 1 )
         } else {
            B( K, K ) = DCMPLX( ZERO, ZERO )
         }

         ALPHA( K ) = A( K, K )
         BETA( K ) = B( K, K )

   60 CONTINUE

   70 CONTINUE

      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of ZTGSEN

      }
