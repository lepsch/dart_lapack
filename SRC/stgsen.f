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
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

      IF( IJOB.LT.0 .OR. IJOB.GT.5 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -16
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STGSEN', -INFO )
         RETURN
      END IF

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      IERR = 0

      WANTP = IJOB.EQ.1 .OR. IJOB.GE.4
      WANTD1 = IJOB.EQ.2 .OR. IJOB.EQ.4
      WANTD2 = IJOB.EQ.3 .OR. IJOB.EQ.5
      WANTD = WANTD1 .OR. WANTD2

      // Set M to the dimension of the specified pair of deflating
      // subspaces.

      M = 0
      PAIR = .FALSE.
      IF( .NOT.LQUERY .OR. IJOB.NE.0 ) THEN
      DO 10 K = 1, N
         IF( PAIR ) THEN
            PAIR = .FALSE.
         } else {
            IF( K.LT.N ) THEN
               IF( A( K+1, K ).EQ.ZERO ) THEN
                  IF( SELECT( K ) ) M = M + 1
               } else {
                  PAIR = .TRUE.
                  IF( SELECT( K ) .OR. SELECT( K+1 ) ) M = M + 2
               END IF
            } else {
               IF( SELECT( N ) ) M = M + 1
            END IF
         END IF
   10 CONTINUE
      END IF

      IF( IJOB.EQ.1 .OR. IJOB.EQ.2 .OR. IJOB.EQ.4 ) THEN
         LWMIN = MAX( 1, 4*N+16, 2*M*(N-M) )
         LIWMIN = MAX( 1, N+6 )
      ELSE IF( IJOB.EQ.3 .OR. IJOB.EQ.5 ) THEN
         LWMIN = MAX( 1, 4*N+16, 4*M*(N-M) )
         LIWMIN = MAX( 1, 2*M*(N-M), N+6 )
      } else {
         LWMIN = MAX( 1, 4*N+16 )
         LIWMIN = 1
      END IF

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN

      IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -22
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -24
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STGSEN', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible.

      IF( M.EQ.N .OR. M.EQ.0 ) THEN
         IF( WANTP ) THEN
            PL = ONE
            PR = ONE
         END IF
         IF( WANTD ) THEN
            DSCALE = ZERO
            DSUM = ONE
            DO 20 I = 1, N
               CALL SLASSQ( N, A( 1, I ), 1, DSCALE, DSUM )
               CALL SLASSQ( N, B( 1, I ), 1, DSCALE, DSUM )
   20       CONTINUE
            DIF( 1 ) = DSCALE*SQRT( DSUM )
            DIF( 2 ) = DIF( 1 )
         END IF
         GO TO 60
      END IF

      // Collect the selected blocks at the top-left corner of (A, B).

      KS = 0
      PAIR = .FALSE.
      DO 30 K = 1, N
         IF( PAIR ) THEN
            PAIR = .FALSE.
         } else {

            SWAP = SELECT( K )
            IF( K.LT.N ) THEN
               IF( A( K+1, K ).NE.ZERO ) THEN
                  PAIR = .TRUE.
                  SWAP = SWAP .OR. SELECT( K+1 )
               END IF
            END IF

            IF( SWAP ) THEN
               KS = KS + 1

               // Swap the K-th block to position KS.
               // Perform the reordering of diagonal blocks in (A, B)
               // by orthogonal transformation matrices and update
               // Q and Z accordingly (if requested):

               KK = K
               IF( K.NE.KS ) CALL STGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, KK, KS, WORK, LWORK, IERR )

               IF( IERR.GT.0 ) THEN

                  // Swap is rejected: exit.

                  INFO = 1
                  IF( WANTP ) THEN
                     PL = ZERO
                     PR = ZERO
                  END IF
                  IF( WANTD ) THEN
                     DIF( 1 ) = ZERO
                     DIF( 2 ) = ZERO
                  END IF
                  GO TO 60
               END IF

               IF( PAIR ) KS = KS + 1
            END IF
         END IF
   30 CONTINUE
      IF( WANTP ) THEN

         // Solve generalized Sylvester equation for R and L
         // and compute PL and PR.

         N1 = M
         N2 = N - M
         I = N1 + 1
         IJB = 0
         CALL SLACPY( 'Full', N1, N2, A( 1, I ), LDA, WORK, N1 )
         CALL SLACPY( 'Full', N1, N2, B( 1, I ), LDB, WORK( N1*N2+1 ), N1 )          CALL STGSYL( 'N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR )

         // Estimate the reciprocal of norms of "projections" onto left
         // and right eigenspaces.

         RDSCAL = ZERO
         DSUM = ONE
         CALL SLASSQ( N1*N2, WORK, 1, RDSCAL, DSUM )
         PL = RDSCAL*SQRT( DSUM )
         IF( PL.EQ.ZERO ) THEN
            PL = ONE
         } else {
            PL = DSCALE / ( SQRT( DSCALE*DSCALE / PL+PL )*SQRT( PL ) )
         END IF
         RDSCAL = ZERO
         DSUM = ONE
         CALL SLASSQ( N1*N2, WORK( N1*N2+1 ), 1, RDSCAL, DSUM )
         PR = RDSCAL*SQRT( DSUM )
         IF( PR.EQ.ZERO ) THEN
            PR = ONE
         } else {
            PR = DSCALE / ( SQRT( DSCALE*DSCALE / PR+PR )*SQRT( PR ) )
         END IF
      END IF

      IF( WANTD ) THEN

         // Compute estimates of Difu and Difl.

         IF( WANTD1 ) THEN
            N1 = M
            N2 = N - M
            I = N1 + 1
            IJB = IDIFJB

            // Frobenius norm-based Difu-estimate.

            CALL STGSYL( 'N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR )

            // Frobenius norm-based Difl-estimate.

            CALL STGSYL( 'N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR )
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

   40       CONTINUE
            CALL SLACN2( MN2, WORK( MN2+1 ), WORK, IWORK, DIF( 1 ), KASE, ISAVE )
            IF( KASE.NE.0 ) THEN
               IF( KASE.EQ.1 ) THEN

                  // Solve generalized Sylvester equation.

                  CALL STGSYL( 'N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR )
               } else {

                  // Solve the transposed variant.

                  CALL STGSYL( 'T', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR )
               END IF
               GO TO 40
            END IF
            DIF( 1 ) = DSCALE / DIF( 1 )

            // 1-norm-based estimate of Difl.

   50       CONTINUE
            CALL SLACN2( MN2, WORK( MN2+1 ), WORK, IWORK, DIF( 2 ), KASE, ISAVE )
            IF( KASE.NE.0 ) THEN
               IF( KASE.EQ.1 ) THEN

                  // Solve generalized Sylvester equation.

                  CALL STGSYL( 'N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR )
               } else {

                  // Solve the transposed variant.

                  CALL STGSYL( 'T', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( 2*N1*N2+1 ), LWORK-2*N1*N2, IWORK, IERR )
               END IF
               GO TO 50
            END IF
            DIF( 2 ) = DSCALE / DIF( 2 )

         END IF
      END IF

   60 CONTINUE

      // Compute generalized eigenvalues of reordered pair (A, B) and
      // normalize the generalized Schur form.

      PAIR = .FALSE.
      DO 70 K = 1, N
         IF( PAIR ) THEN
            PAIR = .FALSE.
         } else {

            IF( K.LT.N ) THEN
               IF( A( K+1, K ).NE.ZERO ) THEN
                  PAIR = .TRUE.
               END IF
            END IF

            IF( PAIR ) THEN

              // Compute the eigenvalue(s) at position K.

               WORK( 1 ) = A( K, K )
               WORK( 2 ) = A( K+1, K )
               WORK( 3 ) = A( K, K+1 )
               WORK( 4 ) = A( K+1, K+1 )
               WORK( 5 ) = B( K, K )
               WORK( 6 ) = B( K+1, K )
               WORK( 7 ) = B( K, K+1 )
               WORK( 8 ) = B( K+1, K+1 )
               CALL SLAG2( WORK, 2, WORK( 5 ), 2, SMLNUM*EPS, BETA( K ), BETA( K+1 ), ALPHAR( K ), ALPHAR( K+1 ), ALPHAI( K ) )
               ALPHAI( K+1 ) = -ALPHAI( K )

            } else {

               IF( SIGN( ONE, B( K, K ) ).LT.ZERO ) THEN

                  // If B(K,K) is negative, make it positive

                  DO 80 I = 1, N
                     A( K, I ) = -A( K, I )
                     B( K, I ) = -B( K, I )
                     IF( WANTQ ) Q( I, K ) = -Q( I, K )
   80             CONTINUE
               END IF

               ALPHAR( K ) = A( K, K )
               ALPHAI( K ) = ZERO
               BETA( K ) = B( K, K )

            END IF
         END IF
   70 CONTINUE

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of STGSEN

      }
