      SUBROUTINE CTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB, ALPHA, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF, WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      bool               WANTQ, WANTZ;
      int                IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK, M, N
      REAL               PL, PR
*     ..
*     .. Array Arguments ..
      bool               SELECT( * );
      int                IWORK( * )
      REAL               DIF( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                IDIFJB
      PARAMETER          ( IDIFJB = 3 )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      bool               LQUERY, SWAP, WANTD, WANTD1, WANTD2, WANTP;
      int                I, IERR, IJB, K, KASE, KS, LIWMIN, LWMIN, MN2, N1, N2
      REAL               DSCALE, DSUM, RDSCAL, SAFMIN
      COMPLEX            TEMP1, TEMP2
*     ..
*     .. Local Arrays ..
      int                ISAVE( 3 )
*     ..
*     .. External Functions ..
      REAL               SROUNDUP_LWORK
      EXTERNAL           SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      REAL               SLAMCH
      EXTERNAL           CLACN2, CLACPY, CLASSQ, CSCAL, CTGEXC, CTGSYL, SLAMCH, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, CONJG, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      IF( IJOB.LT.0 .OR. IJOB.GT.5 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -13
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -15
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTGSEN', -INFO )
         RETURN
      END IF
*
      IERR = 0
*
      WANTP = IJOB.EQ.1 .OR. IJOB.GE.4
      WANTD1 = IJOB.EQ.2 .OR. IJOB.EQ.4
      WANTD2 = IJOB.EQ.3 .OR. IJOB.EQ.5
      WANTD = WANTD1 .OR. WANTD2
*
*     Set M to the dimension of the specified pair of deflating
*     subspaces.
*
      M = 0
      IF( .NOT.LQUERY .OR. IJOB.NE.0 ) THEN
      DO 10 K = 1, N
         ALPHA( K ) = A( K, K )
         BETA( K ) = B( K, K )
         IF( K.LT.N ) THEN
            IF( SELECT( K ) ) M = M + 1
         ELSE
            IF( SELECT( N ) ) M = M + 1
         END IF
   10 CONTINUE
      END IF
*
      IF( IJOB.EQ.1 .OR. IJOB.EQ.2 .OR. IJOB.EQ.4 ) THEN
         LWMIN = MAX( 1, 2*M*(N-M) )
         LIWMIN = MAX( 1, N+2 )
      ELSE IF( IJOB.EQ.3 .OR. IJOB.EQ.5 ) THEN
         LWMIN = MAX( 1, 4*M*(N-M) )
         LIWMIN = MAX( 1, 2*M*(N-M), N+2 )
      ELSE
         LWMIN = 1
         LIWMIN = 1
      END IF
*
      WORK( 1 ) =  SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN
*
      IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -21
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -23
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTGSEN', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( M.EQ.N .OR. M.EQ.0 ) THEN
         IF( WANTP ) THEN
            PL = ONE
            PR = ONE
         END IF
         IF( WANTD ) THEN
            DSCALE = ZERO
            DSUM = ONE
            DO 20 I = 1, N
               CALL CLASSQ( N, A( 1, I ), 1, DSCALE, DSUM )
               CALL CLASSQ( N, B( 1, I ), 1, DSCALE, DSUM )
   20       CONTINUE
            DIF( 1 ) = DSCALE*SQRT( DSUM )
            DIF( 2 ) = DIF( 1 )
         END IF
         GO TO 70
      END IF
*
*     Get machine constant
*
      SAFMIN = SLAMCH( 'S' )
*
*     Collect the selected blocks at the top-left corner of (A, B).
*
      KS = 0
      DO 30 K = 1, N
         SWAP = SELECT( K )
         IF( SWAP ) THEN
            KS = KS + 1
*
*           Swap the K-th block to position KS. Compute unitary Q
*           and Z that will swap adjacent diagonal blocks in (A, B).
*
            IF( K.NE.KS ) CALL CTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, K, KS, IERR )
*
            IF( IERR.GT.0 ) THEN
*
*              Swap is rejected: exit.
*
               INFO = 1
               IF( WANTP ) THEN
                  PL = ZERO
                  PR = ZERO
               END IF
               IF( WANTD ) THEN
                  DIF( 1 ) = ZERO
                  DIF( 2 ) = ZERO
               END IF
               GO TO 70
            END IF
         END IF
   30 CONTINUE
      IF( WANTP ) THEN
*
*        Solve generalized Sylvester equation for R and L:
*                   A11 * R - L * A22 = A12
*                   B11 * R - L * B22 = B12
*
         N1 = M
         N2 = N - M
         I = N1 + 1
         CALL CLACPY( 'Full', N1, N2, A( 1, I ), LDA, WORK, N1 )
         CALL CLACPY( 'Full', N1, N2, B( 1, I ), LDB, WORK( N1*N2+1 ), N1 )
         IJB = 0
         CALL CTGSYL( 'N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR )
*
*        Estimate the reciprocal of norms of "projections" onto
*        left and right eigenspaces
*
         RDSCAL = ZERO
         DSUM = ONE
         CALL CLASSQ( N1*N2, WORK, 1, RDSCAL, DSUM )
         PL = RDSCAL*SQRT( DSUM )
         IF( PL.EQ.ZERO ) THEN
            PL = ONE
         ELSE
            PL = DSCALE / ( SQRT( DSCALE*DSCALE / PL+PL )*SQRT( PL ) )
         END IF
         RDSCAL = ZERO
         DSUM = ONE
         CALL CLASSQ( N1*N2, WORK( N1*N2+1 ), 1, RDSCAL, DSUM )
         PR = RDSCAL*SQRT( DSUM )
         IF( PR.EQ.ZERO ) THEN
            PR = ONE
         ELSE
            PR = DSCALE / ( SQRT( DSCALE*DSCALE / PR+PR )*SQRT( PR ) )
         END IF
      END IF
      IF( WANTD ) THEN
*
*        Compute estimates Difu and Difl.
*
         IF( WANTD1 ) THEN
            N1 = M
            N2 = N - M
            I = N1 + 1
            IJB = IDIFJB
*
*           Frobenius norm-based Difu estimate.
*
            CALL CTGSYL( 'N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR )
*
*           Frobenius norm-based Difl estimate.
*
            CALL CTGSYL( 'N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR )
         ELSE
*
*           Compute 1-norm-based estimates of Difu and Difl using
*           reversed communication with CLACN2. In each step a
*           generalized Sylvester equation or a transposed variant
*           is solved.
*
            KASE = 0
            N1 = M
            N2 = N - M
            I = N1 + 1
            IJB = 0
            MN2 = 2*N1*N2
*
*           1-norm-based estimate of Difu.
*
   40       CONTINUE
            CALL CLACN2( MN2, WORK( MN2+1 ), WORK, DIF( 1 ), KASE, ISAVE )
            IF( KASE.NE.0 ) THEN
               IF( KASE.EQ.1 ) THEN
*
*                 Solve generalized Sylvester equation
*
                  CALL CTGSYL( 'N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR )
               ELSE
*
*                 Solve the transposed variant.
*
                  CALL CTGSYL( 'C', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK, N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR )
               END IF
               GO TO 40
            END IF
            DIF( 1 ) = DSCALE / DIF( 1 )
*
*           1-norm-based estimate of Difl.
*
   50       CONTINUE
            CALL CLACN2( MN2, WORK( MN2+1 ), WORK, DIF( 2 ), KASE, ISAVE )
            IF( KASE.NE.0 ) THEN
               IF( KASE.EQ.1 ) THEN
*
*                 Solve generalized Sylvester equation
*
                  CALL CTGSYL( 'N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR )
               ELSE
*
*                 Solve the transposed variant.
*
                  CALL CTGSYL( 'C', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK, N2, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK, IERR )
               END IF
               GO TO 50
            END IF
            DIF( 2 ) = DSCALE / DIF( 2 )
         END IF
      END IF
*
*     If B(K,K) is complex, make it real and positive (normalization
*     of the generalized Schur form) and Store the generalized
*     eigenvalues of reordered pair (A, B)
*
      DO 60 K = 1, N
         DSCALE = ABS( B( K, K ) )
         IF( DSCALE.GT.SAFMIN ) THEN
            TEMP1 = CONJG( B( K, K ) / DSCALE )
            TEMP2 = B( K, K ) / DSCALE
            B( K, K ) = DSCALE
            CALL CSCAL( N-K, TEMP1, B( K, K+1 ), LDB )
            CALL CSCAL( N-K+1, TEMP1, A( K, K ), LDA )
            IF( WANTQ ) CALL CSCAL( N, TEMP2, Q( 1, K ), 1 )
         ELSE
            B( K, K ) = CMPLX( ZERO, ZERO )
         END IF
*
         ALPHA( K ) = A( K, K )
         BETA( K ) = B( K, K )
*
   60 CONTINUE
*
   70 CONTINUE
*
      WORK( 1 ) =  SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of CTGSEN
*
      END
