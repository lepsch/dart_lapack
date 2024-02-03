      SUBROUTINE SGGSVD( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, IWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), Q( LDQ, * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      bool               WANTQ, WANTU, WANTV;
      int                I, IBND, ISUB, J, NCYCLE;
      REAL               ANORM, BNORM, SMAX, TEMP, TOLA, TOLB, ULP, UNFL
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANGE
      // EXTERNAL LSAME, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGGSVP, STGSJA, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters
*
      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
      WANTQ = LSAME( JOBQ, 'Q' )
*
      INFO = 0
      IF( .NOT.( WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
         INFO = -16
      ELSE IF( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) THEN
         INFO = -18
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -20
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGGSVD', -INFO )
         RETURN
      END IF
*
      // Compute the Frobenius norm of matrices A and B
*
      ANORM = SLANGE( '1', M, N, A, LDA, WORK )
      BNORM = SLANGE( '1', P, N, B, LDB, WORK )
*
      // Get machine precision and set up threshold for determining
     t // he effective numerical rank of the matrices A and B.
*
      ULP = SLAMCH( 'Precision' )
      UNFL = SLAMCH( 'Safe Minimum' )
      TOLA = MAX( M, N )*MAX( ANORM, UNFL )*ULP
      TOLB = MAX( P, N )*MAX( BNORM, UNFL )*ULP
*
      // Preprocessing
*
      CALL SGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, WORK, WORK( N+1 ), INFO )
*
      // Compute the GSVD of two upper "triangular" matrices
*
      CALL STGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, NCYCLE, INFO )
*
      // Sort the singular values and store the pivot indices in IWORK
      // Copy ALPHA to WORK, then sort ALPHA in WORK
*
      CALL SCOPY( N, ALPHA, 1, WORK, 1 )
      IBND = MIN( L, M-K )
      DO 20 I = 1, IBND
*
         // Scan for largest ALPHA(K+I)
*
         ISUB = I
         SMAX = WORK( K+I )
         DO 10 J = I + 1, IBND
            TEMP = WORK( K+J )
            IF( TEMP.GT.SMAX ) THEN
               ISUB = J
               SMAX = TEMP
            END IF
   10    CONTINUE
         IF( ISUB.NE.I ) THEN
            WORK( K+ISUB ) = WORK( K+I )
            WORK( K+I ) = SMAX
            IWORK( K+I ) = K + ISUB
         ELSE
            IWORK( K+I ) = K + I
         END IF
   20 CONTINUE
*
      RETURN
*
      // End of SGGSVD
*
      END
