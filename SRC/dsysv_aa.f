      SUBROUTINE DSYSV_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), B( LDB, * ), WORK( * );
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKMIN, LWKOPT, LWKOPT_SYTRF, LWKOPT_SYTRS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL ILAENV, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DSYTRF_AA, DSYTRS_AA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      LWKMIN = MAX( 1, 2*N, 3*N-2 )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
*
      IF( INFO.EQ.0 ) THEN
         CALL DSYTRF_AA( UPLO, N, A, LDA, IPIV, WORK, -1, INFO )
         LWKOPT_SYTRF = INT( WORK( 1 ) )
         CALL DSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, -1, INFO )
         LWKOPT_SYTRS = INT( WORK( 1 ) )
         LWKOPT = MAX( LWKMIN, LWKOPT_SYTRF, LWKOPT_SYTRS )
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYSV_AA ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      // Compute the factorization A = U**T*T*U or A = L*T*L**T.
*
      CALL DSYTRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
      IF( INFO.EQ.0 ) THEN
*
         // Solve the system A*X = B, overwriting B with X.
*
         CALL DSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
*
      END IF
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
      // End of DSYSV_AA
*
      END
