      SUBROUTINE SSYSV_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKMIN, LWKOPT, LWKOPT_SYTRF, LWKOPT_SYTRS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      REAL               SROUNDUP_LWORK
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, SSYTRS_AA, SSYTRF_AA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      LWKMIN = MAX( 1, 2*N, 3*N-2 )
      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) {
         INFO = -10
      }

      if ( INFO.EQ.0 ) {
         CALL SSYTRF_AA( UPLO, N, A, LDA, IPIV, WORK, -1, INFO )
         LWKOPT_SYTRF = INT( WORK( 1 ) )
         CALL SSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, -1, INFO )
         LWKOPT_SYTRS = INT( WORK( 1 ) )
         LWKOPT = MAX( LWKMIN, LWKOPT_SYTRF, LWKOPT_SYTRS )
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SSYSV_AA', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Compute the factorization A = U**T*T*U or A = L*T*L**T.

      CALL SSYTRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
      if ( INFO.EQ.0 ) {

         // Solve the system A*X = B, overwriting B with X.

         CALL SSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )

      }

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

      RETURN

      // End of SSYSV_AA

      }
