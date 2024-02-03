      SUBROUTINE SSYSV_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, IPIV2, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, NRHS, LDA, LDB, LTB, LWORK, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IPIV2( * );
      REAL               A( LDA, * ), B( LDB, * ), TB( * ), WORK( * )
      // ..

*  =====================================================================
      // ..
      // .. Local Scalars ..
      bool               UPPER, TQUERY, WQUERY;
      int                LWKMIN, LWKOPT;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      REAL               SROUNDUP_LWORK
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSYTRF_AA_2STAGE, SSYTRS_AA_2STAGE, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      WQUERY = ( LWORK == -1 )
      TQUERY = ( LTB == -1 )
      LWKMIN = MAX( 1, N )
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( NRHS < 0 ) {
         INFO = -3
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5
      } else if ( LTB < MAX( 1, 4*N ) && !TQUERY ) {
         INFO = -7
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -11
      } else if ( LWORK < LWKMIN && !WQUERY ) {
         INFO = -13
      }

      if ( INFO == 0 ) {
         ssytrf_aa_2stage(UPLO, N, A, LDA, TB, -1, IPIV, IPIV2, WORK, -1, INFO );
         LWKOPT = MAX( LWKMIN, INT( WORK( 1 ) ) )
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      }

      if ( INFO != 0 ) {
         xerbla('SSYSV_AA_2STAGE', -INFO );
         RETURN
      } else if ( WQUERY || TQUERY ) {
         RETURN
      }

      // Compute the factorization A = U**T*T*U or A = L*T*L**T.

      ssytrf_aa_2stage(UPLO, N, A, LDA, TB, LTB, IPIV, IPIV2, WORK, LWORK, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         ssytrs_aa_2stage(UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, IPIV2, B, LDB, INFO );

      }

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

      RETURN

      // End of SSYSV_AA_2STAGE

      }
