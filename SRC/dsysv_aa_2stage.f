      void dsysv_aa_2stage(UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, IPIV2, B, LDB, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // IMPLICIT NONE

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, NRHS, LDA, LDB, LTB, LWORK, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IPIV2( * );
      double             A( LDA, * ), B( LDB, * ), TB( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER, TQUERY, WQUERY;
      int                LWKMIN, LWKOPT;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSYTRF_AA_2STAGE, DSYTRS_AA_2STAGE, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      WQUERY = ( LWORK == -1 );
      TQUERY = ( LTB == -1 );
      LWKMIN = max( 1, N );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LTB < max( 1, 4*N ) && !TQUERY ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -11;
      } else if ( LWORK < LWKMIN && !WQUERY ) {
         INFO = -13;
      }

      if ( INFO == 0 ) {
         dsytrf_aa_2stage(UPLO, N, A, LDA, TB, -1, IPIV, IPIV2, WORK, -1, INFO );
         LWKOPT = max( LWKMIN, INT( WORK( 1 ) ) );
         WORK( 1 ) = LWKOPT;
      }

      if ( INFO != 0 ) {
         xerbla('DSYSV_AA_2STAGE', -INFO );
         return;
      } else if ( WQUERY || TQUERY ) {
         return;
      }

      // Compute the factorization A = U**T*T*U or A = L*T*L**T.

      dsytrf_aa_2stage(UPLO, N, A, LDA, TB, LTB, IPIV, IPIV2, WORK, LWORK, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         dsytrs_aa_2stage(UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, IPIV2, B, LDB, INFO );

      }

      WORK( 1 ) = LWKOPT;

      return;

      // End of DSYSV_AA_2STAGE

      }
