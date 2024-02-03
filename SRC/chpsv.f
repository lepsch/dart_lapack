      SUBROUTINE CHPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            AP( * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPTRF, CHPTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CHPSV ', -INFO )
         RETURN
      }

      // Compute the factorization A = U*D*U**H or A = L*D*L**H.

      CALL CHPTRF( UPLO, N, AP, IPIV, INFO )
      if ( INFO.EQ.0 ) {

         // Solve the system A*X = B, overwriting B with X.

         CALL CHPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )

      }
      RETURN

      // End of CHPSV

      }
