      SUBROUTINE ZPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AP( * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZPPTRF, ZPPTRS
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
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('ZPPSV ', -INFO );
         RETURN
      }

      // Compute the Cholesky factorization A = U**H *U or A = L*L**H.

      zpptrf(UPLO, N, AP, INFO );
      if ( INFO.EQ.0 ) {

         // Solve the system A*X = B, overwriting B with X.

         zpptrs(UPLO, N, NRHS, AP, B, LDB, INFO );

      }
      RETURN

      // End of ZPPSV

      }
