      SUBROUTINE ZPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AB( LDAB, * ), B( LDB, * )
      // ..

*  =====================================================================

      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZPBTRF, ZPBTRS
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
      } else if ( KD.LT.0 ) {
         INFO = -3
      } else if ( NRHS.LT.0 ) {
         INFO = -4
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -6
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZPBSV ', -INFO )
         RETURN
      }

      // Compute the Cholesky factorization A = U**H *U or A = L*L**H.

      CALL ZPBTRF( UPLO, N, KD, AB, LDAB, INFO )
      if ( INFO.EQ.0 ) {

         // Solve the system A*X = B, overwriting B with X.

         CALL ZPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )

      }
      RETURN

      // End of ZPBSV

      }
