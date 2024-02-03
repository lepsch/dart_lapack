      SUBROUTINE CSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKOPT;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CSYTRF, CSYTRS, CSYTRS2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      LQUERY = ( LWORK == -1 )
      if ( .NOT.LSAME( UPLO, 'U' ) && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( LWORK.LT.1 && .NOT.LQUERY ) {
         INFO = -10
      }

      if ( INFO == 0 ) {
         if ( N == 0 ) {
            LWKOPT = 1
         } else {
            csytrf(UPLO, N, A, LDA, IPIV, WORK, -1, INFO );
            LWKOPT = INT( WORK( 1 ) )
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      }

      if ( INFO != 0 ) {
         xerbla('CSYSV ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Compute the factorization A = U*D*U**T or A = L*D*L**T.

      csytrf(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO );
      if ( INFO == 0 ) {

         // Solve the system A*X = B, overwriting B with X.

         if ( LWORK.LT.N ) {

         // Solve with TRS ( Use Level BLAS 2)

            csytrs(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO );

         } else {

         // Solve with TRS2 ( Use Level BLAS 3)

            csytrs2(UPLO,N,NRHS,A,LDA,IPIV,B,LDB,WORK,INFO );

         }

      }

      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)

      RETURN

      // End of CSYSV

      }
