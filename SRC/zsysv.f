      SUBROUTINE ZSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
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
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKOPT;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZSYTRF, ZSYTRS, ZSYTRS2
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
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            CALL ZSYTRF( UPLO, N, A, LDA, IPIV, WORK, -1, INFO )
            LWKOPT = INT( DBLE( WORK( 1 ) ) )
         END IF
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYSV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      // Compute the factorization A = U*D*U**T or A = L*D*L**T.
*
      CALL ZSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
      IF( INFO.EQ.0 ) THEN
*
         // Solve the system A*X = B, overwriting B with X.
*
         IF ( LWORK.LT.N ) THEN
*
         // Solve with TRS ( Use Level BLAS 2)
*
            CALL ZSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
         ELSE
*
         // Solve with TRS2 ( Use Level BLAS 3)
*
            CALL ZSYTRS2( UPLO,N,NRHS,A,LDA,IPIV,B,LDB,WORK,INFO )
*
         END IF
*
      END IF
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
      // End of ZSYSV
*
      END
