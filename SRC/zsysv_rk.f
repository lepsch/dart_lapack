      SUBROUTINE ZSYSV_RK( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDB, LWORK, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), E( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LWKOPT;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZSYTRF_RK, ZSYTRS_3
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

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
         INFO = -9
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -11
      END IF

      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            CALL ZSYTRF_RK( UPLO, N, A, LDA, E, IPIV, WORK, -1, INFO )
            LWKOPT = INT( DBLE( WORK( 1 ) ) )
         END IF
         WORK( 1 ) = LWKOPT
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYSV_RK ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Compute the factorization A = P*U*D*(U**T)*(P**T) or
      // A = P*U*D*(U**T)*(P**T).

      CALL ZSYTRF_RK( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )

      IF( INFO.EQ.0 ) THEN

         // Solve the system A*X = B with BLAS3 solver, overwriting B with X.

         CALL ZSYTRS_3( UPLO, N, NRHS, A, LDA, E, IPIV, B, LDB, INFO )

      END IF

      WORK( 1 ) = LWKOPT

      RETURN

      // End of ZSYSV_RK

      END
