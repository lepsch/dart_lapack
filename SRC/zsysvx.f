      SUBROUTINE ZSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             BERR( * ), FERR( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, NOFACT;
      int                LWKOPT, NB;
      double             ANORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, ZLANSY;
      // EXTERNAL LSAME, ILAENV, DLAMCH, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACPY, ZSYCON, ZSYRFS, ZSYTRF, ZSYTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      LQUERY = ( LWORK == -1 )
      if ( .NOT.NOFACT && .NOT.LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( UPLO, 'U' ) && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( NRHS < 0 ) {
         INFO = -4
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDAF < MAX( 1, N ) ) {
         INFO = -8
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -11
      } else if ( LDX < MAX( 1, N ) ) {
         INFO = -13
      } else if ( LWORK < MAX( 1, 2*N ) && .NOT.LQUERY ) {
         INFO = -18
      }

      if ( INFO == 0 ) {
         LWKOPT = MAX( 1, 2*N )
         if ( NOFACT ) {
            NB = ILAENV( 1, 'ZSYTRF', UPLO, N, -1, -1, -1 )
            LWKOPT = MAX( LWKOPT, N*NB )
         }
         WORK( 1 ) = LWKOPT
      }

      if ( INFO != 0 ) {
         xerbla('ZSYSVX', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      if ( NOFACT ) {

         // Compute the factorization A = U*D*U**T or A = L*D*L**T.

         zlacpy(UPLO, N, N, A, LDA, AF, LDAF );
         zsytrf(UPLO, N, AF, LDAF, IPIV, WORK, LWORK, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO
            RETURN
         }
      }

      // Compute the norm of the matrix A.

      ANORM = ZLANSY( 'I', UPLO, N, A, LDA, RWORK )

      // Compute the reciprocal of the condition number of A.

      zsycon(UPLO, N, AF, LDAF, IPIV, ANORM, RCOND, WORK, INFO );

      // Compute the solution vectors X.

      zlacpy('Full', N, NRHS, B, LDB, X, LDX );
      zsytrs(UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      zsyrfs(UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND < DLAMCH( 'Epsilon' ) ) INFO = N + 1

      WORK( 1 ) = LWKOPT

      RETURN

      // End of ZSYSVX

      }
