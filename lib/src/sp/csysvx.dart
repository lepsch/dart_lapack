      void csysvx(FACT, UPLO, N, NRHS, final Matrix<double> A, final int LDA, final Matrix<double> AF, final int LDAF, final Array<int> IPIV, final Matrix<double> B, final int LDB, final Matrix<double> X, final int LDX, RCOND, FERR, BERR, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS;
      double               RCOND;
      int                IPIV( * );
      double               BERR( * ), FERR( * ), RWORK( * );
      Complex            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      bool               LQUERY, NOFACT;
      int                LWKOPT, NB;
      double               ANORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               CLANSY, SLAMCH, SROUNDUP_LWORK;
      // EXTERNAL ILAENV, lsame, CLANSY, SLAMCH, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACPY, CSYCON, CSYRFS, CSYTRF, CSYTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      LQUERY = ( LWORK == -1 );
      if ( !NOFACT && !lsame( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( !lsame( UPLO, 'U' ) && !lsame( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDAF < max( 1, N ) ) {
         INFO = -8;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -11;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -13;
      } else if ( LWORK < max( 1, 2*N ) && !LQUERY ) {
         INFO = -18;
      }

      if ( INFO == 0 ) {
         LWKOPT = max( 1, 2*N );
         if ( NOFACT ) {
            NB = ilaenv( 1, 'CSYTRF', UPLO, N, -1, -1, -1 );
            LWKOPT = max( LWKOPT, N*NB );
         }
         WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }

      if ( INFO != 0 ) {
         xerbla('CSYSVX', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      if ( NOFACT ) {

         // Compute the factorization A = U*D*U**T or A = L*D*L**T.

         clacpy(UPLO, N, N, A, LDA, AF, LDAF );
         csytrf(UPLO, N, AF, LDAF, IPIV, WORK, LWORK, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO;
            return;
         }
      }

      // Compute the norm of the matrix A.

      ANORM = CLANSY( 'I', UPLO, N, A, LDA, RWORK );

      // Compute the reciprocal of the condition number of A.

      csycon(UPLO, N, AF, LDAF, IPIV, ANORM, RCOND, WORK, INFO );

      // Compute the solution vectors X.

      clacpy('Full', N, NRHS, B, LDB, X, LDX );
      csytrs(UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      csyrfs(UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

      // Set INFO = N+1 if the matrix is singular to working precision.

      if( RCOND < SLAMCH( 'Epsilon' ) ) INFO = N + 1;

      WORK[1] = SROUNDUP_LWORK(LWKOPT);

      }
