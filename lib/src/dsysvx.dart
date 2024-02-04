      void dsysvx(FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, IWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             A( LDA, * ), AF( LDAF, * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, NOFACT;
      int                LWKMIN, LWKOPT, NB;
      double             ANORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- double             DLAMCH, DLANSY;
      // EXTERNAL lsame, ILAENV, DLAMCH, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACPY, DSYCON, DSYRFS, DSYTRF, DSYTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      LQUERY = ( LWORK == -1 );
      LWKMIN = max( 1, 3*N );
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
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -18;
      }

      if ( INFO == 0 ) {
         LWKOPT = LWKMIN;
         if ( NOFACT ) {
            NB = ILAENV( 1, 'DSYTRF', UPLO, N, -1, -1, -1 );
            LWKOPT = max( LWKOPT, N*NB );
         }
         WORK[1] = LWKOPT;
      }

      if ( INFO != 0 ) {
         xerbla('DSYSVX', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      if ( NOFACT ) {

         // Compute the factorization A = U*D*U**T or A = L*D*L**T.

         dlacpy(UPLO, N, N, A, LDA, AF, LDAF );
         dsytrf(UPLO, N, AF, LDAF, IPIV, WORK, LWORK, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO;
            return;
         }
      }

      // Compute the norm of the matrix A.

      ANORM = DLANSY( 'I', UPLO, N, A, LDA, WORK );

      // Compute the reciprocal of the condition number of A.

      dsycon(UPLO, N, AF, LDAF, IPIV, ANORM, RCOND, WORK, IWORK, INFO );

      // Compute the solution vectors X.

      dlacpy('Full', N, NRHS, B, LDB, X, LDX );
      dsytrs(UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      dsyrfs(UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO );

      // Set INFO = N+1 if the matrix is singular to working precision.

      if( RCOND < DLAMCH( 'Epsilon' ) ) INFO = N + 1;

      WORK[1] = LWKOPT;

      return;
      }