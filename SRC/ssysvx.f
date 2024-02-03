      SUBROUTINE SSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, NOFACT;
      int                LWKMIN, LWKOPT, NB;
      REAL               ANORM
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANSY, SROUNDUP_LWORK
      // EXTERNAL ILAENV, LSAME, SLAMCH, SLANSY, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY, SSYCON, SSYRFS, SSYTRF, SSYTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
      LWKMIN = MAX( 1, 3*N )
      if ( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( NRHS.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDAF.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -11
      } else if ( LDX.LT.MAX( 1, N ) ) {
         INFO = -13
      } else if ( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) {
         INFO = -18
      }

      if ( INFO.EQ.0 ) {
         LWKOPT = LWKMIN
         if ( NOFACT ) {
            NB = ILAENV( 1, 'SSYTRF', UPLO, N, -1, -1, -1 )
            LWKOPT = MAX( LWKOPT, N*NB )
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SSYSVX', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      if ( NOFACT ) {

         // Compute the factorization A = U*D*U**T or A = L*D*L**T.

         CALL SLACPY( UPLO, N, N, A, LDA, AF, LDAF )
         CALL SSYTRF( UPLO, N, AF, LDAF, IPIV, WORK, LWORK, INFO )

         // Return if INFO is non-zero.

         if ( INFO.GT.0 ) {
            RCOND = ZERO
            RETURN
         }
      }

      // Compute the norm of the matrix A.

      ANORM = SLANSY( 'I', UPLO, N, A, LDA, WORK )

      // Compute the reciprocal of the condition number of A.

      CALL SSYCON( UPLO, N, AF, LDAF, IPIV, ANORM, RCOND, WORK, IWORK, INFO )

      // Compute the solution vectors X.

      CALL SLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL SSYTRS( UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO )

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      CALL SSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND.LT.SLAMCH( 'Epsilon' ) ) INFO = N + 1

      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)

      RETURN

      // End of SSYSVX

      }
