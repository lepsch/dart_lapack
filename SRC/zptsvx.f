      SUBROUTINE ZPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             FACT;
      int                INFO, LDB, LDX, N, NRHS;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      double             BERR( * ), D( * ), DF( * ), FERR( * ), RWORK( * );
      COMPLEX*16         B( LDB, * ), E( * ), EF( * ), WORK( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOFACT;
      double             ANORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANHT;
      // EXTERNAL LSAME, DLAMCH, ZLANHT
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, XERBLA, ZCOPY, ZLACPY, ZPTCON, ZPTRFS, ZPTTRF, ZPTTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      if ( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -9
      } else if ( LDX.LT.MAX( 1, N ) ) {
         INFO = -11
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZPTSVX', -INFO )
         RETURN
      }

      if ( NOFACT ) {

         // Compute the L*D*L**H (or U**H*D*U) factorization of A.

         CALL DCOPY( N, D, 1, DF, 1 )
         IF( N.GT.1 ) CALL ZCOPY( N-1, E, 1, EF, 1 )
         CALL ZPTTRF( N, DF, EF, INFO )

         // Return if INFO is non-zero.

         if ( INFO.GT.0 ) {
            RCOND = ZERO
            RETURN
         }
      }

      // Compute the norm of the matrix A.

      ANORM = ZLANHT( '1', N, D, E )

      // Compute the reciprocal of the condition number of A.

      CALL ZPTCON( N, DF, EF, ANORM, RCOND, RWORK, INFO )

      // Compute the solution vectors X.

      CALL ZLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL ZPTTRS( 'Lower', N, NRHS, DF, EF, X, LDX, INFO )

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      CALL ZPTRFS( 'Lower', N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) INFO = N + 1

      RETURN

      // End of ZPTSVX

      }
