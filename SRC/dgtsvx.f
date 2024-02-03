      SUBROUTINE DGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             FACT, TRANS;
      int                INFO, LDB, LDX, N, NRHS;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             B( LDB, * ), BERR( * ), D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), FERR( * ), WORK( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOFACT, NOTRAN;
      String             NORM;
      double             ANORM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANGT;
      // EXTERNAL LSAME, DLAMCH, DLANGT
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGTCON, DGTRFS, DGTTRF, DGTTRS, DLACPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      NOTRAN = LSAME( TRANS, 'N' )
      if ( .NOT.NOFACT .AND. .NOT.LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( NRHS.LT.0 ) {
         INFO = -4
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -14
      } else if ( LDX.LT.MAX( 1, N ) ) {
         INFO = -16
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DGTSVX', -INFO )
         RETURN
      }

      if ( NOFACT ) {

         // Compute the LU factorization of A.

         CALL DCOPY( N, D, 1, DF, 1 )
         if ( N.GT.1 ) {
            CALL DCOPY( N-1, DL, 1, DLF, 1 )
            CALL DCOPY( N-1, DU, 1, DUF, 1 )
         }
         CALL DGTTRF( N, DLF, DF, DUF, DU2, IPIV, INFO )

         // Return if INFO is non-zero.

         if ( INFO.GT.0 ) {
            RCOND = ZERO
            RETURN
         }
      }

      // Compute the norm of the matrix A.

      if ( NOTRAN ) {
         NORM = '1'
      } else {
         NORM = 'I'
      }
      ANORM = DLANGT( NORM, N, DL, D, DU )

      // Compute the reciprocal of the condition number of A.

      CALL DGTCON( NORM, N, DLF, DF, DUF, DU2, IPIV, ANORM, RCOND, WORK, IWORK, INFO )

      // Compute the solution vectors X.

      CALL DLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL DGTTRS( TRANS, N, NRHS, DLF, DF, DUF, DU2, IPIV, X, LDX, INFO )

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      CALL DGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) INFO = N + 1

      RETURN

      // End of DGTSVX

      }
