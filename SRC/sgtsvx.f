      SUBROUTINE SGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             FACT, TRANS;
      int                INFO, LDB, LDX, N, NRHS;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               B( LDB, * ), BERR( * ), D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), FERR( * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOFACT, NOTRAN;
      String             NORM;
      REAL               ANORM
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANGT
      // EXTERNAL LSAME, SLAMCH, SLANGT
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGTCON, SGTRFS, SGTTRF, SGTTRS, SLACPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      NOTRAN = LSAME( TRANS, 'N' )
      if ( .NOT.NOFACT && .NOT.LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRAN && .NOT.LSAME( TRANS, 'T' ) && .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( NRHS < 0 ) {
         INFO = -4
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -14
      } else if ( LDX < MAX( 1, N ) ) {
         INFO = -16
      }
      if ( INFO != 0 ) {
         xerbla('SGTSVX', -INFO );
         RETURN
      }

      if ( NOFACT ) {

         // Compute the LU factorization of A.

         scopy(N, D, 1, DF, 1 );
         if ( N > 1 ) {
            scopy(N-1, DL, 1, DLF, 1 );
            scopy(N-1, DU, 1, DUF, 1 );
         }
         sgttrf(N, DLF, DF, DUF, DU2, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
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
      ANORM = SLANGT( NORM, N, DL, D, DU )

      // Compute the reciprocal of the condition number of A.

      sgtcon(NORM, N, DLF, DF, DUF, DU2, IPIV, ANORM, RCOND, WORK, IWORK, INFO );

      // Compute the solution vectors X.

      slacpy('Full', N, NRHS, B, LDB, X, LDX );
      sgttrs(TRANS, N, NRHS, DLF, DF, DUF, DU2, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      sgtrfs(TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO );

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND < SLAMCH( 'Epsilon' ) ) INFO = N + 1

      RETURN

      // End of SGTSVX

      }
