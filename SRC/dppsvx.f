      SUBROUTINE DPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             AFP( * ), AP( * ), B( LDB, * ), BERR( * ), FERR( * ), S( * ), WORK( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, RCEQU;
      int                I, INFEQU, J;
      double             AMAX, ANORM, BIGNUM, SCOND, SMAX, SMIN, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANSP;
      // EXTERNAL LSAME, DLAMCH, DLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLACPY, DLAQSP, DPPCON, DPPEQU, DPPRFS, DPPTRF, DPPTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      if ( NOFACT || EQUIL ) {
         EQUED = 'N'
         RCEQU = false;
      } else {
         RCEQU = LSAME( EQUED, 'Y' )
         SMLNUM = DLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
      }

      // Test the input parameters.

      if ( .NOT.NOFACT && .NOT.EQUIL && .NOT.LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( UPLO, 'U' ) && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( NRHS < 0 ) {
         INFO = -4
      } else if ( LSAME( FACT, 'F' ) && .NOT. ( RCEQU || LSAME( EQUED, 'N' ) ) ) {
         INFO = -7
      } else {
         if ( RCEQU ) {
            SMIN = BIGNUM
            SMAX = ZERO
            for (J = 1; J <= N; J++) { // 10
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
            } // 10
            if ( SMIN.LE.ZERO ) {
               INFO = -8
            } else if ( N > 0 ) {
               SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
            } else {
               SCOND = ONE
            }
         }
         if ( INFO == 0 ) {
            if ( LDB < MAX( 1, N ) ) {
               INFO = -10
            } else if ( LDX < MAX( 1, N ) ) {
               INFO = -12
            }
         }
      }

      if ( INFO != 0 ) {
         xerbla('DPPSVX', -INFO );
         RETURN
      }

      if ( EQUIL ) {

         // Compute row and column scalings to equilibrate the matrix A.

         dppequ(UPLO, N, AP, S, SCOND, AMAX, INFEQU );
         if ( INFEQU == 0 ) {

            // Equilibrate the matrix.

            dlaqsp(UPLO, N, AP, S, SCOND, AMAX, EQUED );
            RCEQU = LSAME( EQUED, 'Y' )
         }
      }

      // Scale the right-hand side.

      if ( RCEQU ) {
         for (J = 1; J <= NRHS; J++) { // 30
            for (I = 1; I <= N; I++) { // 20
               B( I, J ) = S( I )*B( I, J )
            } // 20
         } // 30
      }

      if ( NOFACT || EQUIL ) {

         // Compute the Cholesky factorization A = U**T * U or A = L * L**T.

         dcopy(N*( N+1 ) / 2, AP, 1, AFP, 1 );
         dpptrf(UPLO, N, AFP, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO
            RETURN
         }
      }

      // Compute the norm of the matrix A.

      ANORM = DLANSP( 'I', UPLO, N, AP, WORK )

      // Compute the reciprocal of the condition number of A.

      dppcon(UPLO, N, AFP, ANORM, RCOND, WORK, IWORK, INFO );

      // Compute the solution matrix X.

      dlacpy('Full', N, NRHS, B, LDB, X, LDX );
      dpptrs(UPLO, N, NRHS, AFP, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      dpprfs(UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO );

      // Transform the solution matrix X to a solution of the original
      // system.

      if ( RCEQU ) {
         for (J = 1; J <= NRHS; J++) { // 50
            for (I = 1; I <= N; I++) { // 40
               X( I, J ) = S( I )*X( I, J )
            } // 40
         } // 50
         for (J = 1; J <= NRHS; J++) { // 60
            FERR( J ) = FERR( J ) / SCOND
         } // 60
      }

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND < DLAMCH( 'Epsilon' ) ) INFO = N + 1

      RETURN

      // End of DPPSVX

      }
