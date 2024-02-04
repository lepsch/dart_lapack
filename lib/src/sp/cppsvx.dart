      void cppsvx(FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      double               RCOND;
      // ..
      // .. Array Arguments ..
      double               BERR( * ), FERR( * ), RWORK( * ), S( * );
      Complex            AFP( * ), AP( * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, RCEQU;
      int                I, INFEQU, J;
      double               AMAX, ANORM, BIGNUM, SCOND, SMAX, SMIN, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANHP, SLAMCH;
      // EXTERNAL lsame, CLANHP, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLACPY, CLAQHP, CPPCON, CPPEQU, CPPRFS, CPPTRF, CPPTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      EQUIL = lsame( FACT, 'E' );
      if ( NOFACT || EQUIL ) {
         EQUED = 'N';
         RCEQU = false;
      } else {
         RCEQU = lsame( EQUED, 'Y' );
         SMLNUM = SLAMCH( 'Safe minimum' );
         BIGNUM = ONE / SMLNUM;
      }

      // Test the input parameters.

      if ( !NOFACT && !EQUIL && !lsame( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( !lsame( UPLO, 'U' ) && !lsame( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( lsame( FACT, 'F' ) && !( RCEQU || lsame( EQUED, 'N' ) ) ) {
         INFO = -7;
      } else {
         if ( RCEQU ) {
            SMIN = BIGNUM;
            SMAX = ZERO;
            for (J = 1; J <= N; J++) { // 10
               SMIN = min( SMIN, S( J ) );
               SMAX = max( SMAX, S( J ) );
            } // 10
            if ( SMIN <= ZERO ) {
               INFO = -8;
            } else if ( N > 0 ) {
               SCOND = max( SMIN, SMLNUM ) / min( SMAX, BIGNUM );
            } else {
               SCOND = ONE;
            }
         }
         if ( INFO == 0 ) {
            if ( LDB < max( 1, N ) ) {
               INFO = -10;
            } else if ( LDX < max( 1, N ) ) {
               INFO = -12;
            }
         }
      }

      if ( INFO != 0 ) {
         xerbla('CPPSVX', -INFO );
         return;
      }

      if ( EQUIL ) {

         // Compute row and column scalings to equilibrate the matrix A.

         cppequ(UPLO, N, AP, S, SCOND, AMAX, INFEQU );
         if ( INFEQU == 0 ) {

            // Equilibrate the matrix.

            claqhp(UPLO, N, AP, S, SCOND, AMAX, EQUED );
            RCEQU = lsame( EQUED, 'Y' );
         }
      }

      // Scale the right-hand side.

      if ( RCEQU ) {
         for (J = 1; J <= NRHS; J++) { // 30
            for (I = 1; I <= N; I++) { // 20
               B[I, J] = S( I )*B( I, J );
            } // 20
         } // 30
      }

      if ( NOFACT || EQUIL ) {

         // Compute the Cholesky factorization A = U**H * U or A = L * L**H.

         ccopy(N*( N+1 ) / 2, AP, 1, AFP, 1 );
         cpptrf(UPLO, N, AFP, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO;
            return;
         }
      }

      // Compute the norm of the matrix A.

      ANORM = CLANHP( 'I', UPLO, N, AP, RWORK );

      // Compute the reciprocal of the condition number of A.

      cppcon(UPLO, N, AFP, ANORM, RCOND, WORK, RWORK, INFO );

      // Compute the solution matrix X.

      clacpy('Full', N, NRHS, B, LDB, X, LDX );
      cpptrs(UPLO, N, NRHS, AFP, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      cpprfs(UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

      // Transform the solution matrix X to a solution of the original
      // system.

      if ( RCEQU ) {
         for (J = 1; J <= NRHS; J++) { // 50
            for (I = 1; I <= N; I++) { // 40
               X[I, J] = S( I )*X( I, J );
            } // 40
         } // 50
         for (J = 1; J <= NRHS; J++) { // 60
            FERR[J] = FERR( J ) / SCOND;
         } // 60
      }

      // Set INFO = N+1 if the matrix is singular to working precision.

      if( RCOND < SLAMCH( 'Epsilon' ) ) INFO = N + 1;

      return;
      }