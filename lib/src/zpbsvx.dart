      void zpbsvx(FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             EQUED, FACT, UPLO;
      int                INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS;
      double             RCOND;
      double             BERR( * ), FERR( * ), RWORK( * ), S( * );
      Complex         AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               EQUIL, NOFACT, RCEQU, UPPER;
      int                I, INFEQU, J, J1, J2;
      double             AMAX, ANORM, BIGNUM, SCOND, SMAX, SMIN, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, ZLANHB;
      // EXTERNAL lsame, DLAMCH, ZLANHB
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZLACPY, ZLAQHB, ZPBCON, ZPBEQU, ZPBRFS, ZPBTRF, ZPBTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      EQUIL = lsame( FACT, 'E' );
      UPPER = lsame( UPLO, 'U' );
      if ( NOFACT || EQUIL ) {
         EQUED = 'N';
         RCEQU = false;
      } else {
         RCEQU = lsame( EQUED, 'Y' );
         SMLNUM = dlamch( 'Safe minimum' );
         BIGNUM = ONE / SMLNUM;
      }

      // Test the input parameters.

      if ( !NOFACT && !EQUIL && !lsame( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( KD < 0 ) {
         INFO = -4;
      } else if ( NRHS < 0 ) {
         INFO = -5;
      } else if ( LDAB < KD+1 ) {
         INFO = -7;
      } else if ( LDAFB < KD+1 ) {
         INFO = -9;
      } else if ( lsame( FACT, 'F' ) && !( RCEQU || lsame( EQUED, 'N' ) ) ) {
         INFO = -10;
      } else {
         if ( RCEQU ) {
            SMIN = BIGNUM;
            SMAX = ZERO;
            for (J = 1; J <= N; J++) { // 10
               SMIN = min( SMIN, S( J ) );
               SMAX = max( SMAX, S( J ) );
            } // 10
            if ( SMIN <= ZERO ) {
               INFO = -11;
            } else if ( N > 0 ) {
               SCOND = max( SMIN, SMLNUM ) / min( SMAX, BIGNUM );
            } else {
               SCOND = ONE;
            }
         }
         if ( INFO == 0 ) {
            if ( LDB < max( 1, N ) ) {
               INFO = -13;
            } else if ( LDX < max( 1, N ) ) {
               INFO = -15;
            }
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZPBSVX', -INFO );
         return;
      }

      if ( EQUIL ) {

         // Compute row and column scalings to equilibrate the matrix A.

         zpbequ(UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFEQU );
         if ( INFEQU == 0 ) {

            // Equilibrate the matrix.

            zlaqhb(UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED );
            RCEQU = lsame( EQUED, 'Y' );
         }
      }

      // Scale the right-hand side.

      if ( RCEQU ) {
         for (J = 1; J <= NRHS; J++) { // 30
            for (I = 1; I <= N; I++) { // 20
               B[I][J] = S( I )*B( I, J );
            } // 20
         } // 30
      }

      if ( NOFACT || EQUIL ) {

         // Compute the Cholesky factorization A = U**H *U or A = L*L**H.

         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 40
               J1 = max( J-KD, 1 );
               zcopy(J-J1+1, AB( KD+1-J+J1, J ), 1, AFB( KD+1-J+J1, J ), 1 );
            } // 40
         } else {
            for (J = 1; J <= N; J++) { // 50
               J2 = min( J+KD, N );
               zcopy(J2-J+1, AB( 1, J ), 1, AFB( 1, J ), 1 );
            } // 50
         }

         zpbtrf(UPLO, N, KD, AFB, LDAFB, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO;
            return;
         }
      }

      // Compute the norm of the matrix A.

      ANORM = ZLANHB( '1', UPLO, N, KD, AB, LDAB, RWORK );

      // Compute the reciprocal of the condition number of A.

      zpbcon(UPLO, N, KD, AFB, LDAFB, ANORM, RCOND, WORK, RWORK, INFO );

      // Compute the solution matrix X.

      zlacpy('Full', N, NRHS, B, LDB, X, LDX );
      zpbtrs(UPLO, N, KD, NRHS, AFB, LDAFB, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      zpbrfs(UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

      // Transform the solution matrix X to a solution of the original
      // system.

      if ( RCEQU ) {
         for (J = 1; J <= NRHS; J++) { // 70
            for (I = 1; I <= N; I++) { // 60
               X[I][J] = S( I )*X( I, J );
            } // 60
         } // 70
         for (J = 1; J <= NRHS; J++) { // 80
            FERR[J] = FERR( J ) / SCOND;
         } // 80
      }

      // Set INFO = N+1 if the matrix is singular to working precision.

      if( RCOND < dlamch( 'Epsilon' ) ) INFO = N + 1;

      }
