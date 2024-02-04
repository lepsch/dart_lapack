      void dpbsvx(FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), BERR( * ), FERR( * ), S( * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, RCEQU, UPPER;
      int                I, INFEQU, J, J1, J2;
      double             AMAX, ANORM, BIGNUM, SCOND, SMAX, SMIN, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLANSB;
      // EXTERNAL lsame, DLAMCH, DLANSB
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLACPY, DLAQSB, DPBCON, DPBEQU, DPBRFS, DPBTRF, DPBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      EQUIL = lsame( FACT, 'E' );
      UPPER = lsame( UPLO, 'U' );
      if ( NOFACT || EQUIL ) {
         EQUED = 'N';
         RCEQU = false;
      } else {
         RCEQU = lsame( EQUED, 'Y' );
         SMLNUM = DLAMCH( 'Safe minimum' );
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
         xerbla('DPBSVX', -INFO );
         return;
      }

      if ( EQUIL ) {

         // Compute row and column scalings to equilibrate the matrix A.

         dpbequ(UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFEQU );
         if ( INFEQU == 0 ) {

            // Equilibrate the matrix.

            dlaqsb(UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED );
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

         // Compute the Cholesky factorization A = U**T *U or A = L*L**T.

         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 40
               J1 = max( J-KD, 1 );
               dcopy(J-J1+1, AB( KD+1-J+J1, J ), 1, AFB( KD+1-J+J1, J ), 1 );
            } // 40
         } else {
            for (J = 1; J <= N; J++) { // 50
               J2 = min( J+KD, N );
               dcopy(J2-J+1, AB( 1, J ), 1, AFB( 1, J ), 1 );
            } // 50
         }

         dpbtrf(UPLO, N, KD, AFB, LDAFB, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO;
            return;
         }
      }

      // Compute the norm of the matrix A.

      ANORM = DLANSB( '1', UPLO, N, KD, AB, LDAB, WORK );

      // Compute the reciprocal of the condition number of A.

      dpbcon(UPLO, N, KD, AFB, LDAFB, ANORM, RCOND, WORK, IWORK, INFO );

      // Compute the solution matrix X.

      dlacpy('Full', N, NRHS, B, LDB, X, LDX );
      dpbtrs(UPLO, N, KD, NRHS, AFB, LDAFB, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      dpbrfs(UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO );

      // Transform the solution matrix X to a solution of the original
      // system.

      if ( RCEQU ) {
         for (J = 1; J <= NRHS; J++) { // 70
            for (I = 1; I <= N; I++) { // 60
               X[I, J] = S( I )*X( I, J );
            } // 60
         } // 70
         for (J = 1; J <= NRHS; J++) { // 80
            FERR[J] = FERR( J ) / SCOND;
         } // 80
      }

      // Set INFO = N+1 if the matrix is singular to working precision.

      if( RCOND < DLAMCH( 'Epsilon' ) ) INFO = N + 1;

      return;
      }