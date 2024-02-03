      SUBROUTINE ZPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      double             BERR( * ), FERR( * ), RWORK( * ), S( * );
      COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, RCEQU, UPPER;
      int                I, INFEQU, J, J1, J2;
      double             AMAX, ANORM, BIGNUM, SCOND, SMAX, SMIN, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANHB;
      // EXTERNAL LSAME, DLAMCH, ZLANHB
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZLACPY, ZLAQHB, ZPBCON, ZPBEQU, ZPBRFS, ZPBTRF, ZPBTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      UPPER = LSAME( UPLO, 'U' )
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
      } else if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( KD < 0 ) {
         INFO = -4
      } else if ( NRHS < 0 ) {
         INFO = -5
      } else if ( LDAB < KD+1 ) {
         INFO = -7
      } else if ( LDAFB < KD+1 ) {
         INFO = -9
      } else if ( LSAME( FACT, 'F' ) && .NOT. ( RCEQU || LSAME( EQUED, 'N' ) ) ) {
         INFO = -10
      } else {
         if ( RCEQU ) {
            SMIN = BIGNUM
            SMAX = ZERO
            for (J = 1; J <= N; J++) { // 10
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
            } // 10
            if ( SMIN <= ZERO ) {
               INFO = -11
            } else if ( N > 0 ) {
               SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
            } else {
               SCOND = ONE
            }
         }
         if ( INFO == 0 ) {
            if ( LDB < MAX( 1, N ) ) {
               INFO = -13
            } else if ( LDX < MAX( 1, N ) ) {
               INFO = -15
            }
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZPBSVX', -INFO );
         RETURN
      }

      if ( EQUIL ) {

         // Compute row and column scalings to equilibrate the matrix A.

         zpbequ(UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFEQU );
         if ( INFEQU == 0 ) {

            // Equilibrate the matrix.

            zlaqhb(UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED );
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

         // Compute the Cholesky factorization A = U**H *U or A = L*L**H.

         if ( UPPER ) {
            for (J = 1; J <= N; J++) { // 40
               J1 = MAX( J-KD, 1 )
               zcopy(J-J1+1, AB( KD+1-J+J1, J ), 1, AFB( KD+1-J+J1, J ), 1 );
            } // 40
         } else {
            for (J = 1; J <= N; J++) { // 50
               J2 = MIN( J+KD, N )
               zcopy(J2-J+1, AB( 1, J ), 1, AFB( 1, J ), 1 );
            } // 50
         }

         zpbtrf(UPLO, N, KD, AFB, LDAFB, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO
            RETURN
         }
      }

      // Compute the norm of the matrix A.

      ANORM = ZLANHB( '1', UPLO, N, KD, AB, LDAB, RWORK )

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
               X( I, J ) = S( I )*X( I, J )
            } // 60
         } // 70
         for (J = 1; J <= NRHS; J++) { // 80
            FERR( J ) = FERR( J ) / SCOND
         } // 80
      }

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND < DLAMCH( 'Epsilon' ) ) INFO = N + 1

      RETURN

      // End of ZPBSVX

      }
