      SUBROUTINE ZCPOSV( UPLO, N, NRHS, A, LDA, B, LDB, X, LDX, WORK, SWORK, RWORK, ITER, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, ITER, LDA, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX            SWORK( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( N, * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      bool               DOITREF;
      const              DOITREF = true ;

      int                ITERMAX;
      const              ITERMAX = 30 ;

      double             BWDMAX;
      const              BWDMAX = 1.0E+00 ;

      COMPLEX*16         NEGONE, ONE
      const              NEGONE = ( -1.0D+00, 0.0D+00 ), ONE = ( 1.0D+00, 0.0D+00 ) ;

      // .. Local Scalars ..
      int                I, IITER, PTSA, PTSX;
      double             ANRM, CTE, EPS, RNRM, XNRM;
      COMPLEX*16         ZDUM

      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZHEMM, ZLACPY, ZLAT2C, ZLAG2C, CLAG2Z, CPOTRF, CPOTRS, XERBLA, ZPOTRF, ZPOTRS
      // ..
      // .. External Functions ..
      int                IZAMAX;
      double             DLAMCH, ZLANHE;
      bool               LSAME;
      // EXTERNAL IZAMAX, DLAMCH, ZLANHE, LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, SQRT
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      INFO = 0
      ITER = 0

      // Test the input parameters.

      if ( .NOT.LSAME( UPLO, 'U' ) && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( NRHS < 0 ) {
         INFO = -3
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDX < MAX( 1, N ) ) {
         INFO = -9
      }
      if ( INFO != 0 ) {
         xerbla('ZCPOSV', -INFO );
         RETURN
      }

      // Quick return if (N == 0).

      if (N == 0) RETURN;

      // Skip single precision iterative refinement if a priori slower
      // than double precision factorization.

      if ( .NOT.DOITREF ) {
         ITER = -1
         GO TO 40
      }

      // Compute some constants.

      ANRM = ZLANHE( 'I', UPLO, N, A, LDA, RWORK )
      EPS = DLAMCH( 'Epsilon' )
      CTE = ANRM*EPS*SQRT( DBLE( N ) )*BWDMAX

      // Set the indices PTSA, PTSX for referencing SA and SX in SWORK.

      PTSA = 1
      PTSX = PTSA + N*N

      // Convert B from double precision to single precision and store the
      // result in SX.

      zlag2c(N, NRHS, B, LDB, SWORK( PTSX ), N, INFO );

      if ( INFO != 0 ) {
         ITER = -2
         GO TO 40
      }

      // Convert A from double precision to single precision and store the
      // result in SA.

      zlat2c(UPLO, N, A, LDA, SWORK( PTSA ), N, INFO );

      if ( INFO != 0 ) {
         ITER = -2
         GO TO 40
      }

      // Compute the Cholesky factorization of SA.

      cpotrf(UPLO, N, SWORK( PTSA ), N, INFO );

      if ( INFO != 0 ) {
         ITER = -3
         GO TO 40
      }

      // Solve the system SA*SX = SB.

      cpotrs(UPLO, N, NRHS, SWORK( PTSA ), N, SWORK( PTSX ), N, INFO );

      // Convert SX back to COMPLEX*16

      clag2z(N, NRHS, SWORK( PTSX ), N, X, LDX, INFO );

      // Compute R = B - AX (R is WORK).

      zlacpy('All', N, NRHS, B, LDB, WORK, N );

      zhemm('Left', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, ONE, WORK, N );

      // Check whether the NRHS normwise backward errors satisfy the
      // stopping criterion. If yes, set ITER=0 and return.

      for (I = 1; I <= NRHS; I++) {
         XNRM = CABS1( X( IZAMAX( N, X( 1, I ), 1 ), I ) )
         RNRM = CABS1( WORK( IZAMAX( N, WORK( 1, I ), 1 ), I ) )
         if (RNRM > XNRM*CTE) GO TO 10;
      }

      // If we are here, the NRHS normwise backward errors satisfy the
      // stopping criterion. We are good to exit.

      ITER = 0
      RETURN

      } // 10

      for (IITER = 1; IITER <= ITERMAX; IITER++) { // 30

         // Convert R (in WORK) from double precision to single precision
         // and store the result in SX.

         zlag2c(N, NRHS, WORK, N, SWORK( PTSX ), N, INFO );

         if ( INFO != 0 ) {
            ITER = -2
            GO TO 40
         }

         // Solve the system SA*SX = SR.

         cpotrs(UPLO, N, NRHS, SWORK( PTSA ), N, SWORK( PTSX ), N, INFO );

         // Convert SX back to double precision and update the current
         // iterate.

         clag2z(N, NRHS, SWORK( PTSX ), N, WORK, N, INFO );

         for (I = 1; I <= NRHS; I++) {
            zaxpy(N, ONE, WORK( 1, I ), 1, X( 1, I ), 1 );
         }

         // Compute R = B - AX (R is WORK).

         zlacpy('All', N, NRHS, B, LDB, WORK, N );

         zhemm('L', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, ONE, WORK, N );

         // Check whether the NRHS normwise backward errors satisfy the
         // stopping criterion. If yes, set ITER=IITER>0 and return.

         for (I = 1; I <= NRHS; I++) {
            XNRM = CABS1( X( IZAMAX( N, X( 1, I ), 1 ), I ) )
            RNRM = CABS1( WORK( IZAMAX( N, WORK( 1, I ), 1 ), I ) )
            if (RNRM > XNRM*CTE) GO TO 20;
         }

         // If we are here, the NRHS normwise backward errors satisfy the
         // stopping criterion, we are good to exit.

         ITER = IITER

         RETURN

         } // 20

      } // 30

      // If we are at this place of the code, this is because we have
      // performed ITER=ITERMAX iterations and never satisfied the
      // stopping criterion, set up the ITER flag accordingly and follow
      // up on double precision routine.

      ITER = -ITERMAX - 1

      } // 40

      // Single-precision iterative refinement failed to converge to a
      // satisfactory solution, so we resort to double precision.

      zpotrf(UPLO, N, A, LDA, INFO );

      if (INFO != 0) RETURN;

      zlacpy('All', N, NRHS, B, LDB, X, LDX );
      zpotrs(UPLO, N, NRHS, A, LDA, X, LDX, INFO );

      RETURN

      // End of ZCPOSV

      }
