      void dtpt03(UPLO, TRANS, DIAG, N, NRHS, AP, SCALE, CNORM, TSCAL, X, LDX, B, LDB, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, N, NRHS;
      double             RESID, SCALE, TSCAL;
      // ..
      // .. Array Arguments ..
      double             AP( * ), B( LDB, * ), CNORM( * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                IX, J, JJ;
      double             BIGNUM, EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DSCAL, DTPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }
      EPS = DLAMCH( 'Epsilon' );
      SMLNUM = DLAMCH( 'Safe minimum' );
      BIGNUM = ONE / SMLNUM;

      // Compute the norm of the triangular matrix A using the column
      // norms already computed by DLATPS.

      TNORM = ZERO;
      if ( LSAME( DIAG, 'N' ) ) {
         if ( LSAME( UPLO, 'U' ) ) {
            JJ = 1;
            for (J = 1; J <= N; J++) { // 10
               TNORM = max( TNORM, TSCAL*ABS( AP( JJ ) )+CNORM( J ) );
               JJ = JJ + J + 1;
            } // 10
         } else {
            JJ = 1;
            for (J = 1; J <= N; J++) { // 20
               TNORM = max( TNORM, TSCAL*ABS( AP( JJ ) )+CNORM( J ) );
               JJ = JJ + N - J + 1;
            } // 20
         }
      } else {
         for (J = 1; J <= N; J++) { // 30
            TNORM = max( TNORM, TSCAL+CNORM( J ) );
         } // 30
      }

      // Compute the maximum over the number of right hand sides of
         // norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 40
         dcopy(N, X( 1, J ), 1, WORK, 1 );
         IX = IDAMAX( N, WORK, 1 );
         XNORM = max( ONE, ABS( X( IX, J ) ) );
         XSCAL = ( ONE / XNORM ) / DBLE( N );
         dscal(N, XSCAL, WORK, 1 );
         dtpmv(UPLO, TRANS, DIAG, N, AP, WORK, 1 );
         daxpy(N, -SCALE*XSCAL, B( 1, J ), 1, WORK, 1 );
         IX = IDAMAX( N, WORK, 1 );
         ERR = TSCAL*ABS( WORK( IX ) );
         IX = IDAMAX( N, X( 1, J ), 1 );
         XNORM = ABS( X( IX, J ) );
         if ( ERR*SMLNUM <= XNORM ) {
            if (XNORM > ZERO) ERR = ERR / XNORM;
         } else {
            if (ERR > ZERO) ERR = ONE / EPS;
         }
         if ( ERR*SMLNUM <= TNORM ) {
            if (TNORM > ZERO) ERR = ERR / TNORM;
         } else {
            if (ERR > ZERO) ERR = ONE / EPS;
         }
         RESID = max( RESID, ERR );
      } // 40

      return;
      }
