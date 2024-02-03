      void ztbt03(UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, SCALE, CNORM, TSCAL, X, LDX, B, LDB, WORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                KD, LDAB, LDB, LDX, N, NRHS;
      double             RESID, SCALE, TSCAL;
      // ..
      // .. Array Arguments ..
      double             CNORM( * );
      Complex         AB( LDAB, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================


      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                IX, J;
      double             EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IZAMAX, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY, ZDSCAL, ZTBMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }
      EPS = DLAMCH( 'Epsilon' );
      SMLNUM = DLAMCH( 'Safe minimum' );

      // Compute the norm of the triangular matrix A using the column
      // norms already computed by ZLATBS.

      TNORM = ZERO;
      if ( LSAME( DIAG, 'N' ) ) {
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 10
               TNORM = max( TNORM, TSCAL*( AB( KD+1, J ) ).abs()+ CNORM( J ) );
            } // 10
         } else {
            for (J = 1; J <= N; J++) { // 20
               TNORM = max( TNORM, TSCAL*( AB( 1, J ) ).abs()+CNORM( J ) );
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
         zcopy(N, X( 1, J ), 1, WORK, 1 );
         IX = IZAMAX( N, WORK, 1 );
         XNORM = max( ONE, ( X( IX, J ) ) ).abs();
         XSCAL = ( ONE / XNORM ) / DBLE( KD+1 );
         zdscal(N, XSCAL, WORK, 1 );
         ztbmv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1 );
         zaxpy(N, DCMPLX( -SCALE*XSCAL ), B( 1, J ), 1, WORK, 1 );
         IX = IZAMAX( N, WORK, 1 );
         ERR = TSCAL*( WORK( IX ) ).abs();
         IX = IZAMAX( N, X( 1, J ), 1 );
         XNORM = ( X( IX, J ) ).abs();
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
