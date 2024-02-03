      SUBROUTINE ZTRT03( UPLO, TRANS, DIAG, N, NRHS, A, LDA, SCALE, CNORM, TSCAL, X, LDX, B, LDB, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      double             RESID, SCALE, TSCAL;
      // ..
      // .. Array Arguments ..
      double             CNORM( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
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
      // EXTERNAL ZAXPY, ZCOPY, ZDSCAL, ZTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0

      if ( N.LE.0 || NRHS.LE.0 ) {
         RESID = ZERO
         RETURN
      }
      EPS = DLAMCH( 'Epsilon' )
      SMLNUM = DLAMCH( 'Safe minimum' )

      // Compute the norm of the triangular matrix A using the column
      // norms already computed by ZLATRS.

      TNORM = ZERO
      if ( LSAME( DIAG, 'N' ) ) {
         for (J = 1; J <= N; J++) { // 10
            TNORM = MAX( TNORM, TSCAL*ABS( A( J, J ) )+CNORM( J ) )
         } // 10
      } else {
         for (J = 1; J <= N; J++) { // 20
            TNORM = MAX( TNORM, TSCAL+CNORM( J ) )
         } // 20
      }

      // Compute the maximum over the number of right hand sides of
         // norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 30
         zcopy(N, X( 1, J ), 1, WORK, 1 );
         IX = IZAMAX( N, WORK, 1 )
         XNORM = MAX( ONE, ABS( X( IX, J ) ) )
         XSCAL = ( ONE / XNORM ) / DBLE( N )
         zdscal(N, XSCAL, WORK, 1 );
         ztrmv(UPLO, TRANS, DIAG, N, A, LDA, WORK, 1 );
         zaxpy(N, DCMPLX( -SCALE*XSCAL ), B( 1, J ), 1, WORK, 1 );
         IX = IZAMAX( N, WORK, 1 )
         ERR = TSCAL*ABS( WORK( IX ) )
         IX = IZAMAX( N, X( 1, J ), 1 )
         XNORM = ABS( X( IX, J ) )
         if ( ERR*SMLNUM.LE.XNORM ) {
            if (XNORM.GT.ZERO) ERR = ERR / XNORM;
         } else {
            if (ERR.GT.ZERO) ERR = ONE / EPS;
         }
         if ( ERR*SMLNUM.LE.TNORM ) {
            if (TNORM.GT.ZERO) ERR = ERR / TNORM;
         } else {
            if (ERR.GT.ZERO) ERR = ONE / EPS;
         }
         RESID = MAX( RESID, ERR )
      } // 30

      RETURN

      // End of ZTRT03

      }
