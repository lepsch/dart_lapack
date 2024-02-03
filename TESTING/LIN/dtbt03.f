      SUBROUTINE DTBT03( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, SCALE, CNORM, TSCAL, X, LDX, B, LDB, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                KD, LDAB, LDB, LDX, N, NRHS;
      double             RESID, SCALE, TSCAL;
      // ..
      // .. Array Arguments ..
      double             AB( LDAB, * ), B( LDB, * ), CNORM( * ), WORK( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                IX, J;
      double             BIGNUM, EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DSCAL, DTBMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0

      if ( N.LE.0 .OR. NRHS.LE.0 ) {
         RESID = ZERO
         RETURN
      }
      EPS = DLAMCH( 'Epsilon' )
      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM

      // Compute the norm of the triangular matrix A using the column
      // norms already computed by DLATBS.

      TNORM = ZERO
      if ( LSAME( DIAG, 'N' ) ) {
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 10
               TNORM = MAX( TNORM, TSCAL*ABS( AB( KD+1, J ) )+ CNORM( J ) )
            } // 10
         } else {
            for (J = 1; J <= N; J++) { // 20
               TNORM = MAX( TNORM, TSCAL*ABS( AB( 1, J ) )+CNORM( J ) )
            } // 20
         }
      } else {
         for (J = 1; J <= N; J++) { // 30
            TNORM = MAX( TNORM, TSCAL+CNORM( J ) )
         } // 30
      }

      // Compute the maximum over the number of right hand sides of
         // norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 40
         dcopy(N, X( 1, J ), 1, WORK, 1 );
         IX = IDAMAX( N, WORK, 1 )
         XNORM = MAX( ONE, ABS( X( IX, J ) ) )
         XSCAL = ( ONE / XNORM ) / DBLE( KD+1 )
         dscal(N, XSCAL, WORK, 1 );
         dtbmv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1 );
         daxpy(N, -SCALE*XSCAL, B( 1, J ), 1, WORK, 1 );
         IX = IDAMAX( N, WORK, 1 )
         ERR = TSCAL*ABS( WORK( IX ) )
         IX = IDAMAX( N, X( 1, J ), 1 )
         XNORM = ABS( X( IX, J ) )
         if ( ERR*SMLNUM.LE.XNORM ) {
            IF( XNORM.GT.ZERO ) ERR = ERR / XNORM
         } else {
            IF( ERR.GT.ZERO ) ERR = ONE / EPS
         }
         if ( ERR*SMLNUM.LE.TNORM ) {
            IF( TNORM.GT.ZERO ) ERR = ERR / TNORM
         } else {
            IF( ERR.GT.ZERO ) ERR = ONE / EPS
         }
         RESID = MAX( RESID, ERR )
      } // 40

      RETURN

      // End of DTBT03

      }
