      SUBROUTINE STPT03( UPLO, TRANS, DIAG, N, NRHS, AP, SCALE, CNORM, TSCAL, X, LDX, B, LDB, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, N, NRHS;
      REAL               RESID, SCALE, TSCAL
      // ..
      // .. Array Arguments ..
      REAL               AP( * ), B( LDB, * ), CNORM( * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                IX, J, JJ;
      REAL               BIGNUM, EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SLAMCH
      // EXTERNAL LSAME, ISAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SSCAL, STPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 .OR. NRHS.LE.0 ) {
         RESID = ZERO
         RETURN
      }
      EPS = SLAMCH( 'Epsilon' )
      SMLNUM = SLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM

      // Compute the norm of the triangular matrix A using the column
      // norms already computed by SLATPS.

      TNORM = ZERO
      if ( LSAME( DIAG, 'N' ) ) {
         if ( LSAME( UPLO, 'U' ) ) {
            JJ = 1
            for (J = 1; J <= N; J++) { // 10
               TNORM = MAX( TNORM, TSCAL*ABS( AP( JJ ) )+CNORM( J ) )
               JJ = JJ + J + 1
            } // 10
         } else {
            JJ = 1
            for (J = 1; J <= N; J++) { // 20
               TNORM = MAX( TNORM, TSCAL*ABS( AP( JJ ) )+CNORM( J ) )
               JJ = JJ + N - J + 1
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
         scopy(N, X( 1, J ), 1, WORK, 1 );
         IX = ISAMAX( N, WORK, 1 )
         XNORM = MAX( ONE, ABS( X( IX, J ) ) )
         XSCAL = ( ONE / XNORM ) / REAL( N )
         sscal(N, XSCAL, WORK, 1 );
         stpmv(UPLO, TRANS, DIAG, N, AP, WORK, 1 );
         saxpy(N, -SCALE*XSCAL, B( 1, J ), 1, WORK, 1 );
         IX = ISAMAX( N, WORK, 1 )
         ERR = TSCAL*ABS( WORK( IX ) )
         IX = ISAMAX( N, X( 1, J ), 1 )
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

      // End of STPT03

      }
