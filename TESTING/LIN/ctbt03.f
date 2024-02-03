      SUBROUTINE CTBT03( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, SCALE, CNORM, TSCAL, X, LDX, B, LDB, WORK, RESID );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                KD, LDAB, LDB, LDX, N, NRHS;
      REAL               RESID, SCALE, TSCAL;
      // ..
      // .. Array Arguments ..
      REAL               CNORM( * );
      COMPLEX            AB( LDAB, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================


      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                IX, J;
      REAL               EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX;
      REAL               SLAMCH;
      // EXTERNAL LSAME, ICAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CSSCAL, CTBMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX, REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }
      EPS = SLAMCH( 'Epsilon' );
      SMLNUM = SLAMCH( 'Safe minimum' );

      // Compute the norm of the triangular matrix A using the column
      // norms already computed by CLATBS.

      TNORM = ZERO;
      if ( LSAME( DIAG, 'N' ) ) {
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 10
               TNORM = MAX( TNORM, TSCAL*ABS( AB( KD+1, J ) )+ CNORM( J ) );
            } // 10
         } else {
            for (J = 1; J <= N; J++) { // 20
               TNORM = MAX( TNORM, TSCAL*ABS( AB( 1, J ) )+CNORM( J ) );
            } // 20
         }
      } else {
         for (J = 1; J <= N; J++) { // 30
            TNORM = MAX( TNORM, TSCAL+CNORM( J ) );
         } // 30
      }

      // Compute the maximum over the number of right hand sides of
         // norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 40
         ccopy(N, X( 1, J ), 1, WORK, 1 );
         IX = ICAMAX( N, WORK, 1 );
         XNORM = MAX( ONE, ABS( X( IX, J ) ) );
         XSCAL = ( ONE / XNORM ) / REAL( KD+1 );
         csscal(N, XSCAL, WORK, 1 );
         ctbmv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1 );
         caxpy(N, CMPLX( -SCALE*XSCAL ), B( 1, J ), 1, WORK, 1 );
         IX = ICAMAX( N, WORK, 1 );
         ERR = TSCAL*ABS( WORK( IX ) );
         IX = ICAMAX( N, X( 1, J ), 1 );
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
         RESID = MAX( RESID, ERR );
      } // 40

      return;

      // End of CTBT03

      }
