      void zsyt01_3(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> AFAC_, final int LDAFAC, final int E, final Array<int> IPIV_, final Matrix<double> C_, final int LDC, final Array<double> RWORK_, final int RESID,) {
  final A = A_.dim();
  final AFAC = AFAC_.dim();
  final IPIV = IPIV_.dim();
  final C = C_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      double             RESID;
      int                IPIV( * );
      double             RWORK( * );
      Complex         A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), E( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I, INFO, J;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, ZLANSY;
      // EXTERNAL lsame, DLAMCH, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASET, ZLAVSY_ROOK, ZSYCONVF_ROOK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // a) Revert to multipliers of L

      zsyconvf_rook(UPLO, 'R', N, AFAC, LDAFAC, E, IPIV, INFO );

      // 1) Determine EPS and the norm of A.

      EPS = dlamch( 'Epsilon' );
      ANORM = ZLANSY( '1', UPLO, N, A, LDA, RWORK );

      // 2) Initialize C to the identity matrix.

      zlaset('Full', N, N, CZERO, CONE, C, LDC );

      // 3) Call ZLAVSY_ROOK to form the product D * U' (or D * L' ).

      zlavsy_rook(UPLO, 'Transpose', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // 4) Call ZLAVSY_ROOK again to multiply by U (or L ).

      zlavsy_rook(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // 5) Compute the difference  C - A .

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
               C[I][J] = C( I, J ) - A( I, J );
            }
         }
      } else {
         for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
               C[I][J] = C( I, J ) - A( I, J );
            }
         }
      }

      // 6) Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = ZLANSY( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / N.toDouble() ) / ANORM ) / EPS;
      }


      // b) Convert to factor of L (or U)

      zsyconvf_rook(UPLO, 'C', N, AFAC, LDAFAC, E, IPIV, INFO );

      }
