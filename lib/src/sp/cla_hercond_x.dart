      double cla_hercond_x(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> AF_, final int LDAF, final Array<int> IPIV_, final int X, final int INFO, final Array<double> _WORK_, final Array<double> RWORK_,) {
  final A = A_.dim();
  final AF = AF_.dim();
  final IPIV = IPIV_.dim();
  final _WORK = _WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                N, LDA, LDAF, INFO;
      int                IPIV( * );
      Complex            A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * );
      double               RWORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                KASE, I, J;
      double               AINVNM, ANORM, TMP;
      bool               UP, UPPER;
      Complex            ZDUM;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CHETRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Statement Functions ..
      double CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1[ZDUM] = ( double( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();

      CLA_HERCOND_X = 0.0;

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( LDAF < max( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('CLA_HERCOND_X', -INFO );
         return;
      }
      UP = false;
      if ( lsame( UPLO, 'U' ) ) UP = true;

      // Compute norm of op(A)*op2(C).

      ANORM = 0.0;
      if ( UP ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            for (J = 1; J <= I; J++) {
               TMP = TMP + CABS1( A( J, I ) * X( J ) );
            }
            for (J = I+1; J <= N; J++) {
               TMP = TMP + CABS1( A( I, J ) * X( J ) );
            }
            RWORK[I] = TMP;
            ANORM = max( ANORM, TMP );
         }
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            for (J = 1; J <= I; J++) {
               TMP = TMP + CABS1( A( I, J ) * X( J ) );
            }
            for (J = I+1; J <= N; J++) {
               TMP = TMP + CABS1( A( J, I ) * X( J ) );
            }
            RWORK[I] = TMP;
            ANORM = max( ANORM, TMP );
         }
      }

      // Quick return if possible.

      if ( N == 0 ) {
         CLA_HERCOND_X = 1.0;
         return;
      } else if ( ANORM == 0.0 ) {
         return;
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0;

      KASE = 0;
      } // 10
      clacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == 2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK[I] = WORK( I ) * RWORK( I );
            }

            if ( UP ) {
               chetrs('U', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               chetrs('L', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by inv(X).

            for (I = 1; I <= N; I++) {
               WORK[I] = WORK( I ) / X( I );
            }
         } else {

            // Multiply by inv(X**H).

            for (I = 1; I <= N; I++) {
               WORK[I] = WORK( I ) / X( I );
            }

            if ( UP ) {
               chetrs('U', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               chetrs('L', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK[I] = WORK( I ) * RWORK( I );
            }
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != 0.0) CLA_HERCOND_X = 1.0 / AINVNM;

      }
