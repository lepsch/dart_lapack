      double sla_porcond(final int UPLO, final int N, final Matrix<double> A, final int LDA, final Matrix<double> AF, final int LDAF, final int CMODE, final int C, final int INFO, final Array<double> _WORK, final Array<int> IWORK) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                N, LDA, LDAF, INFO, CMODE;
      double               A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * );
      int                IWORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                KASE, I, J;
      double               AINVNM, TMP;
      bool               UP;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SPOTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX

      SLA_PORCOND = 0.0;

      INFO = 0;
      if ( N < 0 ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('SLA_PORCOND', -INFO );
         return;
      }

      if ( N == 0 ) {
         SLA_PORCOND = 1.0;
         return;
      }
      UP = false;
      if ( lsame( UPLO, 'U' ) ) UP = true;

      // Compute the equilibration matrix R such that
      // inv(R)*A*C has unit 1-norm.

      if ( UP ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            if ( CMODE == 1 ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( J, I ) * C( J ) );
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) * C( J ) );
               }
            } else if ( CMODE == 0 ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ( A( J, I ) ).abs();
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ( A( I, J ) ).abs();
               }
            } else {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( J ,I ) / C( J ) );
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) / C( J ) );
               }
            }
            WORK[2*N+I] = TMP;
         }
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            if ( CMODE == 1 ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( I, J ) * C( J ) );
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) * C( J ) );
               }
            } else if ( CMODE == 0 ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ( A( I, J ) ).abs();
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ( A( J, I ) ).abs();
               }
            } else {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( I, J ) / C( J ) );
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) / C( J ) );
               }
            }
            WORK[2*N+I] = TMP;
         }
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0;

      KASE = 0;
      } // 10
      slacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == 2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK[I] = WORK( I ) * WORK( 2*N+I );
            }

            if (UP) {
               spotrs('Upper', N, 1, AF, LDAF, WORK, N, INFO );
            } else {
               spotrs('Lower', N, 1, AF, LDAF, WORK, N, INFO );
            }

            // Multiply by inv(C).

            if ( CMODE == 1 ) {
               for (I = 1; I <= N; I++) {
                  WORK[I] = WORK( I ) / C( I );
               }
            } else if ( CMODE == -1 ) {
               for (I = 1; I <= N; I++) {
                  WORK[I] = WORK( I ) * C( I );
               }
            }
         } else {

            // Multiply by inv(C**T).

            if ( CMODE == 1 ) {
               for (I = 1; I <= N; I++) {
                  WORK[I] = WORK( I ) / C( I );
               }
            } else if ( CMODE == -1 ) {
               for (I = 1; I <= N; I++) {
                  WORK[I] = WORK( I ) * C( I );
               }
            }

            if ( UP ) {
               spotrs('Upper', N, 1, AF, LDAF, WORK, N, INFO );
            } else {
               spotrs('Lower', N, 1, AF, LDAF, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK[I] = WORK( I ) * WORK( 2*N+I );
            }
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != 0.0) SLA_PORCOND = ( 1.0 / AINVNM );

      }
