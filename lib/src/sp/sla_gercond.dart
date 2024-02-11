      double sla_gercond(final int TRANS, final int N, final Matrix<double> A, final int LDA, final Matrix<double> AF, final int LDAF, final Array<int> IPIV, final int CMODE, final int C, final int INFO, final Array<double> _WORK, final Array<int> IWORK,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                N, LDA, LDAF, INFO, CMODE;
      int                IPIV( * ), IWORK( * );
      double               A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * );
// ..

// =====================================================================

      // .. Local Scalars ..
      bool               NOTRANS;
      int                KASE, I, J;
      double               AINVNM, TMP;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SGETRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX

      SLA_GERCOND = 0.0;

      INFO = 0;
      NOTRANS = lsame( TRANS, 'N' );
      if ( !NOTRANS && !lsame(TRANS, 'T') && !lsame(TRANS, 'C') ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( LDAF < max( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('SLA_GERCOND', -INFO );
         return;
      }
      if ( N == 0 ) {
         SLA_GERCOND = 1.0;
         return;
      }

      // Compute the equilibration matrix R such that
      // inv(R)*A*C has unit 1-norm.

      if (NOTRANS) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            if ( CMODE == 1 ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) * C( J ) );
               }
            } else if ( CMODE == 0 ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ( A( I, J ) ).abs();
               }
            } else {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) / C( J ) );
               }
            }
            WORK[2*N+I] = TMP;
         }
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            if ( CMODE == 1 ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) * C( J ) );
               }
            } else if ( CMODE == 0 ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ( A( J, I ) ).abs();
               }
            } else {
               for (J = 1; J <= N; J++) {
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
               WORK[I] = WORK(I) * WORK(2*N+I);
            }

            if (NOTRANS) {
               sgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               sgetrs('Transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
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

            if (NOTRANS) {
               sgetrs('Transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               sgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK[I] = WORK( I ) * WORK( 2*N+I );
            }
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != 0.0) SLA_GERCOND = ( 1.0 / AINVNM );

      }
