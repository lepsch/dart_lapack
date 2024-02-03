      double           FUNCTION DLA_PORCOND( UPLO, N, A, LDA, AF, LDAF, CMODE, C, INFO, WORK, IWORK );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, LDA, LDAF, INFO, CMODE;
      double             A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * );
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                KASE, I, J;
      double             AINVNM, TMP;
      bool               UP;
      // ..
      // .. Array Arguments ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DPOTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      DLA_PORCOND = 0.0;

      INFO = 0;
      if ( N < 0 ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('DLA_PORCOND', -INFO );
         return;
      }

      if ( N == 0 ) {
         DLA_PORCOND = 1.0;
         return;
      }
      UP = false;
      if ( LSAME( UPLO, 'U' ) ) UP = true;

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
                  TMP = TMP + ABS( A( J, I ) );
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) );
               }
            } else {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( J ,I ) / C( J ) );
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) / C( J ) );
               }
            }
            WORK( 2*N+I ) = TMP;
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
                  TMP = TMP + ABS( A( I, J ) );
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) );
               }
            } else {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + ABS( A( I, J ) / C( J ) );
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) / C( J ) );
               }
            }
            WORK( 2*N+I ) = TMP;
         }
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0;

      KASE = 0;
      } // 10
      dlacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == 2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * WORK( 2*N+I );
            }

            if (UP) {
               dpotrs('Upper', N, 1, AF, LDAF, WORK, N, INFO );
            } else {
               dpotrs('Lower', N, 1, AF, LDAF, WORK, N, INFO );
            }

            // Multiply by inv(C).

            if ( CMODE == 1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) / C( I );
               }
            } else if ( CMODE == -1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I );
               }
            }
         } else {

            // Multiply by inv(C**T).

            if ( CMODE == 1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) / C( I );
               }
            } else if ( CMODE == -1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I );
               }
            }

            if ( UP ) {
               dpotrs('Upper', N, 1, AF, LDAF, WORK, N, INFO );
            } else {
               dpotrs('Lower', N, 1, AF, LDAF, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * WORK( 2*N+I );
            }
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != 0.0) DLA_PORCOND = ( 1.0 / AINVNM );

      return;

      // End of DLA_PORCOND

      }
