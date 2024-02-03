      REAL cla_gercond_c(TRANS, N, A, LDA, AF, LDAF, IPIV, C, CAPPLY, INFO, WORK, RWORK ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      bool               CAPPLY;
      int                N, LDA, LDAF, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * );
      REAL               C( * ), RWORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               NOTRANS;
      int                KASE, I, J;
      REAL               AINVNM, ANORM, TMP;
      COMPLEX            ZDUM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CGETRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL, AIMAG
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) );
      // ..
      // .. Executable Statements ..
      CLA_GERCOND_C = 0.0;

      INFO = 0;
      NOTRANS = LSAME( TRANS, 'N' );
      if ( !NOTRANS && !LSAME( TRANS, 'T' ) && !LSAME( TRANS, 'C' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( LDAF < max( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('CLA_GERCOND_C', -INFO );
         return;
      }

      // Compute norm of op(A)*op2(C).

      ANORM = 0.0;
      if ( NOTRANS ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            if ( CAPPLY ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + CABS1( A( I, J ) ) / C( J );
               }
            } else {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + CABS1( A( I, J ) );
               }
            }
            RWORK( I ) = TMP;
            ANORM = max( ANORM, TMP );
         }
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            if ( CAPPLY ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + CABS1( A( J, I ) ) / C( J );
               }
            } else {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + CABS1( A( J, I ) );
               }
            }
            RWORK( I ) = TMP;
            ANORM = max( ANORM, TMP );
         }
      }

      // Quick return if possible.

      if ( N == 0 ) {
         CLA_GERCOND_C = 1.0;
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
               WORK( I ) = WORK( I ) * RWORK( I );
            }

            if (NOTRANS) {
               cgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               cgetrs('Conjugate transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by inv(C).

            if ( CAPPLY ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I );
               }
            }
         } else {

            // Multiply by inv(C**H).

            if ( CAPPLY ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I );
               }
            }

            if ( NOTRANS ) {
               cgetrs('Conjugate transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               cgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * RWORK( I );
            }
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != 0.0) CLA_GERCOND_C = 1.0 / AINVNM;

      return;

      // End of CLA_GERCOND_C

      }
