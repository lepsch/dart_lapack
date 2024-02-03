      double zla_gbrcond_x(TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, X, INFO, WORK, RWORK ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                N, KL, KU, KD, KE, LDAB, LDAFB, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ), X( * );
      double             RWORK( * );


// =====================================================================

      // .. Local Scalars ..
      bool               NOTRANS;
      int                KASE, I, J;
      double             AINVNM, ANORM, TMP;
      COMPLEX*16         ZDUM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLACN2, ZGBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) );
      // ..
      // .. Executable Statements ..

      ZLA_GBRCOND_X = 0.0;

      INFO = 0;
      NOTRANS = LSAME( TRANS, 'N' );
      if ( !NOTRANS && !LSAME(TRANS, 'T') && !LSAME( TRANS, 'C' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KL < 0 || KL > N-1 ) {
         INFO = -3;
      } else if ( KU < 0 || KU > N-1 ) {
         INFO = -4;
      } else if ( LDAB < KL+KU+1 ) {
         INFO = -6;
      } else if ( LDAFB < 2*KL+KU+1 ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('ZLA_GBRCOND_X', -INFO );
         return;
      }

      // Compute norm of op(A)*op2(C).

      KD = KU + 1;
      KE = KL + 1;
      ANORM = 0.0;
      if ( NOTRANS ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            DO J = max( I-KL, 1 ), min( I+KU, N );
               TMP = TMP + CABS1( AB( KD+I-J, J) * X( J ) );
            }
            RWORK( I ) = TMP;
            ANORM = max( ANORM, TMP );
         }
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            DO J = max( I-KL, 1 ), min( I+KU, N );
               TMP = TMP + CABS1( AB( KE-I+J, I ) * X( J ) );
            }
            RWORK( I ) = TMP;
            ANORM = max( ANORM, TMP );
         }
      }

      // Quick return if possible.

      if ( N == 0 ) {
         ZLA_GBRCOND_X = 1.0;
         return;
      } else if ( ANORM == 0.0 ) {
         return;
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0;

      KASE = 0;
      } // 10
      zlacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == 2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * RWORK( I );
            }

            if ( NOTRANS ) {
               zgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            } else {
               zgbtrs('Conjugate transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            }

            // Multiply by inv(X).

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) / X( I );
            }
         } else {

            // Multiply by inv(X**H).

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) / X( I );
            }

            if ( NOTRANS ) {
               zgbtrs('Conjugate transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            } else {
               zgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * RWORK( I );
            }
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != 0.0) ZLA_GBRCOND_X = 1.0 / AINVNM;

      return;

      // End of ZLA_GBRCOND_X

      }
