      REAL cla_gbrcond_x(TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, X, INFO, WORK, RWORK ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                N, KL, KU, KD, KE, LDAB, LDAFB, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ), X( * );
      REAL               RWORK( * );
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
      //- bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CGBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ( REAL( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      CLA_GBRCOND_X = 0.0;

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
         xerbla('CLA_GBRCOND_X', -INFO );
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
         CLA_GBRCOND_X = 1.0;
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

            if ( NOTRANS ) {
               cgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            } else {
               cgbtrs('Conjugate transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
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
               cgbtrs('Conjugate transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            } else {
               cgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * RWORK( I );
            }
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != 0.0) CLA_GBRCOND_X = 1.0 / AINVNM;

      return;
      }
