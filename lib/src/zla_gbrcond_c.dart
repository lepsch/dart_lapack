      double zla_gbrcond_c(TRANS, N, KL, KU, final Matrix<double> AB, final int LDAB, final Matrix<double> AFB, final int LDAFB, final Array<int> IPIV, C, CAPPLY, INFO, final Array<double> _WORK, RWORK ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      bool               CAPPLY;
      int                N, KL, KU, KD, KE, LDAB, LDAFB, INFO;
      int                IPIV( * );
      Complex         AB( LDAB, * ), AFB( LDAFB, * ), WORK( * );
      double             C( * ), RWORK( * );


// =====================================================================

      // .. Local Scalars ..
      bool               NOTRANS;
      int                KASE, I, J;
      double             AINVNM, ANORM, TMP;
      Complex         ZDUM;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
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
      CABS1[ZDUM] = ( ZDUM.toDouble() ).abs() + ( DIMAG( ZDUM ) ).abs();
      ZLA_GBRCOND_C = 0.0;

      INFO = 0;
      NOTRANS = lsame( TRANS, 'N' );
      if ( !NOTRANS && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
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
         xerbla('ZLA_GBRCOND_C', -INFO );
         return;
      }

      // Compute norm of op(A)*op2(C).

      ANORM = 0.0;
      KD = KU + 1;
      KE = KL + 1;
      if ( NOTRANS ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            if ( CAPPLY ) {
               for (J = max( I-KL, 1 ); J <= min( I+KU, N ); J++) {
                  TMP = TMP + CABS1( AB( KD+I-J, J ) ) / C( J );
               }
            } else {
               for (J = max( I-KL, 1 ); J <= min( I+KU, N ); J++) {
                  TMP = TMP + CABS1( AB( KD+I-J, J ) );
               }
            }
            RWORK[I] = TMP;
            ANORM = max( ANORM, TMP );
         }
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            if ( CAPPLY ) {
               for (J = max( I-KL, 1 ); J <= min( I+KU, N ); J++) {
                  TMP = TMP + CABS1( AB( KE-I+J, I ) ) / C( J );
               }
            } else {
               for (J = max( I-KL, 1 ); J <= min( I+KU, N ); J++) {
                  TMP = TMP + CABS1( AB( KE-I+J, I ) );
               }
            }
            RWORK[I] = TMP;
            ANORM = max( ANORM, TMP );
         }
      }

      // Quick return if possible.

      if ( N == 0 ) {
         ZLA_GBRCOND_C = 1.0;
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
               WORK[I] = WORK( I ) * RWORK( I );
            }

            if ( NOTRANS ) {
               zgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            } else {
               zgbtrs('Conjugate transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            }

            // Multiply by inv(C).

            if ( CAPPLY ) {
               for (I = 1; I <= N; I++) {
                  WORK[I] = WORK( I ) * C( I );
               }
            }
         } else {

            // Multiply by inv(C**H).

            if ( CAPPLY ) {
               for (I = 1; I <= N; I++) {
                  WORK[I] = WORK( I ) * C( I );
               }
            }

            if ( NOTRANS ) {
               zgbtrs('Conjugate transpose', N, KL, KU, 1, AFB, LDAFB, IPIV,  WORK, N, INFO );
            } else {
               zgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK[I] = WORK( I ) * RWORK( I );
            }
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != 0.0) ZLA_GBRCOND_C = 1.0 / AINVNM;

      }
