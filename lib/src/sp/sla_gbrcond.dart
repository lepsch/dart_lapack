      double sla_gbrcond(TRANS, N, KL, KU, final Matrix<double> AB, final int LDAB, final Matrix<double> AFB, final int LDAFB, final Array<int> IPIV, CMODE, C, INFO, WORK, IWORK ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                N, LDAB, LDAFB, INFO, KL, KU, CMODE;
      int                IWORK( * ), IPIV( * );
      double               AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ), C( * );
// ..

// =====================================================================

      // .. Local Scalars ..
      bool               NOTRANS;
      int                KASE, I, J, KD, KE;
      double               AINVNM, TMP;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SGBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX

      SLA_GBRCOND = 0.0;

      INFO = 0;
      NOTRANS = lsame( TRANS, 'N' );
      if ( !NOTRANS && !lsame(TRANS, 'T') && !lsame(TRANS, 'C') ) {
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
         xerbla('SLA_GBRCOND', -INFO );
         return;
      }
      if ( N == 0 ) {
         SLA_GBRCOND = 1.0;
         return;
      }

      // Compute the equilibration matrix R such that
      // inv(R)*A*C has unit 1-norm.

      KD = KU + 1;
      KE = KL + 1;
      if ( NOTRANS ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
               if ( CMODE == 1 ) {
               for (J = max( I-KL, 1 ); J <= min( I+KU, N ); J++) {
                  TMP = TMP + ABS( AB( KD+I-J, J ) * C( J ) );
               }
               } else if ( CMODE == 0 ) {
                  for (J = max( I-KL, 1 ); J <= min( I+KU, N ); J++) {
                     TMP = TMP + ( AB( KD+I-J, J ) ).abs();
                  }
               } else {
                  for (J = max( I-KL, 1 ); J <= min( I+KU, N ); J++) {
                     TMP = TMP + ABS( AB( KD+I-J, J ) / C( J ) );
                  }
               }
            WORK[2*N+I] = TMP;
         }
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0;
            if ( CMODE == 1 ) {
               for (J = max( I-KL, 1 ); J <= min( I+KU, N ); J++) {
                  TMP = TMP + ABS( AB( KE-I+J, I ) * C( J ) );
               }
            } else if ( CMODE == 0 ) {
               for (J = max( I-KL, 1 ); J <= min( I+KU, N ); J++) {
                  TMP = TMP + ( AB( KE-I+J, I ) ).abs();
               }
            } else {
               for (J = max( I-KL, 1 ); J <= min( I+KU, N ); J++) {
                  TMP = TMP + ABS( AB( KE-I+J, I ) / C( J ) );
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

            if ( NOTRANS ) {
               sgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            } else {
               sgbtrs('Transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
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

            if ( NOTRANS ) {
               sgbtrs('Transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            } else {
               sgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK[I] = WORK( I ) * WORK( 2*N+I );
            }
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != 0.0) SLA_GBRCOND = ( 1.0 / AINVNM );

      }
